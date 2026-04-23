#' @title Operating characteristics for the inserted Shared Keyboard design
#' 
#' @description 
#' Simulate phase I dose-finding trials under the inserted Shared Keyboard design (Insert-SKBD),
#' where the working dose grid can be adaptively refined by inserting new dose levels when
#' posterior evidence suggests that the target toxicity level is not well covered by the
#' current prespecified grid.
#'
#' @details
#' The prespecified dose set \code{dose_set} is internally standardized to \eqn{[0,1]} for
#' kernel construction and numerical stability. At each interim, Insert-SKBD:
#' \enumerate{
#'   \item updates Beta pseudo-posteriors across all current dose levels via kernel borrowing;
#'   \item checks whether dose insertion is warranted near the current dose based on
#'         posterior-evidence thresholds \code{C1} and \code{C2};
#'   \item if insertion occurs, augments the working grid and continues cohort-wise treating;
#'   \item applies overdose elimination and a keyboard-style move rule for escalation/de-escalation.
#' }
#'
#' The final selected dose is the admissible dose whose isotonic-smoothed posterior mean
#' toxicity is closest to \code{target_prob}. If the trial stops early for safety, no dose is selected.
#'
#' @param target_prob Scalar in \eqn{(0,1)}. The target toxicity rate \eqn{\phi}.
#' @param tox_prob Numeric vector in \eqn{[0,1]}. True DLT probabilities at prespecified doses.
#' Must have the same length as \code{dose_set}.
#' @param dose_set Numeric vector of prespecified doses on the original (clinical) scale.
#' Must be strictly increasing with no duplicates.
#' @param n_cohort Positive integer. Number of cohorts per simulated trial.
#' @param cohort_size Positive integer. Number of patients per cohort.
#' @param C1,C2 Scalars in \eqn{(0,1)}. Posterior-evidence thresholds for triggering insertion.
#' Larger values make insertion more conservative.
#' @param start_dose Positive integer. Starting dose index (in \code{1:length(dose_set)}).
#' @param margin_left,margin_right Nonnegative scalars. Define the target key
#' \eqn{(\phi-\epsilon_1,\phi+\epsilon_2)} with \code{key_L = target_prob - margin_left} and
#' \code{key_U = target_prob + margin_right}.
#' @param n_earlystop Positive integer. Maximum number of patients allowed at a single dose
#' before stopping accrual at that dose (operational cap).
#' @param cutoff_elimin Scalar in \eqn{(0,1)}. Overdose elimination cutoff: a dose is eliminated if
#' \eqn{Pr(\pi(d)>\phi \mid \mathcal D)} exceeds \code{cutoff_elimin} (or \code{cutoff_elimin - offset}
#' when \code{extra_safe = TRUE}).
#' @param extra_safe Logical. If \code{TRUE}, use a more conservative elimination cutoff
#' \code{cutoff_elimin - offset}.
#' @param offset Scalar in \eqn{[0,0.5)}. Safety offset used when \code{extra_safe = TRUE}.
#' @param symmetric Logical. If \code{TRUE}, use symmetric kernel borrowing; otherwise allow
#' asymmetric borrowing via \code{k_left} and \code{k_right}.
#' @param k_left,k_right Scalars in \eqn{(0,1)} controlling borrowing strength to the left and right
#' neighbors in the SKBD pseudo-posterior update. Typically \code{k_right > k_left} for safety.
#' @param ref_gap Optional positive scalar. Reference gap used for kernel scaling. Note: in the current
#' implementation, \code{ref_gap} is recomputed internally as \code{min(diff(dose_set_work))} after each
#' insertion/update.
#' @param shared Logical. Reserved for compatibility with the shared keyboard framework.
#' @param light_return Logical. If \code{TRUE}, return only summaries; if \code{FALSE}, also return
#' per-trial details in \code{trial_detail}.
#' @param n_trial Positive integer. Number of simulated trials.
#' @param seed Integer. RNG seed for reproducibility.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{sel_pct_prespec}}{Numeric vector. Selection percentage (%) at each prespecified dose.}
#'   \item{\code{pts_pct_prespec}}{Numeric vector. Patient allocation percentage (%) at each prespecified dose.}
#'   \item{\code{insertion}}{A list summarizing inserted-dose behavior:
#'     \code{sel_pct}, \code{pts_pct}, \code{dose_mean}, \code{dose_sd}, \code{trial_pct},
#'     \code{cohort_mean}, \code{n_median}.}
#'   \item{\code{simdata}}{A data frame storing per-trial, per-dose counts:
#'     \code{Simulation}, \code{Dose}, \code{N} (treated), \code{X} (DLTs), and \code{Selection}.}
#'   \item{\code{sel_dose_idx}}{Integer vector of selected dose indices on the *working* grid
#'     for each trial; \code{-1} indicates early stop / no selection.}
#'   \item{\code{sel_dose}}{Numeric vector of selected dose values on the original scale (NA if none).}
#'   \item{\code{n_insertions}}{Integer vector. Number of insertions performed in each trial.}
#'   \item{\code{insert_at_cohort}}{Integer vector. Cohort index when the first insertion occurs
#'     (\code{0} if no insertion).}
#'   \item{\code{trial_detail}}{(Only if \code{light_return = FALSE}) A list of length \code{n_trial}
#'     with per-trial working grids, counts, elimination flags, insertion info, and early-stop indicator.}
#' }
#'
#' @examples
#' tox_prob <- c(0.14, 0.45, 0.63, 0.74, 0.80)
#' dose_set <- c(5, 15, 25, 35, 45)
#' out <- get_OC_Insert_SKBD(
#'   target_prob = 0.30,
#'   tox_prob = tox_prob,
#'   dose_set = dose_set,
#'   n_trial = 100
#' )
#' out$insertion$trial_pct
#'
#' @importFrom stats median pbeta rbeta rbinom rexp runif sd
#' @importFrom utils tail
#' @export

get_OC_Insert_SKBD <- function(
    target_prob, tox_prob, dose_set,
    n_cohort = 10, cohort_size = 3,
    C1 = 0.6, C2 = 0.6, 
    start_dose = 1,
    margin_left = 0.05, margin_right = 0.05, 
    n_earlystop = 1000, 
    cutoff_elimin = 0.95,
    extra_safe = FALSE, offset = 0.05, 
    symmetric = FALSE, 
    k_left = 0.2, k_right = 0.8, ref_gap = NULL,
    shared = TRUE, # incorporate the keybooard design as a special case
    light_return = TRUE,
    n_trial = 1000, seed = 6
) {
  set.seed(seed)
  
  #----------------------------- argument checks ---------------------------------------------
  # target
  if (!is.numeric(target_prob) || length(target_prob) != 1 || is.na(target_prob) ||
      target_prob <= 0 || target_prob >= 1) {
    stop("`target_prob` must be a single number in (0, 1).")
  }
  if (target_prob < 0.05) warning("`target_prob` is very low (< 0.05); OC may be unstable.")
  if (target_prob > 0.5)  warning("`target_prob` is very high (> 0.5); uncommon in phase I.")
  
  # margins
  if (!is.numeric(margin_left) || !is.numeric(margin_right) ||
      length(margin_left) != 1 || length(margin_right) != 1 ||
      anyNA(c(margin_left, margin_right)) || margin_left < 0 || margin_right < 0) {
    stop("`margin_left` and `margin_right` must be nonnegative scalars.")
  }
  key_L = target_prob - margin_left
  key_U = target_prob + margin_right
  if (key_L <= 0 || key_U >= 1 || key_L >= key_U) {
    stop("Invalid target key: require 0 < target_prob - margin_left < target_prob + margin_right < 1.")
  }
  
  # tox_prob
  if (!is.numeric(tox_prob) || anyNA(tox_prob) || any(tox_prob < 0 | tox_prob > 1)) {
    stop("`tox_prob` must be numeric in [0, 1] with no NA.")
  }
  
  # n_cohort / cohort_size / n_trial
  int1 <- function(x) is.numeric(x) && length(x) == 1 && !is.na(x) && x %% 1 == 0 && x >= 1
  if (!int1(n_cohort))    stop("`n_cohort` must be a positive integer.")
  if (!int1(cohort_size)) stop("`cohort_size` must be a positive integer.")
  if (!int1(n_trial))     stop("`n_trial` must be a positive integer.")
  
  # n_earlystop (here used as posterior-sample size option placeholder; we keep the check)
  if (!int1(n_earlystop)) stop("`n_earlystop` must be a positive integer.")
  if (n_earlystop <= 6) {
    warning("`n_earlystop` is small; Recommended range is 9-18.")
  }
  
  # cutoff / extra safety
  if (!is.numeric(cutoff_elimin) || length(cutoff_elimin) != 1 || is.na(cutoff_elimin) ||
      cutoff_elimin <= 0 || cutoff_elimin >= 1) {
    stop("`cutoff_elimin` must be a scalar in (0, 1).")
  }
  if (!is.numeric(offset) || length(offset) != 1 || is.na(offset) || offset < 0 || offset >= 0.5) {
    stop("`offset` must be a scalar in [0, 0.5).")
  }
  if (extra_safe && (cutoff_elimin - offset) <= 0) {
    stop("When `extra_safe=TRUE`, require `cutoff_elimin - offset > 0`.")
  }
  
  # dose_set vs tox_prob
  if (!is.numeric(dose_set) || anyNA(dose_set)) stop("`dose_set` must be numeric with no NA.")
  if (length(dose_set) != length(tox_prob)) {
    stop("`dose_set` must have the same length as `tox_prob` (i.e., n_dose).")
  }
  if (anyDuplicated(dose_set)) stop("`dose_set` must have unique values (no duplicates).")
  if (is.unsorted(dose_set, strictly = TRUE)) {
    stop("`dose_set` must be strictly increasing (sorted).")
  }
  
  # start_dose index (in prespecified doses)
  if (!int1(start_dose) || start_dose < 1 || start_dose > length(tox_prob)) {
    stop("`start_dose` must be a single integer in 1:n_dose (indexing prespecified doses).")
  }
  
  # kernel borrowing strengths
  if (!is.numeric(k_left) || length(k_left) != 1 || is.na(k_left) || k_left <= 0 || k_left >= 1) {
    stop("`k_left` must be a single number in (0, 1).")
  }
  if (!is.numeric(k_right) || length(k_right) != 1 || is.na(k_right) || k_right <= 0 || k_right >= 1) {
    stop("`k_right` must be a single number in (0, 1).")
  }
  if (!is.null(ref_gap)) {
    if (!is.numeric(ref_gap) || length(ref_gap) != 1 || is.na(ref_gap) || !is.finite(ref_gap) || ref_gap <= 0) {
      stop("`ref_gap` must be a positive finite scalar when provided.")
    }
  }
  
  # Standardize doses to [0,1]
  dose_range = range(dose_set)
  dose_set_std = (dose_set - dose_range[1]) / (dose_range[2] - dose_range[1])
  
  #------------------------------- Simulation containers ------------------------------------------
  n_dose     = length(tox_prob)       # number of dose levels
  n_patients = n_cohort * cohort_size # total number of patients
  
  sel_dose_idx     = integer(n_trial) # selected dose level index (j*)
  sel_dose         = numeric(n_trial) # selected dose value (d_j*)
  
  insert_at_cohort = integer(n_trial) # cohort index when first insertion occurs (0 if none)
  n_insertions     = integer(n_trial) # number of insertions in each trial
  
  simdata = NULL
  
  # hyperparameters for beta prior
  pri_alpha = pri_beta = 1 
  
  # keep details optionally
  trial_detail = if (!light_return) vector("list", n_trial) else NULL
  
  ## get keys
  keys  = get_Key(target_prob, margin_left, margin_right)
  target_key = which.min(abs(keys - key_L))
  
  #--------------------------------- Main simulation loop --------------------------------------------
  for(trial in 1:n_trial){
    
    ## current dose index (starting dose for this trial)
    d = start_dose
    
    ## working dose grid (may expand after insertion)
    dose_set_work = dose_set_std   # current dose set; starts with prespecified doses and may add inserted doses
    tox_prob_work = tox_prob       # true toxicity probs aligned with dose_set_work; expanded when a new dose is inserted
    n_dose_work   = length(tox_prob_work)  # current number of dose levels on the working grid
    
    ## trial-level sufficient statistics on the working grid
    n_dlt_work     = rep(0, n_dose_work)   # number of DLTs observed at each dose level
    n_treated_work = rep(0, n_dose_work)   # number of patients treated at each dose level
    
    ## safety tracking
    is_eliminated = rep(FALSE, n_dose_work)  # dose elimination indicator (TRUE = eliminated / not allowed)
    is_earlystop  = FALSE                    # early termination flag (TRUE = stop the trial for safety)
    
    ## insertion bookkeeping
    n_insertions_trial = 0  # number of dose insertions performed in this trial
    insert_at_cohort[trial] = 0
    
    ## treat sequentially by cohort until total N or early stop
    for (cohort in 1:n_cohort){
      
      # if enough patients have been treated at dose d, exit the loop
      if(n_treated_work[d] >= n_earlystop){
        break
      }
      
      # Treat one cohort at the current dose level d
      n_dlt_work[d] = n_dlt_work[d] + rbinom(1, cohort_size, tox_prob_work[d])
      n_treated_work[d] = n_treated_work[d] + cohort_size
      
      ## BKP posterior update for ALL doses
      ref_gap = min(diff(dose_set_work))
      post_work_insert = post_par_all(
        n_dlt = n_dlt_work, 
        n_treated = n_treated_work, 
        dose_set = dose_set_work, 
        pri_alpha = 0.5, 
        pri_beta = 0.5, 
        symmetric = TRUE, 
        k_left = 0.01,
        k_right = 0.01,
        ref_gap = ref_gap
      )
      
      post_alpha_work_insert = post_work_insert$post_alpha
      post_beta_work_insert  = post_work_insert$post_beta
      
      # Determine if insertion is warranted
      res = insert_check(
        j = d,
        post_alpha = post_alpha_work_insert,
        post_beta  = post_beta_work_insert,
        key_L = key_L,
        key_U = key_U,
        C1 = C1,
        C2 = C2,
        n_treated = n_treated_work
      )
      
      need_insert = res$need_insert
      insert_code = res$insert_code
      LOC_BELOW_MIN = res$LOC_BELOW_MIN
      LOC_ABOVE_MAX = res$LOC_ABOVE_MAX
      
      if(need_insert){
        n_insertions_trial = n_insertions_trial + 1
      }
      
      #-------------- insertion block ----------------------------- 
      if (need_insert && cohort< n_cohort) {
        
        d_old = d  # keep the current dose index at trigger time (needed for boundary ins_true_prob)
        
        ## BELOW MIN: insert half + treat, then choose again + treat (ADM style)
        if (insert_code == LOC_BELOW_MIN) {
          
          # need room for 2 extra cohorts (half-dose cohort + second inserted cohort)
          if (d == 1 && sum(n_treated_work) <= (n_patients - 2 * cohort_size)) {
            # record first insertion cohort (once, and only if insertion happens)
            if (insert_at_cohort[trial] == 0) {
              insert_at_cohort[trial] = sum(n_treated_work) / cohort_size
            }
            
            # insert half of current minimum
            half_dose = (dose_set_work[1] - dose_range[1]/(dose_range[2]-dose_range[1]))/2
            # Equivalent two-step transformation (std -> original -> /2 -> std), shown for clarity:
            # half_dose = (dose_set_work[1] * (dose_range[2] - dose_range[1]) + dose_range[1])/ 2 
            # half_dose = (half_dose -dose_range[1] )/ (dose_range[2] - dose_range[1])
            
            tox_half = ins_true_prob(
              insert_code = LOC_BELOW_MIN,
              d = d_old,
              p = tox_prob_work,
              target = target_prob,
              LOC_BELOW_MIN = LOC_BELOW_MIN,
              LOC_ABOVE_MAX = LOC_ABOVE_MAX
            )
            
            ins = insert_sorted(
              dose_set_work, tox_prob_work, n_dlt_work, n_treated_work, is_eliminated,
              newdose = half_dose, tox_new = tox_half
            )
            dose_set_work  = ins$dose_set_work
            tox_prob_work  = ins$tox_prob_work
            n_dlt_work     = ins$n_dlt_work
            n_treated_work = ins$n_treated_work
            is_eliminated  = ins$is_eliminated
            n_dose_work    = length(dose_set_work)
            
            next # go to treat one cohort at the inserted dose
          }
        } else {
          ## BETWEEN (i,i+1) or ABOVE MAX (metric one-shot)
          
          # Determine eligible interval endpoints (dl, dr)
          if (insert_code == LOC_ABOVE_MAX) {
            dl = tail(dose_set_work, 1)
            dr = (tail(dose_set, 1) * 1.5  - dose_range[1]) / (dose_range[2] - dose_range[1])
          } else {
            i = as.integer(insert_code)       # eligible interval is (i, i+1)
            dl = dose_set_work[i]
            dr = dose_set_work[i + 1]
          }
          
          if ((dr - dl) <= 0.002){
            # cat("trial = ", trial, ", dl = ", dl,  ", dr = ", dr, "\n")
            next
          } 
          
          if (insert_code == LOC_ABOVE_MAX) {
            # Above the current maximum, directly insert at 1.5x of the prespecified max dose.
            newdose = dr
          } else {
            ref_gap = min(diff(dose_set_work))
            sel = choose_newdose(
              dl = dl, dr = dr,
              dose_set_work = dose_set_work,
              n_dlt_work = n_dlt_work,
              n_treated_work = n_treated_work,
              pri_alpha = target_prob,
              pri_beta  = (1 - target_prob),
              key_L = key_L,
              key_U = key_U,
              symmetric = TRUE,
              k_left = 0.2,
              k_right = 0.2,
              ref_gap = ref_gap
            )
            newdose = sel$newdose
          }
          
          if (insert_at_cohort[trial] == 0L) {
            insert_at_cohort[trial] = sum(n_treated_work) / cohort_size
          }
          
          tox_new = ins_true_prob(
            insert_code = insert_code,
            d = d_old,
            p = tox_prob_work,
            target = target_prob,
            LOC_BELOW_MIN = LOC_BELOW_MIN,
            LOC_ABOVE_MAX = LOC_ABOVE_MAX
          )
          
          ins = insert_sorted(
            dose_set_work, tox_prob_work, n_dlt_work, n_treated_work, is_eliminated,
            newdose = newdose, tox_new = tox_new
          )
          dose_set_work  = ins$dose_set_work
          tox_prob_work  = ins$tox_prob_work
          n_dlt_work     = ins$n_dlt_work
          n_treated_work = ins$n_treated_work
          is_eliminated  = ins$is_eliminated
          n_dose_work    = length(dose_set_work)
          
          d = ins$d_idx
          next
        }
      }
      #-------------- end insertion block ----------------------------------------------------------------------
      
      
      #--------------------- SKBD elimination & move rule ----------------------------------------------
      # recompute BKP posterior
      ref_gap = min(diff(dose_set_work))
      post_work = post_par_all(
        n_dlt = n_dlt_work,
        n_treated = n_treated_work,
        dose_set = dose_set_work,
        pri_alpha = pri_alpha,
        pri_beta = pri_beta,
        symmetric = symmetric,
        k_left = k_left,
        k_right = k_right,
        ref_gap = ref_gap
      )
      post_alpha_work = post_work$post_alpha
      post_beta_work  = post_work$post_beta
      
      # whether the trial is overdose
      overdose_prob = 1 - pbeta(target_prob, post_alpha_work[d], post_beta_work[d]) # P(pi_d > target_prob | data)
      
      if(extra_safe){
        is_overdose = overdose_prob > cutoff_elimin - offset
      }else{
        is_overdose = overdose_prob > cutoff_elimin
      }
      
      # stopping rule
      if (n_treated_work[d] >= 3 && is_overdose) {
        is_eliminated[d:n_dose_work] = TRUE
        if (d == 1) {
          is_earlystop = TRUE
          break
        } else {
          ## enforce: eliminated dose cannot be assigned again
          d = d - 1
          next
        }
      }
      
      ## strongest key based on current posterior
      strong_key = get_strongKey(post_alpha_work[d], post_beta_work[d], keys, margin_left, margin_right)
      
      
      ## escalation/de-escalation rule
      if(strong_key < target_key && d < n_dose_work) {
        if (!is_eliminated[d + 1]) {
          d = d + 1
        }
      }else if(strong_key > target_key && d > 1) {
        d = d - 1
      }else{
        d = d
      }
    }
    #-------------------------------- end for --------------------------------------------
    
    #--------------------- final dose selection --------------------------------------------------------
    admissible_set = (n_treated_work > 0) & (!is_eliminated) # adimissble set
    adm_idx <- which(admissible_set)
    
    if(is_earlystop || length(adm_idx) == 0){
      sel_dose_idx[trial] = -1  # no dose should be selected as the MTD
      sel_dose[trial] = NA_real_
    }else{
      ref_gap = min(diff(dose_set_work))
      post_work = post_par_all(
        n_dlt = n_dlt_work, 
        n_treated = n_treated_work, 
        dose_set = dose_set_work, 
        pri_alpha = 0.01, 
        pri_beta = 0.01, 
        symmetric = TRUE, 
        k_left = 0.2,
        k_right = 0.2,
        ref_gap = ref_gap
      )
      post_alpha_work = post_work$post_alpha[adm_idx]
      post_beta_work  = post_work$post_beta[adm_idx]
      
      tox_prob_hat = post_alpha_work / (post_alpha_work + post_beta_work)
      tox_prob_hat = pava(tox_prob_hat)
      tox_prob_hat = tox_prob_hat + seq_along(tox_prob_hat) * 1e-10 # break ties by adding an increasingly small number
      
      sel_dose_idx[trial] = adm_idx[ which.min(abs(tox_prob_hat - target_prob))]
      sel_dose[trial] = dose_set_work[sel_dose_idx[trial]] * (dose_range[2] - dose_range[1]) + dose_range[1]
    }
    
    # save per-trial insertion count
    n_insertions[trial] = n_insertions_trial
    
    # save simdata rows
    temp = data.frame(
      Simulation = rep(trial, n_dose_work),
      Dose = dose_set_work * (dose_range[2] - dose_range[1]) + dose_range[1],
      N = n_treated_work,
      X = n_dlt_work,
      Selection = rep(sel_dose_idx[trial], n_dose_work)
    )
    simdata = rbind(simdata, temp)
    
    if (!light_return) {
      trial_detail[[trial]] = list(
        dose_set_work = dose_set_work,
        tox_prob_work = tox_prob_work,
        n_treated_work = n_treated_work,
        n_dlt_work = n_dlt_work,
        is_eliminated = is_eliminated,
        insert_at_cohort = insert_at_cohort[trial],
        n_insertions = n_insertions_trial,
        early_stop = is_earlystop
      )
    }
  } # end trial loop
  
  
  #--------------------------------- summary ----------------------------------------------------------
  # Percent selected / treated at each *prespecified* dose level
  sel_pct_prespec = numeric(n_dose)  # selection percentage (%)
  pts_pct_prespec = numeric(n_dose)  # patient allocation percentage (%)
  
  for (j in seq_len(n_dose)) {
    sel_pct_prespec[j] = mean(sel_dose == dose_set[j], na.rm = TRUE) * 100
    pts_pct_prespec[j] = sum(simdata$N[simdata$Dose == dose_set[j]]) / (n_trial * n_patients) * 100
  }
  
  # Percent selected / treated at *inserted* doses (complement to 100%)
  ins_sel_pct = 100 - sum(sel_pct_prespec)
  ins_pts_pct = 100 - sum(pts_pct_prespec)
  
  # Inserted dose values among selected doses (exclude NA and prespecified doses)
  ins_dose_vals = sel_dose[!is.na(sel_dose) & !(sel_dose %in% dose_set)]
  ins_dose_mean = mean(ins_dose_vals, na.rm = TRUE)
  ins_dose_sd   = sd(ins_dose_vals, na.rm = TRUE)
  
  # Insertion frequency and timing
  ins_trial_pct = mean(n_insertions > 0) * 100
  ins_cohort_vals = insert_at_cohort[insert_at_cohort > 0]
  ins_cohort_mean = if (length(ins_cohort_vals) > 0) mean(ins_cohort_vals) else NA_real_
  ins_n_median = median(n_insertions)
  
  out <- list(
    # prespecified-dose summaries
    sel_pct_prespec = sel_pct_prespec,
    pts_pct_prespec = pts_pct_prespec,
    
    # insertion summaries
    insertion = list(
      sel_pct     = ins_sel_pct,
      pts_pct     = ins_pts_pct,
      dose_mean   = ins_dose_mean,
      dose_sd     = ins_dose_sd,
      trial_pct   = ins_trial_pct,
      cohort_mean = ins_cohort_mean,
      n_median    = ins_n_median
    ),
    
    # others
    simdata = simdata,
    sel_dose_idx = sel_dose_idx,
    sel_dose = sel_dose,
    n_insertions = n_insertions,
    insert_at_cohort = insert_at_cohort
  )
  
  if (!light_return) out$trial_detail = trial_detail
  return(out)
}
