#' @title Operating Characteristics for the Inserted Shared Keyboard Design
#' 
#' @description
#' A short description...
#' 
#' @param target_prob description
#' 
#' @export

get_OC_Insert_SKBD <- function(
    target_prob, tox_prob, 
    n_cohort, cohort_size,
    dose_set, dose_range,
    C1 = 0.6, C2 = 0.6, 
    start_dose = 1,
    margin_left = 0.05, margin_right = 0.05, 
    n_earlystop = 1000, 
    cutoff_elimin = 0.95,
    extra_safe = FALSE, offset = 0.05, 
    symmetric = FALSE, theta = NULL,
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
  key_L <- target_prob - margin_left
  key_U <- target_prob + margin_right
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
  
  # dose_range interval
  if (!is.numeric(dose_range) || length(dose_range) != 2 || anyNA(dose_range)) {
    stop("`dose_range` must be a numeric vector of length 2: c(lower, upper).")
  }
  if (dose_range[1] >= dose_range[2]) stop("`dose_range[1]` must be < `dose_range[2]`.")
  if (min(dose_set) < dose_range[1] - 1e-12 || max(dose_set) > dose_range[2] + 1e-12) {
    stop("All `dose_set` values must lie inside `dose_range`.")
  }
  
  # start_dose index (in prespecified doses)
  if (!int1(start_dose) || start_dose < 1 || start_dose > length(tox_prob)) {
    stop("`start_dose` must be a single integer in 1:n_dose (indexing prespecified doses).")
  }
  
  # Standardize doses to [0,1]
  dose_set_std <- (dose_set - dose_range[1]) / (dose_range[2] - dose_range[1])
  
  # theta defaults
  if (is.null(theta)) {
    dose_diff_min = min(diff(dose_set_std))
    # defaults (safe-ish): asym borrows more from higher doses => theta1 > theta2
    theta <- if (symmetric) {
      -log(0.5) / dose_diff_min^2
    } else {
      c(-log(0.2), -log(0.8)) / dose_diff_min^2
    }
  }
  if (symmetric) {
    if (!is.numeric(theta) || length(theta) != 1 || is.na(theta) || theta <= 0) {
      stop("When `symmetric=TRUE`, `theta` must be a single positive number.")
    }
  } else {
    if (!is.numeric(theta) || length(theta) != 2 || anyNA(theta) || any(theta <= 0)) {
      stop("When `symmetric=FALSE`, `theta` must be a length-2 positive vector c(theta1, theta2).")
    }
    if (!(theta[1] > theta[2])) {
      stop("For asymmetric kernel, require theta[1] > theta[2] (borrow less from lower doses).")
    }
  }
  
  #------------------------------- Simulation containers ------------------------------------------
  n_dose     = length(tox_prob)       # number of dose levels
  n_patients = n_cohort * cohort_size # total number of patients
  
  sel_dose_idx     = integer(n_trial) # selected dose level index (j*)
  sel_dose         = numeric(n_trial) # selected dose value (d_j*)
  
  insert_at_cohort = integer(n_trial) # cohort index when first insertion occurs (0 if none)
  n_insertions     = integer(n_trial) # number of insertions in each trial
  
  # hyperparameters for beta prior
  pri_alpha = pri_beta = 0.5 
  
  # keep details optionally
  trial_detail = if (!light_return) vector("list", n_trial) else NULL
  
  ## get keys
  key_cutpoints  = get_Key(target_prob, margin_left, margin_right)
  target_key_idx = which.min(abs(key_cutpoints - (target_prob - margin_left)))
  
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
    
    ## treat sequentially by cohort until total N or early stop
    while(TRUE){
      
      # if enough patients have been treated at dose d, exit the loop
      if(n_treated_work[d] >= n_earlystop){
        break
      }
      
      # Treat one cohort at the current dose level d and generate DLT outcomes
      n_dlt_work[d] = n_dlt_work[d] + rbinom(1, cohort_size, tox_prob_work[d])
      n_treated_work[d] = n_treated_work[d] + cohort_size
      
      ## BKP posterior update for ALL working doses
      post_work = post_par_all(
        n_dlt = n_dlt_work, 
        n_treated = n_treated_work, 
        dose_set = dose_set_work, 
        pri_alpha = pri_alpha, 
        pri_beta = pri_beta, 
        symmetric = symmetric, 
        theta = theta)
      post_alpha_work = post_work$post_alpha
      post_beta_work  = post_work$post_beta

      #------------- determine if the dose insertion conditions are satisfied ------------------------
      res = insert_check(
        j = d,
        post_alpha = post_alpha_work,
        post_beta  = post_beta_work,
        key_L = key_L,
        key_U = key_U,
        C1 = C1,
        C2 = C2,
        n_treated = n_treated_work,
        use_pava = TRUE,
        return_probs = FALSE
      )
      
      need_insert = res$need_insert
      insert_code = res$insert_code
      LOC_BELOW_MIN = res$LOC_BELOW_MIN
      LOC_ABOVE_MAX = res$LOC_ABOVE_MAX
      
      #-------------- implement the dose insertion if conditions are satisfied ----------------------------- 
      if (need_insert && sum(n_treated_work) < n_patients) {
        
        tol = 1e-12
        inserted_any = FALSE
        d_old = d  # keep the current dose index at trigger time (needed for boundary ins_true_prob)
        
        # =========================
        # Case A: BELOW MIN  (ADM style: insert half, treat, then choose again)
        # =========================
        if (insert_code == LOC_BELOW_MIN) {
          # ADM-style feasibility: need room for 2 extra cohorts (half-dose cohort + second inserted cohort)
          if (abs(dose_set_work[d_old] - dose_set_work[1]) < tol &&
              sum(n_treated_work) <= (n_patients - 2 * cohort_size)) {
            
            if (insert_at_cohort[trial] == 0) {
              insert_at_cohort[trial] = sum(n_treated_work) / cohort_size
            }
            n_dose_work = n_dose_work + 1
            inserted_any = TRUE
            n_insertions_trial = n_insertions_trial + 1
            
            # insert half of current minimum
            halfdose = dose_set_work[1] / 2
            
            tox_half = ins_true_prob(
              insert_code = LOC_BELOW_MIN,
              d = d_old,
              p = tox_prob_work,
              target = target_prob,
              LOC_BELOW_MIN = LOC_BELOW_MIN,
              LOC_ABOVE_MAX = LOC_ABOVE_MAX
            )
            
            dose_set_work = c(halfdose, dose_set_work)
            tox_prob_work  = c(tox_half, tox_prob_work)
            
            # treat one cohort at the halfdose
            n_dlt_half = rbinom(1, cohort_size, tox_half)
            n_dlt_work     = c(n_dlt_half, n_dlt_work)
            n_treated_work = c(cohort_size, n_treated_work)
            is_eliminated  = c(FALSE, is_eliminated)

            # choose again within (halfdose, old d1) i.e., (dose_set_work[1], dose_set_work[2]) ----
            sel = choose_newdose(
              dl = dose_set_work[1],
              dr = dose_set_work[2],
              dose_set_work = dose_set_work,
              n_dlt_work = n_dlt_work,
              n_treated_work = n_treated_work,
              pri_alpha = pri_alpha,
              pri_beta  = pri_beta,
              key_L = key_L,
              key_U = key_U,
              symmetric = symmetric,
              theta = theta
            )
            newdose = sel$newdose
          } else {
            newdose =  NA
          }
        } else {
          # =========================
          # Case B: BETWEEN (i,i+1) or ABOVE MAX (metric one-shot)
          # =========================
          # Determine eligible interval endpoints (dl, dr)
          if (insert_code == LOC_ABOVE_MAX) {
            dl = dose_set_work[n_dose_work]
            dr = min(1, 1.5*dl) # for safety (ADM-style upper bound)
          } else {
            i = as.integer(insert_code)       # eligible interval is (i, i+1)
            dl = dose_set_work[i]
            dr = dose_set_work[i + 1]
          }

          sel = choose_newdose(
            dl = dl, dr = dr,
            dose_set_work = dose_set_work,
            n_dlt_work = n_dlt_work,
            n_treated_work = n_treated_work,
            pri_alpha = pri_alpha,
            pri_beta  = pri_beta,
            key_L = key_L,
            key_U = key_U,
            symmetric = symmetric,
            theta = theta
          )
          newdose = sel$newdose
        }
        
        # rearrange the dose and toxicity rate vectors to reflect the correct ordering
        if (!is.na(newdose)){
          
          if (insert_at_cohort[trial] == 0) {
            insert_at_cohort[trial] = sum(n_treated_work) / cohort_size
          }
          n_dose_work = n_dose_work + 1
          inserted_any = TRUE
          
          tox_new = ins_true_prob(
            insert_code = insert_code,
            d = d_old,
            p = tox_prob_work,
            target = target_prob,
            LOC_BELOW_MIN = LOC_BELOW_MIN,
            LOC_ABOVE_MAX = LOC_ABOVE_MAX
          )
          
          # treat one cohort at the newdose
          n_dlt_newdose = rbinom(1, cohort_size, tox_new)
          
          if (insert_code == LOC_BELOW_MIN){
            d = 2
            dose_set_work = c(dose_set_work[1], newdose, dose_set_work[-1])
            tox_prob_work = c(tox_prob_work[1], tox_new, tox_prob_work[-1])
            n_dlt_work     = c(n_dlt_work[1], n_dlt_newdose, n_dlt_work[-1])
            n_treated_work = c(n_treated_work[1], cohort_size, n_treated_work[-1])
            is_eliminated  = c(FALSE, is_eliminated)
            
          } else if (insert_code == LOC_ABOVE_MAX){
            d = n_dose_work
            dose_set_work = c(dose_set_work, newdose)
            tox_prob_work = c(tox_prob_work, tox_new)
            n_dlt_work     = c(n_dlt_work, n_dlt_newdose)
            n_treated_work = c(n_treated_work, cohort_size)
            is_eliminated  = c(is_eliminated, FALSE)
            
          } else {
            
            i = as.integer(insert_code)
            d = i + 1
            
            dose_set_work = c(dose_set_work[1:i], newdose, dose_set_work[-(1:i)])
            tox_prob_work = c(tox_prob_work[1:i], tox_new, tox_prob_work[-(1:i)])
            n_dlt_work     = c(n_dlt_work[1:i], n_dlt_newdose, n_dlt_work[-(1:i)])
            n_treated_work = c(n_treated_work[1:i], cohort_size, n_treated_work[-(1:i)])
            is_eliminated  = c(is_eliminated[1:i], FALSE, is_eliminated[1:i])
            
          }
          
          # Recompute BKP posterior for ALL doses after any insertion/treatment
          post_work = post_par_all(
            n_dlt = n_dlt_work, 
            n_treated = n_treated_work, 
            dose_set = dose_set_work, 
            pri_alpha = pri_alpha, 
            pri_beta = pri_beta, 
            symmetric = symmetric, 
            theta = theta)
          post_alpha_work = post_work$post_alpha
          post_beta_work  = post_work$post_beta
        }
      }
      #-------------- end insertion block ----------------------------------------------------------------------
      
      
      ## SKBD design: determine if the current dose should be eliminated
      # whether the trial is overdose
      overdose_prob = 1 - pbeta(target_prob, alpha_post_work[d], beta_post_work[d]) # P(pi_d > target_prob | data)
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
      strong_key = get_strongKey(alpha_post_work[d], beta_post_work[d], keys, margin_left, margin_right)
      
      
      ## escalation and de-escalation rule
      if(strong_key < target_key && d < n_dose_work) {
        if (!is_eliminated[d + 1]) {
          d = d + 1
        }
      }else if(strong_key > target_key && d > 1) {
        d = d - 1
      }else{
        d = d
      }

      # stop when total enrolled patients reaches n_patients
      total = sum(n_treated_work)
      if(total >= n_patients){
        break
      }  
    }
    #-------------------------------- end while --------------------------------------------
    
    # dose selection
    admissible_set = (n_treated_work > 0) & (!is_eliminated) # adimissble set
    adm_idx <- which(admissible_set)
    if(is_earlystop || length(adm_idx) == 0){
      sel_dose_idx[trial] = -1  # no dose should be selected as the MTD
    }else{
      tox_prob_hat = n_dlt_work[adm_idx] / n_treated_work[adm_idx]
      tox_prob_hat = pava(tox_prob_hat)
      
      # break ties by adding an increasingly small number
      tox_prob_hat = tox_prob_hat + seq_along(tox_prob_hat) * 1e-10
      sel_dose_idx[trial] = adm_idx[ which.min(abs(tox_prob_hat - target_prob)) ]
    }
    
    # save the results  
    n_insertions[trial] = n_insertions_trial
    sel_dose[trial] = ifelse(sel_dose_idx[trial] == 0, NA, dose_set_work[sel_dose_idx[trial]])
    temp_simdata = cbind(rep(trial, length(n_dlt_work)), 
                         dose_set_work, 
                         n_treated_work, 
                         n_dlt_work, 
                         rep(sel_dose_idx[trial], length(n_dose_work)))
    simdata = rbind(simdata, temp_simdata)
  }
  #--------------------------------- end main simulation loop ----------------------------------------
  
  
  ## output results
  simdata = data.frame(simdata)
  names(simdata) = c("Simulation", "Dose", "N", "X", "Selection")
  selpercent = rep(0, n_dose) # selection percentage at prespecified doses
  ptspercent = rep(0, n_dose) # percentage of patients at prespecified doses
  for (i in 1:n_dose) {
    selpercent[i] = sum(sel_dose == dose_set_std[i], na.rm=TRUE) / n_trial * 100
    ptspercent[i] = sum(simdata[which(simdata$Dose == dose_set_std[i]), 3]) / (n_trial * n_patients) *100
  }
  ins.select = 100 - sum(selpercent) # selection percentage at inserted dose
  ins.pts = 100 - sum(ptspercent)    # percentage of patients at inserted dose
  
  ins.dose = NULL
  for (i in 1:n_trial){
    if (sel_dose[i] %in% dose_set_std) {
      ins.dose = ins.dose
    } else {
      ins.dose = c(ins.dose, sel_dose[i])
    }
  }
  ins.mean = mean(ins.dose, na.rm=TRUE)
  ins.sd = sd(ins.dose, na.rm=TRUE)
  ins.percent = (1000 - sum(n_insertions == 0)) / n_trial * 100 
  cohort.mean = mean(insert_at_cohort[which(insert_at_cohort!=0)])
  insTimes.median <- median(n_insertions)
  cat ("selection percentage at each prespecified dose level (%):\n")
  cat (formatC(selpercent,digits=2, format="f"),sep="  ","\n")
  cat ("percentage of patients treated at prespecified dose level (%):\n")
  cat (formatC(ptspercent,digits=2, format="f"),sep="  ","\n")
  cat ("mean of the inserted dose (SD):",formatC(ins.mean,digits=2,format="f"),"(",formatC(ins.sd, digits=2, format="f"),") \n")
  cat ("selection percentage of the inserted dose (%):",formatC(ins.select,digits=2, format="f"),"\n")
  cat ("percentage of patients treated at the inserted doses (%):",formatC(ins.pts, digits=2, format="f"),"\n")
  cat ("percentage of trials with dose insertion (%):",formatC(ins.percent,digits=1, format="f"),"\n")
  cat ("The average of cohorts after which trial insertion is likely to take place:",formatC(cohort.mean,digits=1,format="f"),"\n")
  cat ("Median number of insertion times:", formatC(insTimes.median,digits=1,format="f"),"\n")
}