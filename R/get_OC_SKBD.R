#' Operating Characteristics for the Shared Keyboard Design
#'
#' Simulate Phase I dose-finding trials under the shared keyboard design (SKBD)
#' and summarize operating characteristics (e.g., PCS, PCA, overdosing risk).
#' SKBD follows the keyboard decision rules (strongest key vs. target key) but
#' constructs posterior toxicity at the current dose using a kernel-weighted
#' Beta update (Beta Kernel Process-style borrowing) when \code{shared=TRUE}.
#'
#' @param target_prob A scalar in \eqn{(0,1)} giving the target DLT probability \eqn{\phi}.
#' @param tox_prob A numeric vector of length \code{n_dose} giving the true DLT
#'   probabilities at each prespecified dose level (used for simulation).
#' @param n_cohort Integer; number of cohorts per simulated trial.
#' @param cohort_size Integer; number of patients per cohort.
#' @param dose_set Numeric vector of length \code{length(tox_prob)} giving dose
#'   locations (must be ordered increasingly). Used for dose standardization and
#'   kernel construction when \code{shared=TRUE}.
#' @param start_dose Integer in \code{1:n_dose}; starting dose index for each trial.
#' @param margin_left Nonnegative scalar; left tolerance \eqn{\varepsilon_1} for the
#'   target key \eqn{(\phi-\varepsilon_1,\phi+\varepsilon_2)}.
#' @param margin_right Nonnegative scalar; right tolerance \eqn{\varepsilon_2} for the
#'   target key \eqn{(\phi-\varepsilon_1,\phi+\varepsilon_2)}.
#' @param n_earlystop Integer; if the number treated at the current dose reaches
#'   \code{n_earlystop}, the trial loop stops (useful for preventing excessive
#'   sampling at one dose in simulations).
#' @param cutoff_elimin Scalar in \eqn{(0,1)}; overdose-control cutoff. A dose is
#'   eliminated if \eqn{Pr(\pi_j>\phi\mid\mathcal{D})} exceeds this cutoff and
#'   \code{n_j >= 3}.
#' @param extra_safe Logical; if \code{TRUE}, uses a more conservative elimination
#'   threshold \code{cutoff_elimin - offset}.
#' @param offset Nonnegative scalar; safety offset used only when \code{extra_safe=TRUE}.
#' @param shared Logical; if \code{TRUE}, use kernel-weighted borrowing across doses
#'   to construct the posterior at the current dose. If \code{FALSE}, reduces to
#'   local updating (keyboard-style) via an identity kernel.
#' @param symmetric Logical; kernel type when \code{shared=TRUE}. If \code{TRUE},
#'   uses a symmetric Gaussian kernel; otherwise uses an asymmetric Gaussian kernel.
#' @param theta Kernel parameter(s). If \code{NULL}, calibrated from the minimum
#'   standardized dose spacing. For symmetric kernels, \code{theta} is a scalar; for
#'   asymmetric kernels, \code{theta} is a length-2 vector \code{c(theta_left, theta_right)}.
#' @param light_return Logical; if \code{TRUE}, do not store/return individual-level
#'   dose/DLT paths to reduce memory usage. If \code{FALSE}, returns
#'   \code{dose_Paths} and \code{DLT_Paths}.
#' @param n_trial Integer; number of simulated trials.
#' @param seed Integer; random seed for reproducibility.
#'
#' @details
#' In each simulated trial, cohorts are treated sequentially. After each cohort,
#' the posterior toxicity probability at the current dose is updated:
#' \itemize{
#' \item If \code{shared=FALSE}, the update uses only the data at the current dose.
#' \item If \code{shared=TRUE}, the update uses normalized kernel weights over the
#' set of doses with observed data, borrowing more from nearby doses (and optionally
#' asymmetrically from higher vs. lower doses).
#' }
#' The design identifies the strongest key (toxicity interval with largest posterior
#' probability) and applies the standard keyboard escalation/de-escalation rule relative
#' to the target key. An overdose-control rule eliminates overly toxic doses and all
#' higher doses. At the end of the trial, the MTD is selected by applying isotonic
#' regression (via \code{pava()}, assumed available in the package) to observed rates
#' among admissible doses and choosing the dose closest to \code{target_prob}.
#'
#' @return A list with components:
#' \itemize{
#' \item \code{PCS}: percent correct selection (in \%).
#' \item \code{PCA}: percent correct allocation (in \%).
#' \item \code{above_MTD}: percent treated above the MTD (in \%).
#' \item \code{ROD60}, \code{ROD80}: risk of overdosing (in \%), defined as the percent of trials
#'   with >60\% or >80\% patients treated above the MTD.
#' \item \code{dose_select}: length-\code{n_trial} vector of selected MTD indices; \code{-1} indicates no MTD.
#' \item \code{select_percent}: selection percentages for each dose and \code{-1}.
#' \item \code{n_patient_mean}: mean number of patients per trial.
#' \item \code{n_DLT}: mean number of DLTs per trial.
#' \item \code{Y}, \code{N}: \code{n_trial x n_dose} matrices of DLT counts and treated counts.
#' \item \code{dose_Paths}, \code{DLT_Paths}: only returned when \code{light_return=FALSE};
#'   \code{n_trial x (n_cohort*cohort_size)} matrices recording individual-level assignments and outcomes.
#' }
#'
#' @references
#' Yan, F., Mandrekar, S. J., and Yuan, Y. (2017).
#' Keyboard: A Novel Bayesian Toxicity Probability Interval Design for Phase I Clinical Trials.
#' \emph{Clinical Cancer Research}, \bold{23}(15), 3994--4003. doi:10.1158/1078-0432.CCR-17-0220.
#'
#' @examples
#' \dontrun{
#' res <- get_OC_SKBD(
#'   target_prob = 0.30,
#'   tox_prob = c(0.05, 0.12, 0.30, 0.45, 0.60),
#'   n_cohort = 10, cohort_size = 3
#' )
#' str(res)
#' }
#'
#' @export

get_OC_SKBD <- function(
    target_prob, tox_prob, 
    n_cohort, cohort_size,
    dose_set = 1:length(tox_prob), 
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
  
  ## ---- argument checks ----
  if (!is.numeric(target_prob) || length(target_prob) != 1 || target_prob <= 0 || target_prob >= 1) {
    stop("`target_prob` must be a single number in (0, 1).")
  }
  
  if (target_prob < 0.05) {
    warning("`target_prob` is very low (< 0.05); ",
            "operating characteristics may be unstable.")
  }
  
  if (target_prob > 0.5) {
    warning("`target_prob` is very high (> 0.5); ",
            "this is uncommon in dose-finding designs.")
  }
  
  if (!is.numeric(tox_prob) || any(tox_prob < 0 | tox_prob > 1)) {
    stop("`tox_prob` must be numeric in [0, 1].")
  }
  
  if (offset >= 0.5) {
    stop("`offset` must be smaller than 0.5.")
  }
  
  if (extra_safe && (cutoff_elimin - offset) <= 0) {
    stop("When `extra_safe=TRUE`, require `cutoff_elimin - offset > 0`.")
  }
  
  if (n_earlystop <= 6) {
    warning(
      "`n_earlystop` is too small to ensure good operating characteristics. ",
      "Recommended range is 9-18."
    )
  }
  
  ## check: dose_set length matches tox_prob
  if (length(dose_set) != length(tox_prob)) {
    stop("`dose_set` must have the same length as `tox_prob` (i.e., n_dose).")
  }
  
  ## check: start_dose is a valid dose index
  if (!is.numeric(start_dose) || length(start_dose) != 1 ||
      is.na(start_dose) || start_dose %% 1 != 0 ||
      start_dose < 1 || start_dose > length(tox_prob)) {
    stop("`start_dose` must be a single integer in 1:n_dose.")
  }
  
  n_dose = length(tox_prob)                                   # number of dose
  n_patient = n_cohort * cohort_size                          # number of patients
  Y = N = matrix(NA, nrow = n_trial, ncol = n_dose)           # Y: number of DLT response at each dose; N: number of patients at each dose
  if (!light_return) {
    dose_Paths = matrix(NA, nrow = n_trial, ncol = n_patient) # record the dose path of the dose assignment for all patients
    DLT_Paths = matrix(NA, nrow = n_trial, ncol = n_patient)  # record the DLT path of the dose assignment for all patients
  } else {
    dose_Paths <- NULL
    DLT_Paths  <- NULL
  }
  dose_select = rep(NA, n_trial)                              # MTD selection                           
  
  
  ## proir setting: non-informative
  # r0 = 3
  # alpha_pri = rep((r0-2)*target_prob + 1, n_dose)
  # beta_pri = r0 - alpha_pri
  alpha_pri = beta_pri = rep(1, n_dose)

  ## get keys
  keys = get_Key(target_prob, margin_left, margin_right)
  target_key = which.min(abs(keys - (target_prob - margin_left)))
  
  ## dose standardization
  dose_set_std = (dose_set - min(dose_set)) / (max(dose_set) - min(dose_set))
  
  ## get kernel
  if(shared){
    if(is.null(theta)){
      dose_diff_min = min(diff(dose_set_std))
      if(symmetric){
        theta = -log(0.8) / dose_diff_min^2
      }else{
        # theta = c(-log(0.3 - target_prob / 2), -log(0.7 + target_prob / 2)) / dose_diff_min^2
        # theta = c(-log(target_prob), -log(1 - target_prob)) / dose_diff_min^2
        theta = c(-log(0.2), -log(0.8)) / dose_diff_min^2
      }
    }
    ker_vals = matrix(0, nrow = n_dose, ncol = n_dose)
    for (i in 1:n_dose) {
      ker_vals[i, ] = kernel(dose_set_std[i], dose_set_std, symmetric, theta)
    }
  } else {
    ker_vals = diag(n_dose)
  }

  #------------------------ begin trials -------------------------#
  for (trial in 1:n_trial) {
    ## initial for each trial
    d = start_dose
    y = n = rep(0, n_dose)
    is_earlystop = FALSE               # Whether to stop the design early
    is_eliminated = rep(FALSE, n_dose) 
    alpha_post = alpha_pri
    beta_post  = beta_pri
    
    #-------------------------------- begin one trial --------------------------------------#
    for (cohort in 1:n_cohort) {
      if (n[d] >= n_earlystop) {
        # is_earlystop = TRUE
        break
      }
      
      # DLT = rbinom(cohort_size, 1, tox_prob[d]) 
      DLT = runif(cohort_size) < tox_prob[d]   # DLT in one cohort
      y[d] = y[d] + sum(DLT)                   # the number of DLT responses
      n[d] = n[d] + cohort_size                # the number of patients
      
      ## record the dose and DLT path
      if (!light_return) {
        path_index = ((cohort-1)*cohort_size+1) : (cohort*cohort_size)
        DLT_Paths[trial, path_index] = DLT
        dose_Paths[trial, path_index] = d
      }
      
      # weight = ker_vals[d, ] / sum(ker_vals[d, ])
      # alpha_post[d] = alpha_pri[d] + sum(weight * y)
      # beta_post[d] = beta_pri[d] + sum(weight * (n - y))
      weight = ker_vals[d, n>0] / sum(ker_vals[d, n>0])
      alpha_post[d] = alpha_pri[d] + sum(weight * y[n>0])
      beta_post[d] = beta_pri[d] + sum(weight * (n[n>0] - y[n>0]))
      
      # curve(dbeta(x, shape1 = alpha_post[d], shape2 = beta_post[d]), 0.001, 0.999)
      # abline(v = c(target_prob - margin_left, target_prob, target_prob + margin_right), lty = 2)
      
      # whether the trial is overdose
      overdose_prob = 1 - pbeta(target_prob, alpha_post[d], beta_post[d]) # P(pi_d > target_prob | data)
      if(extra_safe){
        is_overdose = overdose_prob > cutoff_elimin - offset
      }else{
        is_overdose = overdose_prob > cutoff_elimin
      }
      
      # stopping rule
      if (n[d] >= 3 && is_overdose) {
        is_eliminated[d:n_dose] = TRUE
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
      strong_key = get_strongKey(alpha_post[d], beta_post[d], keys, margin_left, margin_right)
      
      
      ## escalation and de-escalation rule
      if(strong_key < target_key && d < n_dose) {
        if (!is_eliminated[d + 1]) {
          d = d + 1
        }
      }else if(strong_key > target_key && d > 1) {
        d = d - 1
      }else{
        d = d
      }
    } 
    #---------------------------- end one trial -------------------------------------------#
    
    ## record data
    Y[trial, ] = y
    N[trial, ] = n
    
    ## Maximum Tolerated Dose (MTD) Selection
    admissible_set = (n > 0) & (!is_eliminated) # adimissble set
    adm_idx <- which(admissible_set)
    if(is_earlystop || length(adm_idx) == 0){
      dose_select[trial] = -1  # no dose should be selected as the MTD
    }else{
      # ## poster mean and variance of toxicity probabilities using beta(0.05, 0.05) as the prior
      # tox_prob_hat = (y[admissible_set] + 0.05) / (n[admissible_set] + 0.1)
      # tox_prob_hat_var = (y[admissible_set] + 0.05) * (n[admissible_set] - y[admissible_set] + 0.05)/(
      #   (n[admissible_set] + 0.1)^2 * (n[admissible_set] + 0.1 + 1))
      # ## perform the isotonic transformation using PAVA
      # tox_prob_hat = pava(tox_prob_hat, wt = tox_prob_hat_var)

      # if(shared){
      #   tox_prob_hat = rep(0, length(admissible_set))
      #   # get kernel
      #   theta_PCS = -log(0.2) / dose_diff_min^2
      #   for (i in 1:length(admissible_set)) {
      #     ker_vals_PCS = kernel(dose_set_std[i], dose_set_std, 
      #                          symmetric=TRUE, theta_PCS)
      #     weight = ker_vals_PCS[n>0] / sum(ker_vals_PCS[n>0])
      #     # alpha_pri = beta_pri = 1
      #     tox_prob_hat[i] = sum(weight * y[n>0]) / sum(weight * n[n>0])
      #   }
      # } else {
      #   # same with the original keyboard design
      #   tox_prob_hat = y[admissible_set] / n[admissible_set]
      # }
      
      # same with the original keyboard design
      tox_prob_hat = y[adm_idx] / n[adm_idx]
      
      # # whether or not monotonic of the estimated toxicity probability
      # is_monotonic = (all(diff(tox_prob_hat) > 0))
      
      tox_prob_hat = pava(tox_prob_hat)
      # break ties by adding an increasingly small number
      tox_prob_hat = tox_prob_hat + seq_along(tox_prob_hat) * 1e-10
      
      # select dose closest to the target_prob as the MTD
      dose_select[trial] = adm_idx[ which.min(abs(tox_prob_hat - target_prob)) ]
    }# end MTD selection 
  }
  #------------------------ end trials -------------------------#
  
  #------------------------ begin performance metrics -----------------------#
  ## 1. accuracy
  # 1.1 percentage of correct selection (PCS)
  select_percent = as.vector(table(factor(dose_select, levels = c(1:n_dose, -1)))) / n_trial * 100
  
  is_all_ovedose = tox_prob[1] > target_prob + margin_right # ref. 18_Zhou_SiM & 18_Zhou_CCR
  if(is_all_ovedose){
    target_dose = n_dose + 1  # no dose should be selected as the MTD
  }else{
    target_dose = which.min(abs(tox_prob - target_prob))
  }
  PCS = select_percent[target_dose]
  
  # 1.2 percentage of correct allocation (PCA)
  n_patient_each_dose = colMeans(N)
  n_patient_mean = sum(N) / n_trial
  if(is_all_ovedose){
    PCA = 0
  }else{
    PCA = n_patient_each_dose[target_dose] / n_patient_mean * 100
  }
  
  ## 2. safty
  # 2.1 percentage of patients treated above the MTD (above_MTD)
  above_MTD = sum(n_patient_each_dose[tox_prob > (target_prob + margin_right)] / n_patient_mean * 100)
  
  # 2.2 risk of overdosing (ROD)
  ROD60 = mean(rowSums(N[, tox_prob > (target_prob + margin_right), drop = FALSE]) > 0.6 * rowSums(N)) * 100
  ROD80 = mean(rowSums(N[, tox_prob > (target_prob + margin_right), drop = FALSE]) > 0.8 * rowSums(N)) * 100
  
  # 2.3 the number of DLT
  n_DLT = mean(rowSums(Y))
  
  #------------------------ end performance metrics -----------------------#
  
  
  out = list(PCS = PCS, PCA = PCA, 
             above_MTD = above_MTD, ROD60 = ROD60, ROD80 = ROD80,
             dose_select = dose_select, select_percent = select_percent, 
             n_patient_mean = n_patient_mean, n_DLT = n_DLT,
             Y = Y, N = N)
  if (!light_return) {
    out$dose_Paths <- dose_Paths
    out$DLT_Paths  <- DLT_Paths
  }
  
  return(out)
}

