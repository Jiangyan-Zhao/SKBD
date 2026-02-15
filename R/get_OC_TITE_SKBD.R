#' @title Operating Characteristics for the Time-to-Event Shared Keyboard Design
#' 
#' @description
#' A short description...
#' 
#' @param target_prob The target dose-limiting toxicity (DLT) rate
#' @param k_left Numeric scalar in `(0,1)`. Left-side neighbor borrowing strength passed to `kernel_fun()`.
#' @param k_right Numeric scalar in `(0,1)`. Right-side neighbor borrowing strength passed to `kernel_fun()`
#'        (and used for symmetric borrowing when `symmetric=TRUE`).
#' @param ref_gap Optional positive scalar. Reference spacing passed to `kernel_fun()`. If `NULL`,
#'        kernel defaults to the minimum adjacent spacing in `dose_set`.
#' @param prior_p a vector of length 3, which specifies the prior probability 
#'                that the time to toxicity lies inside the time interval 
#'                (0,\code{tau}/3), (\code{tau}/3,2*\code{tau}/3), (2*\code{tau}/3,1). 
#'                The default value is \code{prior_p=c(1/3,1/3,1/3)}. 
#' @param accrual the accrual rate, i.e., the number of patients accrued in 1 unit of time
#' @export

get_OC_TITE_SKBD <- function(target_prob, tox_prob, 
                        n_cohort, cohort_size,
                        symmetric = FALSE,
                        k_left = 0.2, k_right = 0.8, ref_gap = NULL,
                        shared = TRUE, # incorporate the keybooard design
                        dose_set = 1:length(tox_prob), 
                        n_earlystop = 1000, start_dose = 1,
                        margin_left = 0.05, margin_right = 0.05, 
                        cutoff_elimin = 0.95,
                        extraSafe = FALSE, offset = 0.05, 
                        tau = 3, accrual = 2, alpha = 0.5,
                        piecewise = FALSE, prior_p = rep(1/3,3), 
                        dist_DLT = c("weibull", "loglogistic", "uniform"),
                        dist_enter = c("exp", "uniform"),
                        n_trial = 1000, seed = 6) {
  set.seed(seed)
  
  dist_DLT = match.arg(dist_DLT)
  dist_enter = match.arg(dist_enter)
  
  ## ---- argument checks ----
  if (target_prob <= 0 || target_prob >= 1) {
    stop("`target_prob` must be in (0, 1).")
  }
  
  if (target_prob < 0.05) {
    warning("`target_prob` is very low (< 0.05); ",
            "operating characteristics may be unstable.")
  }
  
  if (target_prob > 0.5) {
    warning("`target_prob` is very high (> 0.5); ",
            "this is uncommon in dose-finding designs.")
  }
  
  if (offset >= 0.5) {
    stop("`offset` must be smaller than 0.5.")
  }
  
  if (n_earlystop <= 6) {
    warning(
      "`n_earlystop` is too small to ensure good operating characteristics. ",
      "Recommended range is 9–18."
    )
  }
  
  if (tau <= 0) {
    stop("`tau` (maximum follow-up time) must be positive.")
  }
  
  if (accrual <= 0) {
    stop("`accrual` (patient accrual rate) must be positive.")
  }
  
  if (length(prior_p) != 3 || any(prior_p < 0)) {
    stop("`prior_p` must be a non-negative vector of length 3.")
  }
  
  if (abs(sum(prior_p) - 1) > 1e-8) {
    warning("`prior_p` does not sum to 1; it will be normalized internally.")
    prior_p = prior_p / sum(prior_p)
  }
  
  if (!all(dist_DLT %in% c("uniform", "weibull", "loglogistic"))) {
    stop("`dist_DLT` must be one of 'uniform', 'weibull', or 'loglogistic'.")
  }
  
  if (!all(dist_enter %in% c("uniform", "exp"))) {
    stop("`dist_enter` must be either 'uniform' or 'exp'.")
  }
  
  margin_length = margin_left + margin_right      # the length of target key
  
  #-------------------------- begin setting -------------------------------#
  n_dose = length(tox_prob)                                   # number of dose
  n_patient = n_cohort * cohort_size                          # number of patients
  
  Y = N = matrix(NA, nrow = n_trial, ncol = n_dose)           # Y: number of DLT response at each dose; N: number of patients at each dose
  dose_Paths = matrix(NA, nrow = n_trial, ncol = n_patient)   # record the dose path of the dose assignment for all patients
  DLT_Paths = matrix(NA, nrow = n_trial, ncol = n_patient)    # record the DLT path of the dose assignment for all patients
  dose_select = rep(NA, n_trial)                              # MTD selection                           
  
  duration = rep(0, n_trial)
  
  ## proir setting: non-informative
  # r0 = 2
  # pri_alpha = rep((r0-2)*target_prob + 1, n_dose)
  # pri_beta = r0 - pri_alpha
  pri_alpha = pri_beta = rep(1, n_dose)
  post_alpha = post_beta = rep(1, n_dose)
  
  ## dose standardization
  dose_set_std = (dose_set - min(dose_set)) / (max(dose_set) - min(dose_set))
  
  ## get kernel
  if(shared){
    ker_vals = matrix(0, nrow = n_dose, ncol = n_dose)
    for (i in 1:n_dose) {
      ker_vals[i, ] = kernel_fun(
        dose = dose_set_std[i],
        dose_set = dose_set_std,
        symmetric = symmetric,
        k_left = k_left,
        k_right = k_right,
        ref_gap = ref_gap
      )
    }
  } else {
    ker_vals = diag(n_dose)
  }
  #-------------------------- end setting -------------------------------#
  
  
  #------------------------ begin trials -------------------------#
  for (trial in 1:n_trial) {
    ## initial for each trial
    d = start_dose
    is_earlystop = FALSE               # Whether to stop the design early
    eliminate = rep(FALSE, n_dose) 
    
    ## recode the times
    time_DLT = rep(NA, n_patient)        # t_i: the acutal time to DLT
    time_enroll = rep(NA, n_patient)     # the enrollment time
    time_decision = 0                    # the decision making time

    ## recode the effective data
    y_tilde = rep(0, n_dose)        # \tlide{y}_j: # patients experienced DLT by the decision marking time
    m = rep(0, n_dose)              # m_j: # patients completed the DLT assessment without experiencing a DLT
    m_tilde = rep(0, n_dose)        # \tlide{m}_j: the effective number of non-DLTs
    n_tilde = rep(0, n_dose)        # \tlide{n}_j = \tlide{y}_j + \tlide{m}_j: ESS (effective sample size)
    n_ascertaining = rep(0, n_dose) # r_j: # ascertaining patients
    n_pending = rep(0, n_dose)      # c_j: # pending patients
    y = rep(0, n_dose)              # # DLT responses
    n = rep(0, n_dose)              # # patients
    
    #-------------------------------- begin one trial --------------------------------------#
    for(cohort in 1:n_cohort) {
      current_cohort_index = ((cohort-1)*cohort_size+1) : (cohort*cohort_size)
      
      # generate data for the new patient
      obs_cohort = gen_tite(n = cohort_size, pi = tox_prob[d], 
                            dist = dist_DLT, alpha = alpha, tau = tau)
      y[d] = y[d] + sum(obs_cohort$DLT)
      n[d] = n[d] + cohort_size
      
      ## record the dose and DLT path
      DLT_Paths[trial, current_cohort_index] = obs_cohort$DLT
      dose_Paths[trial, current_cohort_index] = rep(d, cohort_size)
      
      # check whether the trial should be early terminated
      if (n[d] >= n_earlystop) {
        # is_earlystop = TRUE
        break
      }
      
      # record the times
      if(dist_enter == "exp"){
        time_enroll[current_cohort_index] = time_decision +
          cumsum(c(0, rexp(n = cohort_size - 1, rate = accrual)))
      } else if(dist_enter == "uniform"){
        time_enroll[current_cohort_index] = time_decision + 
          cumsum(c(0, runif(n = cohort_size - 1, min = 0, max = 2 / accrual)))
      }
      time_decision = max(time_enroll, na.rm = TRUE) # No decision required, just record
      if(cohort == n_cohort){
        duration[trial] = time_decision + tau
      }
      
      time_DLT[current_cohort_index] = time_enroll[current_cohort_index] + obs_cohort$time_DLT
      
      # suspension rule: Forpatient safety, we require that dose escalation is not allowed until 
      # at least two patients have completed the DLT assessment at the current dose level. 
      suspension = TRUE
      while (suspension) {
        suspension = FALSE
        if(dist_enter == "exp"){
          time_decision = time_decision + rexp(n = 1, rate = accrual)
        } else if(dist_enter == "uniform"){
          time_decision = time_decision + runif(n = 1, min = 0, max = 2 / accrual)
        } 
        
        # u_i: the observed follow-up time
        time_follow_up = time_decision - time_enroll 
        
        # Determine whether the assessment has been completed.
        delta = (time_DLT <=  time_decision) | (time_follow_up > tau)
        
        # update the effective dataset
        n_dose_exist = max(dose_Paths[trial, ], na.rm = TRUE)
        for (d_j in 1:n_dose_exist) {
          d_j_idx = which(dose_Paths[trial, ] == d_j)
          
          DLT_j = DLT_Paths[trial, d_j_idx]
          delta_j = delta[d_j_idx]
          
          y_tilde[d_j] = sum(delta_j * DLT_j)
          m[d_j] = sum(delta_j * (1 - DLT_j))
          n_ascertaining[d_j] = sum(delta_j)
          n_pending[d_j] = n[d_j] - n_ascertaining[d_j]

          omega = rep(1, n[d_j])
          ## uniform distribution
          time_follow_up_j = time_follow_up[d_j_idx]         # u_i
          omega[!delta_j] = time_follow_up_j[!delta_j] / tau # omega_i = u_i / tau
          if(piecewise){
            ## piecewise uniform distribution
            omega[!delta_j] = (3 * prior_p[1] * omega[!delta_j]) * (
              time_follow_up_j[!delta_j] <= tau / 3
              ) + (
                prior_p[1] - prior_p[2] + 3 * prior_p[2] * omega[!delta_j]
                ) * (
                  tau / 3 < time_follow_up_j[!delta_j] & time_follow_up_j[!delta_j] <= 2 / 3 * tau
                  ) + (
                    prior_p[1] + prior_p[2] - 2 * prior_p[3] + 3 * prior_p[3] * omega[!delta_j]
                    ) * (
                      2 / 3 * tau < time_follow_up_j[!delta_j] & time_follow_up_j[!delta_j] <= tau
                      )
          }
          
          m_tilde[d_j] = m[d_j] + sum((1 - delta_j) * omega)
          n_tilde[d_j] = y_tilde[d_j] + m_tilde[d_j]
        }
        
        # update the posterior parameter
        for (d_j in 1:n_dose_exist) {
          weight = ker_vals[d_j, ] / sum(ker_vals[d_j, n > 0])
          # weight = ker_vals[d_j, ] / sum(ker_vals[d_j, ])
          
          post_alpha[d_j] = pri_alpha[d_j] + sum(weight * y_tilde)
          
          # post_beta[d_j] = pri_beta[d_j] + sum(weight * m_tilde)
          
          m_1 = m; m_1[d_j] = m_tilde[d_j] # Only share ascertained data
          post_beta[d_j] = pri_beta[d_j] + sum(weight * m_1)
          
        }
        
        # determine which dose level should be eliminated
        for(d_j in 1:n_dose_exist){
          overdose_prob = 1 - pbeta(target_prob, post_alpha[d_j], post_beta[d_j])
          if(extraSafe){
            is_overdose = overdose_prob > cutoff_elimin - offset
          } else {
            is_overdose = overdose_prob > cutoff_elimin
          }
          if (n_ascertaining[d_j] >= 3 & is_overdose){ # stopping rule
            eliminate[d_j:n_dose] = TRUE
            if (d_j == 1) {
              is_earlystop = TRUE
            }
            break # Break out of the for loop
          }
        }
        
        # check whether the current dose level should be eliminated
        if(eliminate[d]) {
          d = which(eliminate == TRUE)[1] - 1
          if(d == 0){
            is_earlystop = TRUE
          }
          break # Break out of the while loop
        }
        
        #check whether the current dose is toxic based on observed data
        target_key_prob = pbeta(target_prob + margin_right, post_alpha[d], post_beta[d]) - 
          pbeta(target_prob - margin_left, post_alpha[d], post_beta[d])
        
        left_key_prob = pbeta(target_prob - margin_left, post_alpha[d], post_beta[d]) - 
          pbeta(target_prob - margin_left - margin_length, post_alpha[d], post_beta[d])
        
        right_key_prob = pbeta(target_prob + margin_right + margin_length, post_alpha[d], post_beta[d]) - 
          pbeta(target_prob + margin_right, post_alpha[d], post_beta[d])
        
        # make dose assignment decisions
        is_escalation = (left_key_prob > target_key_prob) && (left_key_prob > right_key_prob) && (d < n_dose)
        is_de_escalation = (right_key_prob > target_key_prob) && (right_key_prob > left_key_prob) && (d > 1)
        if(is_escalation){
          if(n_ascertaining[d] < 2){ # suspension rule
            suspension = TRUE
          }else{
            if (!eliminate[d + 1]) {
              d = d + 1
            }
          }
        }else if (is_de_escalation){
          d = d - 1
        }
      } # end while
      
      if(is_earlystop){
        break
      }
    }
    #---------------------------- end one trial -------------------------------------------#

    ## record data
    Y[trial, ] = y
    N[trial, ] = n
    
    ## Maximum Tolerated Dose (MTD) Selection
    admissble_set = (n > 0) & (!eliminate) # adimissble set
    if(is_earlystop){
      dose_select[trial] = -1  # no dose should be selected as the MTD
    }else{
      # same with the original keyboard design
      tox_prob_hat = y[admissble_set] / n[admissble_set]
      tox_prob_hat = pava(tox_prob_hat)
      # break ties by adding an increasingly small number
      tox_prob_hat = tox_prob_hat + (1:length(tox_prob_hat)) * 1e-10 
      
      # select dose closest to the target_prob as the MTD
      dose_select[trial] = which.min(abs(tox_prob_hat - target_prob))
    }# end MTD selection 
  }
  #------------------------ end trials -------------------------#
  
  #------------------------ begin performance metrics -----------------------#
  ## 1. accuracy
  # 1.1 percentage of correct selection (PCS)
  select_percent = as.vector(table(factor(dose_select, levels = c(1:n_dose, -1)))) / n_trial * 100
  is_all_ovedose = tox_prob[1] > target_prob + 0.1 # ref. 18_Zhou_SiM & 18_Zhou_CCR
  if(is_all_ovedose){
    target_dose = n_dose + 1  # no dose should be selected as the MTD
  }else{
    target_dose = which.min(abs(tox_prob - target_prob))
  }
  PCS = select_percent[target_dose]
  
  # 1.2 percentage of patients treated at the MTD (TR_MTD)
  n_patient_each_dose = colMeans(N)
  n_patient_mean = sum(N) / n_trial
  if(is_all_ovedose){
    TR_MTD = 0
  }else{
    TR_MTD = n_patient_each_dose[target_dose] / n_patient_mean * 100
  }
  
  
  ## 2. safty
  # 2.1 percentage of patients treated above the MTD (TR_aboveMTD)
  TR_aboveMTD = sum(n_patient_each_dose[tox_prob > (target_prob + margin_right)] / n_patient_mean * 100)
  
  # 2.2 risk of overdosing 
  overdosing60 = mean(rowSums(N[, tox_prob > (target_prob + margin_right), drop = FALSE]) > 0.6 * n_patient) * 100
  overdosing80 = mean(rowSums(N[, tox_prob > (target_prob + margin_right), drop = FALSE]) > 0.8 * n_patient) * 100
  
  # duration time 
  duration_mean = mean(duration)
  #------------------------ end performance metrics -----------------------#
  
  out = list(PCS = PCS, TR_MTD = TR_MTD, TR_aboveMTD = TR_aboveMTD, 
             overdosing60 = overdosing60, overdosing80 = overdosing80,
             dose_select = dose_select, select_percent = select_percent, 
             n_patient_mean = n_patient_mean, duration_mean = duration_mean,
             Y = Y, N = N, dose_Paths = dose_Paths, DLT_Paths = DLT_Paths)
  return(out)
}







# target_prob = 0.3
# tox_prob = c(0.01,0.12,0.30,0.41,0.55)
# n_cohort = 10
# cohort_size = 3
# symmetric = FALSE
# k_left = 0.2
# k_right = 0.8
# shared = TRUE # incorporate the keybooard design
# dose_set = 1:length(tox_prob)
# n_earlystop = 1000
# start_dose = 1
# margin_left = 0.05
# margin_right = 0.05
# cutoff_elimin = 0.95
# extraSafe = FALSE
# offset = 0.05
# tau = 3
# accrual = 2
# alpha = 0.5
# piecewise = FALSE
# prior_p = rep(1/3,3)
# dist_DLT = "weibull"
# dist_enter = "exp"
# n_trial = 1000
# seed = 6
