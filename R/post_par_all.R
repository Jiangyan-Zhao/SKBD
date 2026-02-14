#' @title BKP posterior update for ALL doses
#' 
#' @description
#' A short description...
#' 
#' 
#' @examples
#' n_dlt = c(0, 0, 1, 2, 0)
#' n_treated = c(3, 3, 6, 3, 0)
#' dose_set = seq(0, 1, length.out = 5)
#' pri_alpha = 1
#' pri_beta = 1 
#' symmetric = TRUE
#' theta = 0.7
#' post_par_all(n_dlt, n_treated, dose_set, pri_alpha, pri_beta, symmetric, theta)


post_par_all = function(
    n_dlt, n_treated, dose_set, 
    pri_alpha, pri_beta, symmetric, theta
) {
  n_dose = length(dose_set)
  obs_idx = which(n_treated > 0)
  
  post_alpha = rep(pri_alpha, n_dose)
  post_beta  = rep(pri_beta,  n_dose)
  
  if (length(obs_idx) == 0) return(list(post_alpha = post_alpha, post_beta = post_beta))
  
  dose_obs = dose_set[obs_idx]
  y_obs    = n_dlt[obs_idx]
  m_obs    = n_treated[obs_idx]
  z_obs    = m_obs - y_obs
  
  K = W = matrix(0, nrow = n_dose, ncol = length(obs_idx))
  for (j in 1:n_dose) {
    K[j, ] = kernel(dose_set[j], dose_obs, symmetric, theta)
    km = max(K[j, ])
    if (!is.finite(km) || km <= 0) km = 1
    K[j, ] = K[j, ] / km
    den = sum(K[j, ])
    if (!is.finite(den) || den <= 1e-12) den = 1e-12
    W[j, ] = K[j, ] / den
  }
  post_alpha = pri_alpha + drop(W %*% y_obs)
  post_beta  = pri_beta  + drop(W %*% z_obs)
  list(post_alpha = post_alpha, post_beta = post_beta)
}
