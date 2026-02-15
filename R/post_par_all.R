#' Compute BKP/SKBD pseudo-posterior Beta parameters for all doses on a working grid
#'
#' @description
#' `post_par_all()` computes the BKP/SKBD pseudo-posterior Beta parameters
#' \eqn{\alpha_j, \beta_j} for **every** dose level \eqn{d_j} on the current working dose grid.
#' The update borrows information only from dose levels with accrued data (i.e., `n_treated > 0`)
#' via kernel-normalized weights.
#'
#' @details
#' Let \eqn{S=\{s: n_s > 0\}} denote the set of observed (treated) dose indices.
#' For each target dose \eqn{d_j} and each observed dose \eqn{d_s, s\in S}, define kernel
#' similarity \eqn{k(d_j, d_s)} and normalize across \eqn{S} to obtain weights
#' \deqn{w_{js} = \frac{k(d_j,d_s)}{\sum_{t\in S} k(d_j,d_t)}.}
#' The pseudo-posterior Beta parameters are then
#' \deqn{\alpha_j = \alpha_0 + \sum_{s\in S} w_{js} y_s,\qquad
#'       \beta_j  = \beta_0  + \sum_{s\in S} w_{js} (m_s-y_s),}
#' where \eqn{y_s} is the number of DLTs and \eqn{m_s} is the number treated at dose \eqn{d_s}.
#'
#' **Numerical stabilizations**
#' \itemize{
#'   \item Each row of the kernel matrix is rescaled by its row maximum `km = max(K[j, ])`
#'         (if `km` is non-finite or non-positive, it is set to 1). This reduces underflow
#'         when kernel values are extremely small.
#'   \item The row sum `den` used for normalization is lower-bounded by `1e-12` to avoid
#'         division by zero.
#' }
#'
#' @param n_dlt Numeric/integer vector. Number of observed DLTs at each working dose.
#'   Must have the same length as `dose_set`.
#' @param n_treated Numeric/integer vector. Number treated at each working dose.
#'   Must have the same length as `dose_set`.
#' @param dose_set Numeric vector. Working dose grid (often standardized to `[0, 1]`), assumed sorted.
#' @param pri_alpha Numeric scalar. Prior alpha \eqn{\alpha_0} for Beta pseudo-posterior.
#' @param pri_beta Numeric scalar. Prior beta \eqn{\beta_0} for Beta pseudo-posterior.
#' @param symmetric Logical. Passed to `kernel()` to indicate symmetric vs asymmetric kernel mode.
#' @param theta Numeric. Passed to `kernel()` as kernel hyperparameter(s).
#'   Length 1 for symmetric, length 2 for asymmetric (per your `kernel()` definition).
#'
#' @return A list with two numeric vectors (both of length `length(dose_set)`):
#' \describe{
#'   \item{post_alpha}{Pseudo-posterior alpha parameters at each working dose.}
#'   \item{post_beta}{Pseudo-posterior beta parameters at each working dose.}
#' }
#'
#' @examples
#' # Example with 3 doses and data accrued at doses 1 and 2
#' dose_set = c(0.25, 0.50, 0.75)
#' n_treated = c(3, 6, 0)
#' n_dlt = c(0, 2, 0)
#'
#' out = post_par_all(
#'   n_dlt = n_dlt,
#'   n_treated = n_treated,
#'   dose_set = dose_set,
#'   pri_alpha = 0.5,
#'   pri_beta = 0.5,
#'   symmetric = FALSE,
#'   theta = c(10, 5)
#' )
#' out$post_alpha
#' out$post_beta
#'
#' @export
post_par_all = function(
    n_dlt, n_treated, dose_set,
    pri_alpha, pri_beta, symmetric, theta
) {
  
  # ---- Basic setup ----
  n_dose = length(dose_set)
  obs_idx = which(n_treated > 0)
  
  # Default: prior only (if no data accrued yet)
  post_alpha = rep(pri_alpha, n_dose)
  post_beta  = rep(pri_beta,  n_dose)
  
  if (length(obs_idx) == 0) return(list(post_alpha = post_alpha, post_beta = post_beta))
  
  # ---- Extract observed-dose quantities once ----
  dose_obs = dose_set[obs_idx]
  y_obs    = n_dlt[obs_idx]
  m_obs    = n_treated[obs_idx]
  z_obs    = m_obs - y_obs
  
  # ---- Build kernel matrix K and normalized weight matrix W ----
  K = W = matrix(0, nrow = n_dose, ncol = length(obs_idx))
  
  for (j in 1:n_dose) {
    
    # Kernel similarities between target dose_set[j] and all observed doses
    K[j, ] = kernel(dose_set[j], dose_obs, symmetric, theta)
    
    # Row-wise rescaling to reduce underflow when values are extremely small
    km = max(K[j, ])
    if (!is.finite(km) || km <= 0) km = 1
    K[j, ] = K[j, ] / km
    
    # Normalize to weights; protect against zero/invalid denominators
    den = sum(K[j, ])
    if (!is.finite(den) || den <= 1e-12) den = 1e-12
    W[j, ] = K[j, ] / den
  }
  
  # ---- Weighted pseudo-count update to Beta parameters ----
  post_alpha = pri_alpha + drop(W %*% y_obs)
  post_beta  = pri_beta  + drop(W %*% z_obs)
  
  return(list(post_alpha = post_alpha, post_beta = post_beta))
}
