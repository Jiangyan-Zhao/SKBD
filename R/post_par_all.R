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
#' @param symmetric Logical. Passed to `kernel_fun()` to indicate symmetric vs asymmetric kernel mode.
#' @param k_left Numeric scalar in `(0,1)`. Passed to `kernel_fun()` as left-side
#'   neighbor borrowing strength (used when `dose_set < dose`).
#' @param k_right Numeric scalar in `(0,1)`. Passed to `kernel_fun()` as right-side
#'   neighbor borrowing strength (used when `dose_set >= dose`, and as the symmetric
#'   borrowing strength when `symmetric = TRUE`).
#' @param ref_gap Optional positive scalar. Passed to `kernel_fun()` as the reference
#'   spacing used to infer decay from `k_left`/`k_right`.
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
#'   k_left = 0.2,
#'   k_right = 0.8
#' )
#' out$post_alpha
#' out$post_beta
#'
#' @export
post_par_all = function(
    n_dlt, n_treated, dose_set,
    pri_alpha, pri_beta, symmetric,
    k_left = 0.2, k_right = 0.8, ref_gap = NULL, ker_scale = TRUE
) {
  
  n_dose = length(dose_set)
  
  # Default: prior only (if no data accrued yet)
  post_alpha = rep(pri_alpha, n_dose)
  post_beta  = rep(pri_beta,  n_dose)
  
  obs_idx = which(n_treated > 0)
  
  if (length(obs_idx) == 0) return(list(post_alpha = post_alpha, post_beta = post_beta))
  
  # ---- Build kernel matrix K and normalized weight matrix W ----
  for (d in obs_idx) {
    
    # Kernel similarities between target dose_set[j] and all observed doses
    ker_vals = kernel_fun(
      dose_set[d], dose_set,
      symmetric = symmetric,
      k_left = k_left,
      k_right = k_right,
      ref_gap = ref_gap
    )
    if(ker_scale){
      weight = ker_vals[obs_idx] / sum(ker_vals[obs_idx])
    }else{
      weight = ker_vals[obs_idx]
    }
    
    post_alpha[d] = pri_alpha + sum(weight * n_dlt[obs_idx])
    post_beta[d] = pri_beta + sum(weight * (n_treated[obs_idx] - n_dlt[obs_idx]))
  }
  
  return(list(post_alpha = post_alpha, post_beta = post_beta))
}
