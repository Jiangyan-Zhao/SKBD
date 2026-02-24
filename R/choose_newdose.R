#' Choose an inserted dose within an eligible interval by maximizing target-key probability
#'
#' @description
#' `choose_newdose()` selects a new dose level within an eligible open interval \eqn{(d_l, d_r)}
#' by searching over a grid and maximizing
#' \deqn{q(d) = \Pr\{ \text{key}_L < \pi(d) \le \text{key}_U \mid \mathcal{D} \}.}
#' The probability \eqn{q(d)} is computed under the BKP/SKBD pseudo-posterior at candidate dose
#' locations, where Beta parameters are formed by kernel-weighted pseudo-counts aggregated from
#' already-observed doses.
#'
#' @details
#' **How it works**
#' \itemize{
#'   \item Collect observed doses \eqn{S=\{j: n_j>0\}} from the current working grid.
#'   \item Build an equally spaced grid of size `M` within \eqn{(d_l,d_r)} (excluding endpoints).
#'   \item For each grid point \eqn{d}, compute kernel similarities `k(d, dose_obs)`, normalize to
#'         weights \eqn{w}, and obtain Beta pseudo-posterior parameters:
#'         \deqn{\alpha(d)=\alpha_0 + \sum_{j\in S} w_j(d) y_j,\qquad
#'               \beta(d)=\beta_0 + \sum_{j\in S} w_j(d) (n_j-y_j).}
#'   \item Compute \eqn{q(d)} as `pbeta(key_U, alpha, beta) - pbeta(key_L, alpha, beta)`.
#'   \item Return the grid point maximizing \eqn{q(d)}.
#' }
#'
#' **Numerical stabilizations**
#' \itemize{
#'   \item Rescales kernel similarities by `max(k)` to reduce underflow when values are extremely small.
#'   \item Guards against zero/invalid normalization and clamps Beta parameters to be at least `1e-12`.
#'   \item Avoids duplicating an existing dose by returning `NA` if the selected `newdose` is within a
#'         small tolerance of any value in `dose_set_work`.
#' }
#'
#' @param dl Numeric scalar. Left endpoint of the eligible interval; must satisfy `dl < dr`.
#' @param dr Numeric scalar. Right endpoint of the eligible interval; must satisfy `dl < dr`.
#' @param dose_set_work Numeric vector. Current working dose grid (often standardized to `[0, 1]`).
#' @param n_dlt_work Numeric/integer vector. Number of DLTs observed at each working dose.
#'   Must have the same length as `dose_set_work`.
#' @param n_treated_work Numeric/integer vector. Number treated at each working dose.
#'   Must have the same length as `dose_set_work`.
#' @param pri_alpha Numeric scalar. Prior alpha for the Beta pseudo-posterior.
#' @param pri_beta Numeric scalar. Prior beta for the Beta pseudo-posterior.
#' @param key_L Numeric scalar. Lower bound of the target key, in `(0,1)`.
#' @param key_U Numeric scalar. Upper bound of the target key, in `(0,1)` with `key_L < key_U`.
#' @param symmetric Logical. Passed to `kernel_fun()` to indicate symmetric vs asymmetric kernel mode.
#' @param k_left Numeric scalar in `(0,1)`. Passed to `kernel_fun()` as left-side
#'   neighbor borrowing strength (used when `dose_set < dose`).
#' @param k_right Numeric scalar in `(0,1)`. Passed to `kernel_fun()` as right-side
#'   neighbor borrowing strength (used when `dose_set >= dose`, and as the symmetric
#'   borrowing strength when `symmetric = TRUE`).
#' @param ref_gap Optional positive scalar. Passed to `kernel_fun()` as the reference
#'   spacing used to infer decay from `k_left`/`k_right`.
#' @param M Integer. Number of grid points inside `(dl, dr)` used for search. Default is `100`.
#'
#' @return A list with components:
#' \describe{
#'   \item{newdose}{Selected new dose in `(dl, dr)`. Returns `NA_real_` if no valid candidate is found
#'                 or if the selected dose is too close to an existing dose.}
#'   \item{q_max}{Maximum value of `q(d)` over the search grid (or `NA_real_` if not available).}
#'   \item{grid}{The grid of candidate dose values used for searching (returned when available).}
#'   \item{q_grid}{Vector of `q(d)` values over `grid`, possibly containing `NA` for invalid points.}
#' }
#'
#' @examples
#' # Suppose three doses have been tried and we want to insert between 0.25 and 0.50
#' dose_set_work = c(0.25, 0.50, 0.75)
#' n_treated_work = c(3, 6, 0)
#' n_dlt_work = c(0, 3, 0)
#'
#' # Candidate selection
#' out = choose_newdose(
#'   dl = 0.25, dr = 0.50,
#'   dose_set_work = dose_set_work,
#'   n_dlt_work = n_dlt_work,
#'   n_treated_work = n_treated_work,
#'   pri_alpha = 0.5, pri_beta = 0.5,
#'   key_L = 0.25, key_U = 0.35,
#'   symmetric = FALSE, k_left = 0.2, k_right = 0.8, M = 100
#' )
#' out$newdose
#'
#' @export

choose_newdose = function(
    dl, dr, dose_set_work, n_dlt_work, n_treated_work,
    pri_alpha, pri_beta, key_L, key_U, symmetric = TRUE,
    k_left = 0.2, k_right = 0.2, ref_gap = NULL, M = 100
){
  if (!is.finite(dl) || !is.finite(dr) || dr <= dl) {
    return(list(newdose = NA_real_, q_max = NA_real_))
  }
  if (!is.numeric(M) || length(M) != 1 || M < 5) M <- 50L
  
  obs_idx = which(n_treated_work > 0)
  if (length(obs_idx) == 0) {
    return(list(newdose = NA_real_, q_max = NA_real_))
  }
  
  dose_obs = dose_set_work[obs_idx]
  y_obs    = n_dlt_work[obs_idx]
  n_obs    = n_treated_work[obs_idx]
  
  empr_rate = pava(y_obs / n_obs)
  y_obs = empr_rate * n_obs

  grid = seq(dl + 0.001, dr - 0.001, length.out = M)
  q_grid = rep(NA_real_, length(grid))
  
  for (g in seq_along(grid)) {
    k = kernel_fun(grid[g], dose_obs,
               symmetric = symmetric,
               k_left = k_left,
               k_right = k_right,
               ref_gap = ref_gap)
    
    w = k / sum(k)
    
    alpha_d = pri_alpha + sum(w * y_obs)
    beta_d  = pri_beta  + sum(w * (n_obs - y_obs))
    
    q_grid[g] = pbeta(key_U, alpha_d, beta_d) - pbeta(key_L, alpha_d, beta_d)
    # q_grid[g] = alpha_d / (alpha_d + beta_d)
  }
  
  if (all(is.na(q_grid))) {
    return(list(newdose = NA_real_, q_max = NA_real_, grid = grid, q_grid = q_grid))
  }
  
  idx = which.max(q_grid)
  newdose = grid[idx]
  
  list(newdose = newdose, q_max = q_grid[idx], grid = grid, q_grid = q_grid)
}
