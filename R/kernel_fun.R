#' Asymmetric/symmetric kernel with interpretable neighbor borrowing strengths
#'
#' @description
#' `kernel_fun()` computes kernel similarities between a target dose value `dose` and a set of dose
#' locations `dose_set` (typically on a standardized scale, e.g., `[0,1]`).
#'
#' Unlike the conventional parameterization that requires users to specify decay parameters
#' `theta` directly, this kernel is parameterized by **interpretable neighbor kernel values**:
#' `k_left` and `k_right`, which represent the desired kernel similarity at a reference distance
#' `ref_gap` on the left and right sides, respectively. The decay parameters are inferred internally via
#' \deqn{\theta = -\log(k) / \text{ref\_gap}^2.}
#'
#' @details
#' The kernel takes the form \eqn{k(d,d_s)=\exp\{-\theta(d-d_s)^2\}}.
#'
#' \itemize{
#'   \item **Asymmetric mode** (`symmetric = FALSE`):
#'     \eqn{\theta_1} is used for dose locations strictly below `dose` (`dose_set < dose`),
#'     and \eqn{\theta_2} is used for dose locations at or above `dose` (`dose_set >= dose`),
#'     where
#'     \deqn{\theta_1 = -\log(k_{\text{left}})/\text{ref\_gap}^2,\qquad
#'           \theta_2 = -\log(k_{\text{right}})/\text{ref\_gap}^2.}
#'     If `k_right > k_left`, then \eqn{\theta_2 < \theta_1}, meaning the kernel decays more slowly
#'     on the right side (stronger borrowing from higher doses), matching the typical safety-oriented
#'     SKBD/BKP borrowing direction.
#'
#'   \item **Symmetric mode** (`symmetric = TRUE`):
#'     A single \eqn{\theta} is inferred using `k_right` (or its default) at distance `ref_gap`.
#' }
#'
#' The default `ref_gap` is set to the minimum adjacent spacing in `dose_set`, i.e.,
#' `min(diff(dose_set))`. This makes `k_left`/`k_right` interpretable as the desired borrowing strength
#' at the smallest dose spacing on the current grid.
#'
#' @param dose Numeric scalar. Target dose value (standardized dose level, not an index).
#' @param dose_set Numeric vector. Dose locations at which to evaluate similarity (standardized doses).
#'   Typically strictly increasing.
#' @param symmetric Logical. If `TRUE`, use a single symmetric decay; if `FALSE` (default), use an
#'   asymmetric decay with separate left/right decay rates.
#' @param k_left Numeric scalar in `(0,1)`. Desired kernel value at distance `ref_gap` for dose locations
#'   strictly below `dose` (`dose_set < dose`). Default is `0.2` when omitted.
#' @param k_right Numeric scalar in `(0,1)`. Desired kernel value at distance `ref_gap` for dose locations
#'   at or above `dose` (`dose_set >= dose`). Default is `0.8` when omitted.
#' @param ref_gap Numeric scalar. Reference distance used to interpret `k_left` and `k_right`.
#'   If `NULL`, defaults to `min(diff(dose_set))`. Must be positive and finite.
#'
#' @return Numeric vector. Kernel similarities `k` of length `length(dose_set)`, with entries in `(0,1]`.
#'
#' @examples
#' dose_set = c(0.25, 0.50, 0.75)
#'
#' # Asymmetric borrowing: weaker from the left (0.2), stronger from the right (0.8)
#' kernel_fun(dose = 0.50, dose_set = dose_set, symmetric = FALSE, k_left = 0.2, k_right = 0.8)
#'
#' # Symmetric borrowing using k_right as the reference similarity
#' kernel_fun(dose = 0.50, dose_set = dose_set, symmetric = TRUE, k_right = 0.5)
#'
#' # Fix the reference distance explicitly (instead of using min(diff(dose_set)))
#' kernel_fun(dose = 0.50, dose_set = dose_set, symmetric = FALSE,
#'            k_left = 0.2, k_right = 0.8, ref_gap = 0.25)
#'
#' @noRd
kernel_fun <- function(
    dose, dose_set,
    symmetric = FALSE,
    k_left = NULL, k_right = NULL,
    ref_gap = NULL
){
  # dose: scalar standardized dose (not index)
  # dose_set: vector of standardized dose locations (not indices)
  
  # ---- Validate neighbor kernel values (defaults if missing) ----
  # Require at least one of k_left/k_right; fill missing side with defaults
  if (is.null(k_left) && is.null(k_right)) {
    stop("please provide `k_left` and/or `k_right` (both in (0,1)).")
  }
  if (is.null(k_left))  k_left  = 0.2
  if (is.null(k_right)) k_right = 0.8
  
  # Ensure both neighbor kernel values are scalars in (0,1)
  if (!is.numeric(k_left)  || length(k_left)  != 1 || !(0 < k_left  && k_left  < 1)) {
    stop("`k_left` must be a scalar in (0,1).")
  }
  if (!is.numeric(k_right) || length(k_right) != 1 || !(0 < k_right && k_right < 1)) {
    stop("`k_right` must be a scalar in (0,1).")
  }
  
  # ---- Choose reference distance (gap) ----
  # Default: smallest adjacent spacing in dose_set (interprets k_left/k_right as neighbor borrowing)
  if (is.null(ref_gap)) {
    if (length(dose_set) < 2) stop("need at least two distinct dose points to infer ref_gap.")
    ref_gap = min(diff(dose_set))
  }
  if (!is.numeric(ref_gap) || length(ref_gap) != 1 || !is.finite(ref_gap) || ref_gap <= 0) {
    stop("`ref_gap` must be a positive finite scalar.")
  }
  
  # ---- Infer decay parameter(s) theta internally from neighbor kernel values ----
  # exp(-theta * ref_gap^2) = k  =>  theta = -log(k) / ref_gap^2
  if (symmetric) {
    
    # Symmetric kernel: use one theta based on k_right (default 0.8 if omitted)
    k0 = k_right
    theta = -log(k0) / (ref_gap^2)
    
    # Compute symmetric Gaussian kernel values
    k = exp(-theta * (dose - dose_set)^2)
    
  } else {
    
    # Asymmetric kernel: different decay on left and right sides of `dose`
    theta1 = -log(k_left)  / (ref_gap^2)  # left side (dose_set < dose): typically steeper decay
    theta2 = -log(k_right) / (ref_gap^2)  # right side (dose_set >= dose): typically slower decay
    
    # Assign side-specific theta by comparing dose_set locations to dose
    n_dose = length(dose_set)
    n_dose_left = sum(dose_set < dose)
    theta_set = c(rep(theta1, n_dose_left), rep(theta2, n_dose - n_dose_left))
    
    # Compute asymmetric Gaussian kernel values
    k = exp(-theta_set * (dose - dose_set)^2)
  }
  
  return(k)
}
