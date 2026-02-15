#' Assign "true" toxicity probability for an inserted dose (simulation-only)
#'
#' @description
#' `ins_true_prob()` is a **simulation-only** helper used in operating-characteristics (OC)
#' evaluation of adaptive dose insertion. Given the current *true* dose-toxicity vector `p`
#' on the working grid, it assigns a "true toxicity probability" to the newly inserted dose
#' according to simple ADM-style rules.
#'
#' In a real trial, the true toxicity probabilities are unknown, so this function is **not**
#' used in practice.
#'
#' @details
#' The function supports three insertion locations encoded by `insert_code`:
#' \itemize{
#'   \item **Between** \eqn{(i, i+1)}: `insert_code = i` with `i >= 1`.
#'     If \eqn{p_i < \phi < p_{i+1}}, set \eqn{p_{\text{new}} = \phi}; otherwise set
#'     \eqn{p_{\text{new}} = (p_i + p_{i+1})/2}.
#'   \item **Below minimum**: `insert_code = LOC_BELOW_MIN` (default `-1L`).
#'     If \eqn{\phi < p_d}, set \eqn{p_{\text{new}} = \phi}; otherwise set \eqn{p_{\text{new}} = p_d/2}.
#'   \item **Above maximum**: `insert_code = LOC_ABOVE_MAX` (default `-2L`).
#'     If \eqn{\phi > p_d}, set \eqn{p_{\text{new}} = \phi}; otherwise set \eqn{p_{\text{new}} = p_d + 0.1}.
#' }
#'
#' For numerical stability in simulations, the returned probability is clamped to
#' \eqn{(\epsilon, 1-\epsilon)} with \eqn{\epsilon = 10^{-6}}.
#'
#' @param insert_code Integer. Insertion location code:
#'   `i >= 1` for inserting between doses `(i, i+1)`;
#'   `LOC_BELOW_MIN` for inserting below the minimum dose;
#'   `LOC_ABOVE_MAX` for inserting above the maximum dose.
#' @param d Integer. Current dose index at which insertion is triggered (used for boundary cases).
#'   Must satisfy `1 <= d <= length(p)`.
#' @param p Numeric vector. Current *true* toxicity probabilities aligned with the working dose grid
#'   (simulation-only). Must have values in `[0, 1]`.
#' @param target Numeric scalar. Target toxicity rate \eqn{\phi} in `(0, 1)`.
#' @param LOC_BELOW_MIN Integer. Code representing insertion below minimum dose. Default `-1L`.
#' @param LOC_ABOVE_MAX Integer. Code representing insertion above maximum dose. Default `-2L`.
#'
#' @return Numeric scalar. The assigned true toxicity probability for the inserted dose, clamped to
#' `(\eqn{\epsilon}, 1-\eqn{\epsilon})` with `\eqn{\epsilon}=1e-6`.
#'
#' @examples
#' # Between (2, 3)
#' p = c(0.05, 0.15, 0.30, 0.45)
#' ins_true_prob(insert_code = 2L, d = 2L, p = p, target = 0.25)
#'
#' # Below minimum
#' ins_true_prob(insert_code = -1L, d = 1L, p = p, target = 0.25)
#'
#' # Above maximum
#' ins_true_prob(insert_code = -2L, d = length(p), p = p, target = 0.25)
#'
#' @export

ins_true_prob = function(
    insert_code, d, p, target,
    LOC_BELOW_MIN = -1L, LOC_ABOVE_MAX = -2L
) {
  
  J = length(p)
  if (J < 1) stop("p must have positive length.")
  if (!(0 < target && target < 1)) stop("target must be in (0, 1).")
  if (d < 1 || d > J) stop("d out of range.")
  
  p_new = NA_real_
  
  # ----------------------------------------------------------
  # Case 1) Insert BETWEEN (i, i+1)
  #   if p[i] < target < p[i+1] then p_new = target
  #   else p_new = (p[i] + p[i+1]) / 2
  # ----------------------------------------------------------
  if (!is.na(insert_code) && insert_code >= 1) {
    
    i = as.integer(insert_code)
    if (i < 1 || i >= J) stop("insert_code (between) out of range.")
    
    p_new = ifelse (p[i] < target && target < p[i + 1], target, (p[i] + p[i + 1]) / 2)
    
    # ----------------------------------------------------------
    # Case 2) Insert BELOW minimum:
    #   if target < p[d] then p_new = target
    #   else p_new = p[d] / 2
    #
    # In your Scheme A logic, this case should occur when d == 1 (lowest dose),
    # but we keep it general and use p[d] exactly as ADM does.
    # ----------------------------------------------------------
  } else if (!is.na(insert_code) && insert_code == LOC_BELOW_MIN) {
    
    p_new = ifelse (target < p[d], target, p[d] / 2)
    
    # ----------------------------------------------------------
    # Case 3) Insert ABOVE maximum: 
    #   if target > p[d] then p_new = target
    #   else p_new = p[d] + 0.1
    #
    # This is how ADM generates a higher true tox probability for extrapolation upward.
    # ----------------------------------------------------------
  } else if (!is.na(insert_code) && insert_code == LOC_ABOVE_MAX) {
    
    p_new = ifelse (target > p[d], target, p[d] + 0.1)
    
  } else {
    stop("unsupported insert_code.")
  }
  
  # Clamp into (0,1) for safety in simulation
  eps_prob = 1e-6
  p_new = max(min(p_new, 1 - eps_prob), eps_prob)
  
  return(p_new)
}
