#' @title Check Whether Dose Insertion Should Be Triggered
#'
#' @description
#' Determine whether a new dose should be inserted around the current dose level
#' based on posterior interval probabilities under the insertion rule.
#'
#' @details
#' This function implements the insertion trigger by comparing posterior
#' probabilities of being below, within, or above the target key at each
#' currently available dose level.
#'
#' For dose \eqn{j}, define
#' \deqn{
#' \Pr(\pi_j \le \mathrm{key\_L} \mid \mathcal D),
#' }
#' \deqn{
#' \Pr(\pi_j > \mathrm{key\_U} \mid \mathcal D),
#' }
#' and
#' \deqn{
#' \Pr(\mathrm{key\_L} < \pi_j \le \mathrm{key\_U} \mid \mathcal D),
#' }
#' where \eqn{\pi_j} is the toxicity probability at dose \eqn{j}.
#'
#' To stabilize the insertion decision under finite samples, the raw posterior
#' probabilities are adjusted by isotonic regression using the pool-adjacent-
#' violators algorithm (PAVA): the probability of being below the target key is
#' enforced to be non-increasing across doses, and the probability of being
#' above the target key is enforced to be non-decreasing across doses.
#'
#' The function then checks whether the evidence around the current dose
#' supports inserting a new dose between two adjacent doses or outside the
#' current dose range. Interior insertion is triggered when the left dose has
#' sufficiently large posterior probability of being below the target key and
#' the right dose has sufficiently large posterior probability of being above
#' the target key. Boundary insertion is treated similarly at the lower and
#' upper ends of the dose range.
#'
#' @param j Positive integer. Index of the current dose level.
#' @param post_alpha Numeric vector of posterior alpha parameters for the Beta
#' distributions at all current dose levels.
#' @param post_beta Numeric vector of posterior beta parameters for the Beta
#' distributions at all current dose levels. Must have the same length as
#' \code{post_alpha}.
#' @param key_L Scalar in \eqn{(0,1)}. Lower boundary of the target key.
#' @param key_U Scalar in \eqn{(0,1)}. Upper boundary of the target key, with
#' \code{key_L < key_U}.
#' @param C1 Scalar in \eqn{[0,1]}. Threshold for the posterior probability of
#' being below the target key.
#' @param C2 Scalar in \eqn{[0,1]}. Threshold for the posterior probability of
#' being above the target key.
#' @param n_treated Integer vector giving the number of treated patients at each
#' current dose level. Must have the same length as \code{post_alpha}.
#' @param LOC_BELOW_MIN Integer code used to indicate insertion below the
#' minimum current dose. Default is \code{-1L}.
#' @param LOC_ABOVE_MAX Integer code used to indicate insertion above the
#' maximum current dose. Default is \code{-2L}.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{need_insert}}{Logical. Whether dose insertion is triggered.}
#'   \item{\code{insert_code}}{Integer insertion code. \code{NA} means no
#'   insertion; \code{LOC_BELOW_MIN} means insertion below the minimum current
#'   dose; \code{LOC_ABOVE_MAX} means insertion above the maximum current dose;
#'   a positive integer \code{i} means insertion between doses \code{i} and
#'   \code{i + 1}.}
#'   \item{\code{LOC_BELOW_MIN}}{The code used for lower-boundary insertion.}
#'   \item{\code{LOC_ABOVE_MAX}}{The code used for upper-boundary insertion.}
#' }
#'
#' @examples
#' post_alpha <- c(1.5, 2.0, 3.5, 5.0)
#' post_beta  <- c(8.5, 7.0, 5.5, 4.0)
#' n_treated  <- c(3, 3, 6, 6)
#'
#' insert_check(
#'   j = 2,
#'   post_alpha = post_alpha,
#'   post_beta = post_beta,
#'   key_L = 0.15,
#'   key_U = 0.25,
#'   C1 = 0.6,
#'   C2 = 0.6,
#'   n_treated = n_treated
#' )
#'
#' @importFrom stats pbeta
#' @noRd

insert_check = function(
    j, post_alpha, post_beta, key_L, key_U, C1, C2, n_treated,
    LOC_BELOW_MIN = -1L, LOC_ABOVE_MAX = -2L
) {
  
  J = length(post_alpha)
  
  if (length(post_beta) != J) stop("post_alpha and post_beta must have same length.")
  if (length(n_treated) != J) stop("n_treated must have same length as post_alpha.")
  if (j < 1 || j > J) stop("j must be in 1..J.")
  if (!(0 < key_L && key_L < key_U && key_U < 1)) stop("require 0 < key_L < key_U < 1.")
  if (!(0 <= C1 && C1 <= 1 && 0 <= C2 && C2 <= 1)) stop("C1,C2 must be in [0,1].")
  
  # -----------------------------
  # Step A: raw region probabilities at each discrete dose
  # -----------------------------
  prob_under = pbeta(key_L, post_alpha, post_beta)
  prob_over  = 1 - pbeta(key_U, post_alpha, post_beta)
  prob_target = 1 - prob_under - prob_over
  
  # -----------------------------
  # Step B: ADM-style PAVA monotonic adjustment (optional)
  # -----------------------------
  
  prob_under_adj = rev(pava(rev(prob_under)))   # enforce non-increasing
  prob_over_adj  = pava(prob_over)              # enforce non-decreasing
  prob_target_adj = pmax(0, 1 - prob_under_adj - prob_over_adj)

  
  # -----------------------------
  # Step C: Eq.(5) trigger using adjusted probabilities
  # -----------------------------
  insert_code = NA_integer_
  
  # Left interval: (d_{j-1}, d_j)  <=> l = j
  if (j == 1) {
    # Boundary (d0, d1): pi_0 = 0 => Pr(pi_0 < key_L) = 1, so only need right bracket evidence.
    left_ok = (n_treated[1] > 0) && (prob_over_adj[1] > C2)
    if (left_ok) insert_code = LOC_BELOW_MIN
  } else {
    left_ok = (n_treated[j - 1] > 0) && (n_treated[j] > 0) &&
      (prob_under_adj[j - 1] > C1) && (prob_over_adj[j] > C2)
    if (left_ok) insert_code = as.integer(j - 1)
  }
  
  # Right interval: (d_j, d_{j+1}) <=> l = j+1
  if (is.na(insert_code)) { # It is impossible to satisfy simultaneously.
    if (j == J) {
      # Boundary (dJ, dJ+1): pi_{J+1} = 1 => Pr(pi_{J+1} > key_U) = 1, so only need left bracket evidence.
      right_ok = (n_treated[J] > 0) && (prob_under_adj[J] > C1)
      if (right_ok) insert_code = LOC_ABOVE_MAX
    } else {
      right_ok = (n_treated[j] > 0) && (n_treated[j + 1] > 0) &&
        (prob_under_adj[j] > C1) && (prob_over_adj[j + 1] > C2)
      if (right_ok) insert_code = as.integer(j)
    }
  }
  
  need_insert = !is.na(insert_code)
  
  out = list(
    need_insert = need_insert,
    insert_code = insert_code,
    LOC_BELOW_MIN = LOC_BELOW_MIN,
    LOC_ABOVE_MAX = LOC_ABOVE_MAX
  )

  return(out)
}
