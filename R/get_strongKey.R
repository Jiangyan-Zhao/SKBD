#' @title Get the strongest key
#'
#' @description
#' Given a Beta posterior \code{Beta(alpha_Post, beta_Post)} for the toxicity probability at a
#' dose, compute the posterior probability mass in each keyboard key interval defined by
#' \code{keys} and return the index of the key with the largest (adjusted) probability.
#'
#' @details
#' The input \code{keys} is an increasing vector of cutoffs in \eqn{[0,1]} that defines
#' \eqn{K = length(keys)-1} adjacent key intervals \eqn{(keys[i], keys[i+1])}.
#'
#' This function applies a boundary-width adjustment for the first and last keys (which may be
#' wider than the target key due to truncation at 0 and 1). Specifically, for \code{i = 1} or
#' \code{i = K}, the probability mass is multiplied by
#' \deqn{\frac{\epsilon_1+\epsilon_2}{keys[i+1]-keys[i]},}
#' where \eqn{\epsilon_1 = margin_left} and \eqn{\epsilon_2 = margin_right}.
#'
#' @param alpha_Post Positive scalar. Posterior shape parameter \eqn{\alpha} of the Beta distribution.
#' @param beta_Post Positive scalar. Posterior shape parameter \eqn{\beta} of the Beta distribution.
#' @param keys Numeric vector of increasing cutoffs in \eqn{[0,1]} defining the keyboard keys,
#' typically produced by \code{\link{get_Key}}.
#' @param margin_left Nonnegative scalar. Left half-width \eqn{\epsilon_1} of the target key.
#' @param margin_right Nonnegative scalar. Right half-width \eqn{\epsilon_2} of the target key.
#'
#' @return
#' An integer in \code{1:(length(keys)-1)} giving the index of the strongest key (the key with the
#' maximum adjusted posterior probability mass).
#'
#' @examples
#' keys <- get_Key(0.30, 0.05, 0.05)
#' get_strongKey(alpha_Post = 3, beta_Post = 7, keys = keys)
#'
#' @noRd
get_strongKey <- function(
    alpha_Post, beta_Post, keys, 
    margin_left = 0.05, margin_right = 0.05
) {
  n_keys = length(keys) - 1  
  key_prob = rep(0, n_keys)
  for (i in 1:n_keys) {
    comp = 1
    if (i == 1 || i == n_keys) {
      comp = (margin_left + margin_right) / (keys[i+1] - keys[i])
    }
    key_prob[i] = pbeta(keys[i+1], alpha_Post, beta_Post) - 
      pbeta(keys[i], alpha_Post, beta_Post)
    key_prob[i] = key_prob[i] * comp
  }
  strong_key = which.max(key_prob)
  return(strong_key)
}