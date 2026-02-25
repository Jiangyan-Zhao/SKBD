#' @title get cutoffs for keys
#' 
#' @description
#' Construct the key cutoffs used by keyboard-type dose-finding rules. The target key is
#' \eqn{(\phi-\epsilon_1,\phi+\epsilon_2)}, and additional keys with the same width
#' \eqn{\epsilon_1+\epsilon_2} are extended to the left and right until reaching 0 and 1.
#'
#' @details
#' This function returns an increasing vector of unique cutoffs in \eqn{[0,1]} that partition
#' the unit interval into equal-width keys. The cutoffs always include \code{0} and \code{1},
#' as well as the target-key boundaries \code{target_DLT - margin_left} and
#' \code{target_DLT + margin_right}.
#'
#' @param target_DLT Scalar in \eqn{(0,1)}. Target toxicity rate \eqn{\phi}.
#' @param margin_left Nonnegative scalar. Left half-width \eqn{\epsilon_1} of the target key.
#' @param margin_right Nonnegative scalar. Right half-width \eqn{\epsilon_2} of the target key.
#'
#' @return
#' A numeric vector of sorted cutoffs in \eqn{[0,1]} (including \code{0} and \code{1}).
#'
#' @examples
#' get_Key(0.30, 0.05, 0.05)
#' get_Key(0.20, 0.03, 0.07)
#'
#' @export
get_Key <- function(
    target_DLT, 
    margin_left = 0.05, 
    margin_right = 0.05
){
  target_left = target_DLT - margin_left          # the left margin of target key
  target_right = target_DLT + margin_right        # the right margin of target key
  margin_length = margin_left + margin_right      # the length of target key
  key_left = target_left - (floor(target_left / margin_length):1) * margin_length
  key_right = target_right + (1:floor((1-target_right) / margin_length)) * margin_length
  keys = unique(c(0, key_left, target_left, target_right, key_right, 1))
  return(keys)
}

