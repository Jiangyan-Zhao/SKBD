#' @title Pool-Adjacent-Violators Algorithm
#'
#' @description
#' Perform isotonic regression under a monotone nondecreasing constraint using
#' the pool-adjacent-violators algorithm (PAVA).
#'
#' @details
#' This function takes a numeric vector \code{x} and returns its isotonic
#' projection under the constraint
#' \deqn{x_1 \le x_2 \le \cdots \le x_n.}
#'
#' When adjacent elements violate monotonicity, they are pooled and replaced by
#' their weighted average. This process is repeated until the entire sequence is
#' monotone nondecreasing.
#'
#' The function is used in the package to stabilize estimated dose-toxicity
#' curves before final dose selection or before applying insertion rules.
#'
#' @param x Numeric vector to be monotonized.
#' @param wt Numeric vector of nonnegative weights with the same length as
#' \code{x}. Defaults to equal weights. Larger weights give more influence to
#' the corresponding elements when adjacent violating blocks are pooled.
#'
#' @return
#' A numeric vector of the same length as \code{x}, giving the isotonic
#' regression fit under a nondecreasing constraint.
#'
#' @examples
#' x <- c(0.10, 0.25, 0.20, 0.40)
#' pava(x)
#'
#' wt <- c(1, 1, 3, 1)
#' pava(x, wt = wt)
#'
#' @noRd

pava <- function(x, wt = rep(1, length(x))) {
  n = length(x)
  if (n <= 1) {
    return(x)
  }
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets = 1:n
  repeat {
    viol = (as.vector(diff(x)) < 0)
    if (!(any(viol))) {
      break
    }
    i = min((1:(n - 1))[viol])
    lvl1 = lvlsets[i]
    lvl2 = lvlsets[i + 1]
    ilvl = (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] = sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] = lvl1
  }
  return(x)
}