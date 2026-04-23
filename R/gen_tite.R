#' @title Generate Time-to-DLT Data for TITE Simulations
#'
#' @description
#' Generate patient-level toxicity outcomes and time-to-DLT data under a
#' prespecified time-to-event model for use in TITE dose-finding simulations.
#'
#' @details
#' This function simulates \eqn{n} patients with marginal DLT probability
#' \eqn{\pi} within an assessment window of length \code{tau}. Three
#' time-to-event models are supported:
#'
#' \itemize{
#'   \item \code{"weibull"}: DLT times are generated from a Weibull distribution.
#'   \item \code{"loglogistic"}: DLT times are generated from a log-logistic distribution.
#'   \item \code{"uniform"}: DLT occurrence is first generated from a Bernoulli
#'   distribution with probability \code{pi}, and conditional on DLT, the event
#'   time is sampled uniformly on \eqn{(0, \tau)}.
#' }
#'
#' For \code{"weibull"} and \code{"loglogistic"}, the parameter \code{alpha}
#' controls the fraction of DLT events expected to occur in the first half of
#' the assessment window. Specifically, \code{alpha * pi} is treated as the DLT
#' probability by time \eqn{\tau/2}. For patients without a DLT within
#' \code{tau}, the corresponding event time is recorded as \code{Inf}.
#'
#' @param n Positive integer. Number of patients to simulate.
#' @param pi Scalar in \eqn{(0,1)}. Marginal probability of experiencing a DLT
#' within the assessment window \code{tau}.
#' @param dist Character string specifying the time-to-DLT distribution.
#' Must be one of \code{"weibull"}, \code{"loglogistic"}, or \code{"uniform"}.
#' @param alpha Scalar in \eqn{(0,1)} controlling the early-event proportion for
#' \code{"weibull"} and \code{"loglogistic"} models. Under these models,
#' \code{alpha * pi} is the DLT probability by time \eqn{\tau/2}. This argument
#' is ignored when \code{dist = "uniform"}.
#' @param tau Positive scalar. Length of the DLT assessment window.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{DLT}}{A binary vector of length \code{n}, where 1 indicates a
#'   DLT within \code{tau} and 0 otherwise.}
#'   \item{\code{time_DLT}}{A numeric vector of length \code{n} giving the
#'   event time for each patient. Patients without a DLT within \code{tau} are
#'   assigned \code{Inf}.}
#'   \item{\code{n_DLT}}{Total number of DLTs, equal to \code{sum(DLT)}.}
#' }
#'
#' @examples
#' set.seed(1)
#'
#' out1 <- gen_tite(
#'   n = 6,
#'   pi = 0.25,
#'   dist = "weibull",
#'   alpha = 0.5,
#'   tau = 3
#' )
#' out1$n_DLT
#' out1$DLT
#' out1$time_DLT
#'
#' out2 <- gen_tite(
#'   n = 6,
#'   pi = 0.25,
#'   dist = "uniform",
#'   tau = 3
#' )
#' out2$n_DLT
#'
#' @importFrom stats rbinom runif
#' @export

gen_tite = function(n, pi, 
                    dist = c("weibull", "loglogistic", "uniform"), 
                    alpha = 0.5, tau = 3) {
  ## ---- argument checks ----
  stopifnot(
    n > 0,
    pi > 0, pi < 1,
    alpha > 0, alpha < 1,
    tau > 0
  )
  dist = match.arg(dist)
  
  #---------------------- subroutines ----------------------
  rweib_tite = function(n, pi, pi_half, tau) {
    ## solve parameters for Weibull given pi=1-S(T) and pi_half=1-S(T/2)
    # S(t) = \exp(-\lambda t^\alpha)
    alpha = log(log(1 - pi) / log(1 - pi_half)) / log(2)
    lambda = -log(1 - pi) / (tau^alpha)
    # the inverse transformation sampling method of the Weibull distribution
    t = (-log(runif(n))/ lambda)^(1 / alpha)
    return(t)
  }
  
  rllogis_tite = function(n, pi, pi_half, tau) {
    ## solve parameters for log-logistic given pi=1-S(T) and pi_half=1-S(T/2)
    # S(t) = \frac{1}{1 + \lambda t^\alpha}
    alpha = log((1 / (1 - pi) - 1)/(1 / (1 - pi_half) - 1)) / log(2)
    lambda = (1 / (1 - pi) - 1)/(tau^alpha)
    t = ((1 / runif(n) -1) / lambda)^(1 / alpha)
    return(t)
  }
  #---------------------- end of subroutines ----------------------
  
  #----------------------- initialization ------------------------
  DLT = rep(0, n)   # x_i: the binary toxicity outcome {0, 1}
  time_DLT = rep(0, n) # t_i: the actual time to DLT
  
  #--------------------- data generation ------------------------
  if(dist == "weibull") {
    pi_half = alpha * pi   # alpha*100% event in (0, 1/2T)
    time_DLT = rweib_tite(n, pi, pi_half, tau)
    DLT = as.integer(time_DLT <= tau)
    time_DLT[DLT == 0] = Inf
  }else if (dist == "loglogistic") {
    pi_half = alpha * pi  # alpha*100% event in (0, 1/2T)
    time_DLT = rllogis_tite(n, pi, pi_half, tau)
    DLT = as.integer(time_DLT <= tau)
    time_DLT[DLT == 0] = Inf
  } else if(dist == "uniform") {  # 50% event in (0, 1/2T)
    DLT = rbinom(n, size = 1, prob = pi)
    time_DLT[DLT == 0] = Inf
    time_DLT[DLT == 1] = runif(sum(DLT), 0, tau)
  }
  
  return(list(DLT = DLT, time_DLT = time_DLT, n_DLT = sum(DLT)))
}
