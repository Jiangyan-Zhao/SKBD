## generative the time-to-event data
# n: Sample size (number of patients)
# pi: toxicity probability
# dist: the underlying distribution of the time to toxicity outcomes.
# alpha: a number from (0,1) that controls alpha*100% events in (0, 1/2T). 
#              The default is \code{alpha=0.5}. 
# tau: the DLT assessment window

#' @noRd

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
