
## Kernel Function
#' @noRd
kernel <- function(dose, dose_set, symmetric = FALSE, theta = NULL){
  # dose--std dose(not index)
  # dose_set--std test dose

  if(symmetric){
    k = exp(-theta * (dose - dose_set)^2)
  }else{
    n_dose = length(dose_set)
    n_dose_left = sum(dose_set < dose)
    
    theta_set = c(rep(theta[1], n_dose_left), rep(theta[2], n_dose - n_dose_left))
    
    k = exp(-theta_set * (dose - dose_set)^2)
  }
  
  return(k)
}

## get cutoffs for keys
#' @noRd
get_Key <- function(target_DLT, margin_left = 0.05, margin_right = 0.05){
  target_left = target_DLT - margin_left          # the left margin of target key
  target_right = target_DLT + margin_right        # the right margin of target key
  margin_length = margin_left + margin_right      # the length of target key
  keyLeft = target_left - (floor(target_left / margin_length):1) * margin_length
  KeyRight = target_right + (1:floor((1-target_right) / margin_length)) * margin_length
  keys = unique(c(0, keyLeft, target_left, target_right, KeyRight, 1))
  return(keys)
}

## get stongest key
#' @noRd
get_strongKey <- function(alpha_Post, beta_Post, keys, 
                          margin_left = 0.05, margin_right = 0.05) {
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

## The Pooled Adjacent Violators Algorithm
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

## the pseudo uniform algorithm
#' @noRd
PUA = function(dose_set, target_prob, n_scenarios = 1000, 
               margin_left = 0.05, margin_right = 0.05, seed = 6){
  set.seed(seed)
  n_dose = length(dose_set)
  tox_prob = matrix(NA, nrow = n_scenarios, ncol = n_dose)
  
  # browser()
  s = 1
  while(s <= n_scenarios){
    ## Step 1: select j uniformly
    d = sample(1:n_dose, 1)
    
    ## Step 2: sample M ~ Beta(max{J-j, 0.5}, 1)
    M = rbeta(1, max(n_dose - d, 0.5), 1)
    
    ## Step 3: rejective sampling
    B = target_prob + (1 - target_prob) * M # upper bound B
    
    n_repeat = 0
    repeat {
      n_repeat = n_repeat + 1
      if(n_repeat > 100) break
      
      tox_prob[s, ] = sort(runif(n_dose, 0, B))
      if(d != which.min(abs(tox_prob[s, ] - target_prob))) next
      
      gap_left  = if(d > 1) tox_prob[s, d] - tox_prob[s, d - 1] else NA
      gap_right = if(d < n_dose) tox_prob[s, d + 1] - tox_prob[s, d] else NA
      
      # well definited MTD, see 21_Zhou_PS_iBOIN
      definited_condition = (
        tox_prob[s, d] >= target_prob - margin_left &&
          tox_prob[s, d] <= target_prob + margin_right &&
          (is.na(gap_left)  || (gap_left  > 0.05 && gap_left  < 0.3)) &&
          (is.na(gap_right) || (gap_right > 0.05 && gap_right < 0.3))
      )
      
      if(definited_condition) {
        s = s + 1
        break
      }
    }
    
    # cat("tox_prob[", s - 1, ", ] = (", sep = "")
    # cat(sprintf("%.3f", tox_prob[s - 1, ]), sep = ", ")
    # cat(")\n", sep = "")
  }
  return(tox_prob)
}


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
