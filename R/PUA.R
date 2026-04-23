## the pseudo uniform algorithm
#' @export
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