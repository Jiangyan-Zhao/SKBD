## Kernel Function

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
