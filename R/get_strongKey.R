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