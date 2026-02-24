## get cutoffs for keys
#' @noRd
get_Key <- function(
    target_DLT, 
    margin_left = 0.05, 
    margin_right = 0.05
){
  target_left = target_DLT - margin_left          # the left margin of target key
  target_right = target_DLT + margin_right        # the right margin of target key
  margin_length = margin_left + margin_right      # the length of target key
  keyLeft = target_left - (floor(target_left / margin_length):1) * margin_length
  KeyRight = target_right + (1:floor((1-target_right) / margin_length)) * margin_length
  keys = unique(c(0, keyLeft, target_left, target_right, KeyRight, 1))
  return(keys)
}

