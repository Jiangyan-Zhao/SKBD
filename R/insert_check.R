# ============================================================
# insert_check()
#   BKP-based insertion trigger exactly matching Eq.(5) in Section 2.4.
#
# Core trigger (Eq. 5):
#   For l in {j, j+1}:
#     Pr(pi_{l-1} < key_L | D) > C1    AND    Pr(pi_l > key_U | D) > C2
#
# "Evaluate only when both bracketing doses have accrued information":
#   Check Eq.(5) only if n_{l-1} > 0 and n_l > 0 (paper rule),
#   with boundary handling using auxiliary pi_0=0 and pi_{J+1}=1.
#   In code:
#     - For below-min interval (d0, d1): left event is always TRUE if key_L>0; require n1>0.
#     - For above-max interval (dJ, dJ+1): right event is always TRUE if key_U<1; require nJ>0.
#
# Scheme A output:
#   insert_code = NA         : no insertion
#   insert_code = -1L        : below minimum (d0, d1)
#   insert_code = -2L        : above maximum (dJ, dJ+1)
#   insert_code = i >= 1     : between (i, i+1)
# ============================================================

insert_check = function(
    j,
    post_alpha,
    post_beta,
    key_L,
    key_U,
    C1,
    C2,
    n_treated,
    LOC_BELOW_MIN = -1L,
    LOC_ABOVE_MAX = -2L
) {
  
  J = length(post_alpha)
  
  if (length(post_beta) != J) stop("post_alpha and post_beta must have same length.")
  if (length(n_treated) != J) stop("n_treated must have same length as post_alpha.")
  if (j < 1 || j > J) stop("j must be in 1..J.")
  if (!(0 < key_L && key_L < key_U && key_U < 1)) stop("require 0 < key_L < key_U < 1.")
  if (!(0 <= C1 && C1 <= 1 && 0 <= C2 && C2 <= 1)) stop("C1,C2 must be in [0,1].")
  
  insert_code = NA_integer_
  
  # Helper probabilities at a *discrete* dose level k using BKP posterior Beta(post_alpha[k], post_beta[k])
  prob_under = function(k) {
    pbeta(key_L, post_alpha[k], post_beta[k])  # Pr(pi_k < key_L)
  }
  prob_over = function(k) {
    1 - pbeta(key_U, post_alpha[k], post_beta[k])  # Pr(pi_k > key_U)
  }
  
  # ----------------------------------------------------------
  # Check left adjacent interval: (d_{j-1}, d_j) corresponds to l = j
  # Condition: Pr(pi_{j-1} < key_L) > C1  AND  Pr(pi_j > key_U) > C2
  # ----------------------------------------------------------
  left_ok = FALSE
  
  if (j == 1) {
    # Boundary interval (d0, d1): treat pi_0 = 0 => Pr(pi_0 < key_L) = 1 (since key_L>0)
    # Only require the right bracket has accrued info (n1>0) and is "strongly over"
    left_ok = (n_treated[1] > 0) && (prob_over(1) > C2)
    if (left_ok) insert_code = LOC_BELOW_MIN
    
  } else {
    # Interior interval: require both bracketing doses have data
    left_ok = (n_treated[j - 1] > 0) && (n_treated[j] > 0) &&
      (prob_under(j - 1) > C1) && (prob_over(j) > C2)
    if (left_ok) insert_code = as.integer(j - 1)  # between (j-1, j)
  }
  
  # ----------------------------------------------------------
  # Check right adjacent interval: (d_j, d_{j+1}) corresponds to l = j+1
  # Condition: Pr(pi_j < key_L) > C1  AND  Pr(pi_{j+1} > key_U) > C2
  # ----------------------------------------------------------
  right_ok = FALSE
  
  if (is.na(insert_code)) {
    if (j == J) {
      # Boundary interval (dJ, dJ+1): treat pi_{J+1} = 1 => Pr(pi_{J+1} > key_U) = 1 (since key_U<1)
      # Only require the left bracket has accrued info (nJ>0) and is "strongly under"
      right_ok = (n_treated[J] > 0) && (prob_under(J) > C1)
      if (right_ok) insert_code = LOC_ABOVE_MAX
      
    } else {
      # Interior interval: require both bracketing doses have data
      right_ok = (n_treated[j] > 0) && (n_treated[j + 1] > 0) &&
        (prob_under(j) > C1) && (prob_over(j + 1) > C2)
      if (right_ok) insert_code = as.integer(j)  # between (j, j+1)
    }
  }
  
  need_insert = !is.na(insert_code)
  
  return(list(
    need_insert = need_insert,
    insert_code = insert_code,
    LOC_BELOW_MIN = LOC_BELOW_MIN,
    LOC_ABOVE_MAX = LOC_ABOVE_MAX
  ))
}
