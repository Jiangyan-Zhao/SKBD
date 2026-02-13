# ============================================================
# insert_check()  [BKP + ADM-style PAVA monotonicity]
#   BKP-based insertion trigger (Eq.(5) logic) + monotonic stabilization.
#
# Step A: compute raw region probabilities at each dose:
#   under[j]  = Pr(pi_j <= key_L | D) = pbeta(key_L; post_alpha[j], post_beta[j])
#   over[j]   = Pr(pi_j >  key_U | D) = 1 - pbeta(key_U; post_alpha[j], post_beta[j])
#   target[j] = Pr(key_L < pi_j <= key_U | D) = pbeta(key_U)-pbeta(key_L)
#
# Step B: enforce monotonicity via PAVA (ADM style):
#   under_pava  : non-increasing with dose index
#   over_pava   : non-decreasing with dose index
#   target_pava : complement, clamped at >=0
#
# Step C: apply Eq.(5) trigger using (optionally) PAVA-adjusted under/over.
#
# Output scheme (Scheme A):
#   insert_code = NA         : no insertion
#   insert_code = -1L        : below minimum
#   insert_code = -2L        : above maximum
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
    use_pava = TRUE,          # NEW: whether to apply ADM-style monotonicity
    return_probs = TRUE,       # NEW: return prob vectors for debugging
    LOC_BELOW_MIN = -1L,
    LOC_ABOVE_MAX = -2L
) {
  
  J = length(post_alpha)
  
  if (length(post_beta) != J) stop("post_alpha and post_beta must have same length.")
  if (length(n_treated) != J) stop("n_treated must have same length as post_alpha.")
  if (j < 1 || j > J) stop("j must be in 1..J.")
  if (!(0 < key_L && key_L < key_U && key_U < 1)) stop("require 0 < key_L < key_U < 1.")
  if (!(0 <= C1 && C1 <= 1 && 0 <= C2 && C2 <= 1)) stop("C1,C2 must be in [0,1].")
  
  if (use_pava) {
    if (!exists("pava", mode = "function")) {
      stop("insert_check(): use_pava=TRUE but pava() not found. Please define pava() first.")
    }
  }
  
  # -----------------------------
  # Step A: raw region probabilities at each discrete dose
  # -----------------------------
  prob_under = pbeta(key_L, post_alpha, post_beta)
  prob_over  = 1 - pbeta(key_U, post_alpha, post_beta)
  prob_target = (pbeta(key_U, post_alpha, post_beta) - pbeta(key_L, post_alpha, post_beta))
  
  # -----------------------------
  # Step B: ADM-style PAVA monotonic adjustment (optional)
  # -----------------------------
  if (use_pava) {
    prob_under_adj = rev(pava(rev(prob_under)))   # enforce non-increasing
    prob_over_adj  = pava(prob_over)              # enforce non-decreasing
    prob_target_adj = pmax(0, 1 - prob_under_adj - prob_over_adj)
  } else {
    prob_under_adj = prob_under
    prob_over_adj  = prob_over
    prob_target_adj = prob_target
  }
  
  # -----------------------------
  # Step C: Eq.(5) trigger using adjusted probabilities
  # -----------------------------
  insert_code = NA_integer_
  
  # Left interval: (d_{j-1}, d_j)  <=> l = j
  if (j == 1) {
    # Boundary (d0, d1): pi_0 = 0 => Pr(pi_0 < key_L) = 1, so only need right bracket evidence.
    left_ok = (n_treated[1] > 0) && (prob_over_adj[1] > C2)
    if (left_ok) insert_code = LOC_BELOW_MIN
  } else {
    left_ok = (n_treated[j - 1] > 0) && (n_treated[j] > 0) &&
      (prob_under_adj[j - 1] > C1) && (prob_over_adj[j] > C2)
    if (left_ok) insert_code = as.integer(j - 1)
  }
  
  # Right interval: (d_j, d_{j+1}) <=> l = j+1
  if (is.na(insert_code)) {
    if (j == J) {
      # Boundary (dJ, dJ+1): pi_{J+1} = 1 => Pr(pi_{J+1} > key_U) = 1, so only need left bracket evidence.
      right_ok = (n_treated[J] > 0) && (prob_under_adj[J] > C1)
      if (right_ok) insert_code = LOC_ABOVE_MAX
    } else {
      right_ok = (n_treated[j] > 0) && (n_treated[j + 1] > 0) &&
        (prob_under_adj[j] > C1) && (prob_over_adj[j + 1] > C2)
      if (right_ok) insert_code = as.integer(j)
    }
  }
  
  need_insert = !is.na(insert_code)
  
  out = list(
    need_insert = need_insert,
    insert_code = insert_code,
    LOC_BELOW_MIN = LOC_BELOW_MIN,
    LOC_ABOVE_MAX = LOC_ABOVE_MAX
  )
  
  if (return_probs) {
    out$prob_under_raw = prob_under
    out$prob_over_raw  = prob_over
    out$prob_target_raw = prob_target
    
    out$prob_under = prob_under_adj
    out$prob_over  = prob_over_adj
    out$prob_target = prob_target_adj
  }
  
  return(out)
}
