# ============================================================
# ins_true_prob()
#   Simulation-only: assign the "true toxicity probability" to the inserted dose.
#   Real trial: you do NOT have tox_prob_work; skip this.
#
# Inputs:
#   insert_code   : Scheme A code
#                  - LOC_BELOW_MIN (-1L): insert below minimum
#                  - LOC_ABOVE_MAX (-2L): insert above maximum
#                  - i >= 1             : insert between (i, i+1)
#   d             : current dose index when insertion is triggered (used for boundary cases)
#   p             : current true tox prob vector aligned with dose_set_work (simulation only)
#   target        : target toxicity rate (phi)
#   LOC_BELOW_MIN : default -1L
#   LOC_ABOVE_MAX : default -2L
# ============================================================

ins_true_prob = function(
    insert_code, d, p, target,
    LOC_BELOW_MIN = -1L, LOC_ABOVE_MAX = -2L
) {
  
  J = length(p)
  if (J < 1) stop("p must have positive length.")
  if (!(0 < target && target < 1)) stop("target must be in (0, 1).")
  if (d < 1 || d > J) stop("d out of range.")
  
  p_new = NA_real_
  
  # ----------------------------------------------------------
  # Case 1) Insert BETWEEN (i, i+1)
  #   if p[i] < target < p[i+1] then p_new = target
  #   else p_new = (p[i] + p[i+1]) / 2
  # ----------------------------------------------------------
  if (!is.na(insert_code) && insert_code >= 1) {
    
    i = as.integer(insert_code)
    if (i < 1 || i >= J) stop("insert_code (between) out of range.")
    
    p_new = ifelse (p[i] < target && target < p[i + 1], target, (p[i] + p[i + 1]) / 2)
    
    # ----------------------------------------------------------
    # Case 2) Insert BELOW minimum:
    #   if target < p[d] then p_new = target
    #   else p_new = p[d] / 2
    #
    # In your Scheme A logic, this case should occur when d == 1 (lowest dose),
    # but we keep it general and use p[d] exactly as ADM does.
    # ----------------------------------------------------------
  } else if (!is.na(insert_code) && insert_code == LOC_BELOW_MIN) {
    
    p_new = ifelse (target < p[d], target, p_new = p[d] / 2)
    
    # ----------------------------------------------------------
    # Case 3) Insert ABOVE maximum: 
    #   if target > p[d] then p_new = target
    #   else p_new = p[d] + 0.1
    #
    # This is how ADM generates a higher true tox probability for extrapolation upward.
    # ----------------------------------------------------------
  } else if (!is.na(insert_code) && insert_code == LOC_ABOVE_MAX) {
    
    p_new = ifelse (target > p[d], target, p_new = p[d] + 0.1)
    
  } else {
    stop("unsupported insert_code.")
  }
  
  # Clamp into (0,1) for safety in simulation
  eps_prob = 1e-6
  p_new = max(min(p_new, 1 - eps_prob), eps_prob)
  
  return(p_new)
}
