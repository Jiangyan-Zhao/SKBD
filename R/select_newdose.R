# ============================================================
# propose_newdose_bkp()
#   Choose inserted dose d^dagger within eligible interval I by maximizing q(d).
#
# Inputs:
#   insert_code    : Scheme A code from insert_check_bkp()
#   dose_set_work  : working dose grid on standardized scale in [d0, dJ1] (typically [0,1])
#   n_dlt_work     : observed DLT counts at each working dose
#   n_treated_work : number treated at each working dose
#   key_L, key_U   : target key bounds
#   kernel_fun     : function(d, ds, ...) returning nonnegative kernel values k(d, ds)
#                   - d is scalar, ds is vector
#   M              : grid size (paper suggests e.g. 200)
#   d0, dJ1        : allowable standardized endpoints (paper: d0=0, dJ+1=1)
#
# Output:
#   newdose        : selected inserted dose (standardized), or NA_real_ if cannot propose
#   q_max          : max q(d)
# ============================================================

propose_newdose_bkp = function(
    insert_code,
    dose_set_work,
    n_dlt_work,
    n_treated_work,
    key_L,
    key_U,
    kernel_fun,
    M = 200,
    d0 = 0,
    dJ1 = 1,
    LOC_BELOW_MIN = -1L,
    LOC_ABOVE_MAX = -2L,
    ...
) {
  
  J = length(dose_set_work)
  
  if (length(n_dlt_work) != J || length(n_treated_work) != J) {
    stop("propose_newdose_bkp(): n_dlt_work and n_treated_work must align with dose_set_work.")
  }
  
  # Observed set S = {s: n_s > 0}, as required by the paper
  S = which(n_treated_work > 0)
  if (length(S) == 0) {
    return(list(newdose = NA_real_, q_max = NA_real_))
  }
  
  ds_obs = dose_set_work[S]
  y_obs  = n_dlt_work[S]
  n_obs  = n_treated_work[S]
  
  # ----------------------------------------------------------
  # Determine eligible interval I = (dl, dr) from insert_code
  #   below_min : (d0, d1)
  #   above_max : (dJ, dJ1)
  #   between i : (d_i, d_{i+1})
  # ----------------------------------------------------------
  dl = NA_real_
  dr = NA_real_
  
  if (is.na(insert_code)) {
    return(list(newdose = NA_real_, q_max = NA_real_))
  }
  
  if (insert_code == LOC_BELOW_MIN) {
    dl = d0
    dr = dose_set_work[1]
    
  } else if (insert_code == LOC_ABOVE_MAX) {
    dl = dose_set_work[J]
    dr = dJ1
    
  } else if (insert_code >= 1) {
    i = as.integer(insert_code)
    if (i < 1 || i >= J) stop("propose_newdose_bkp(): invalid insert_code for between().")
    dl = dose_set_work[i]
    dr = dose_set_work[i + 1]
  }
  
  if (!is.finite(dl) || !is.finite(dr) || dr <= dl) {
    return(list(newdose = NA_real_, q_max = NA_real_))
  }
  
  # Avoid endpoints so we don't “re-insert” an existing dose
  eps = 1e-8
  grid = seq(dl + eps, dr - eps, length.out = M)
  
  # ----------------------------------------------------------
  # Compute q(d) on grid:
  #   1) weights w(d; ds) = k(d, ds) / sum_u k(d, du)
  #   2) alpha(d) = 1 + sum_s w * y_s
  #      beta(d)  = 1 + sum_s w * (n_s - y_s)
  #   3) q(d) = pbeta(key_U, alpha, beta) - pbeta(key_L, alpha, beta)
  # ----------------------------------------------------------
  q = rep(NA_real_, length(grid))
  
  for (g in 1:length(grid)) {
    d = grid[g]
    
    k = kernel_fun(d, ds_obs, ...)
    k = pmax(0, as.numeric(k))
    
    if (all(k == 0)) {
      # If kernel degenerates to all zeros, fall back to equal weights on S
      w = rep(1 / length(S), length(S))
    } else {
      w = k / sum(k)
    }
    
    alpha_d = 1 + sum(w * y_obs)
    beta_d  = 1 + sum(w * (n_obs - y_obs))
    
    q[g] = pbeta(key_U, alpha_d, beta_d) - pbeta(key_L, alpha_d, beta_d)
  }
  
  idx = which.max(q)
  return(list(newdose = grid[idx], q_max = q[idx]))
}
