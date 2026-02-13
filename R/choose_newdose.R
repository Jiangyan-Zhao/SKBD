# ============================================================
# choose_newdose()
#   Choose newdose in an eligible interval (dl, dr) by maximizing q(d) on a grid.
#
# Inputs:
#   dl, dr           : interval endpoints (dl < dr)
#   dose_set_work    : current working dose grid
#   n_dlt_work       : DLT counts on working grid
#   n_treated_work   : treated counts on working grid
#   M                : grid size for searching (e.g., 100)
#
# Output:
#   list(newdose = ..., q_max = ..., grid = ..., q_grid = ...)
# ============================================================

choose_newdose = function(
    dl, dr, dose_set_work, n_dlt_work, n_treated_work,
    pri_alpha, pri_beta, key_L, key_U,
    symmetric, theta, M = 100
){
  if (!is.finite(dl) || !is.finite(dr) || dr <= dl) {
    return(list(newdose = NA_real_, q_max = NA_real_))
  }
  if (!is.numeric(M) || length(M) != 1 || M < 5) M <- 50L
  
  obs_idx = which(n_treated_work > 0)
  if (length(obs_idx) == 0) {
    return(list(newdose = NA_real_, q_max = NA_real_))
  }
  
  dose_obs = dose_set_work[obs_idx]
  y_obs    = n_dlt_work[obs_idx]
  n_obs    = n_treated_work[obs_idx]
  
  eps = 1e-10
  if ((dr - dl) <= 2*eps) return(list(newdose = NA_real_, q_max = NA_real_))
  grid = seq(dl + eps, dr - eps, length.out = M)
  
  q_grid = rep(NA_real_, length(grid))
  for (g in seq_along(grid)) {
    k = kernel(grid[g], dose_obs, symmetric, theta)
    km = max(k)
    if (!is.finite(km) || km <= 0) next
    k = k / km
    den = sum(k)
    if (!is.finite(den) || den <= 0) next
    w = k / den
    
    alpha_d = max(pri_alpha + sum(w * y_obs), 1e-12)
    beta_d  = max(pri_beta  + sum(w * (n_obs - y_obs)), 1e-12)
    
    q_grid[g] = pbeta(key_U, alpha_d, beta_d) - pbeta(key_L, alpha_d, beta_d)
  }
  
  if (all(is.na(q_grid))) return(list(newdose = NA_real_, q_max = NA_real_, grid = grid, q_grid = q_grid))
  idx = which.max(q_grid)
  cand = grid[idx]
  
  tol = 1e-12
  if (any(abs(dose_set_work - cand) < tol)) cand <- NA_real_
  
  list(newdose = cand, q_max = q_grid[idx], grid = grid, q_grid = q_grid)
}
