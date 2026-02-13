get_OC_Insert_SKBD <- function(
    target_prob, tox_prob,
    n_cohort, cohort_size,
    dose_set, dose_range,
    C1 = 0.6, C2 = 0.6,
    start_dose = 1,
    margin_left = 0.05, margin_right = 0.05,
    n_earlystop = 1000,
    cutoff_elimin = 0.95,
    extra_safe = FALSE, offset = 0.05,
    symmetric = FALSE, theta = NULL,
    shared = TRUE,          # kept for compatibility; decision rule below is probability-dominance
    light_return = TRUE,
    n_trial = 1000, seed = 6
) {
  set.seed(seed)
  
  # ============================================================
  # Helper functions (self-contained)
  # ============================================================
  
  int1 <- function(x) is.numeric(x) && length(x) == 1 && !is.na(x) && x %% 1 == 0 && x >= 1
  
  # PAVA (nondecreasing)
  pava <- function(y) stats::isoreg(y)$yf
  
  # Asymmetric / symmetric kernel on standardized dose in [0,1]
  kernel <- function(dj, ds, symmetric, theta) {
    dd <- (dj - ds)
    if (symmetric) {
      exp(-theta * dd^2)
    } else {
      # ds < dj uses theta1 (more decay), ds >= dj uses theta2 (less decay)
      th1 <- theta[1]; th2 <- theta[2]
      out <- ifelse(ds < dj, exp(-th1 * dd^2), exp(-th2 * dd^2))
      out
    }
  }
  
  # choose_newdose(): grid search argmax q(d)
  choose_newdose <- function(
    dl, dr, dose_set_work, n_dlt_work, n_treated_work,
    pri_alpha, pri_beta, key_L, key_U, symmetric, theta, M = 100
  ) {
    if (!is.finite(dl) || !is.finite(dr) || dr <= dl) {
      return(list(newdose = NA_real_, q_max = NA_real_, grid = NA, q_grid = NA))
    }
    obs_idx <- which(n_treated_work > 0)
    if (length(obs_idx) == 0) {
      return(list(newdose = NA_real_, q_max = NA_real_, grid = NA, q_grid = NA))
    }
    
    dose_obs <- dose_set_work[obs_idx]
    y_obs    <- n_dlt_work[obs_idx]
    n_obs    <- n_treated_work[obs_idx]
    
    eps <- 1e-10
    if ((dr - dl) <= 2 * eps) {
      return(list(newdose = NA_real_, q_max = NA_real_, grid = NA, q_grid = NA))
    }
    
    grid <- seq(dl + eps, dr - eps, length.out = M)
    q_grid <- rep(NA_real_, length(grid))
    
    for (g in seq_along(grid)) {
      k <- kernel(grid[g], dose_obs, symmetric, theta)
      
      # stabilize normalization (avoid underflow -> all zeros)
      km <- max(k)
      if (!is.finite(km) || km <= 0) next
      k <- k / km
      
      den <- sum(k)
      if (!is.finite(den) || den <= 1e-12) next
      w <- k / den
      
      alpha_d <- max(pri_alpha + sum(w * y_obs), 1e-12)
      beta_d  <- max(pri_beta  + sum(w * (n_obs - y_obs)), 1e-12)
      
      q_grid[g] <- stats::pbeta(key_U, alpha_d, beta_d) - stats::pbeta(key_L, alpha_d, beta_d)
    }
    
    if (all(is.na(q_grid))) {
      return(list(newdose = NA_real_, q_max = NA_real_, grid = grid, q_grid = q_grid))
    }
    
    idx <- which.max(q_grid)
    newdose <- grid[idx]
    
    # tolerance de-dup
    tol <- 1e-12
    if (any(abs(dose_set_work - newdose) < tol)) newdose <- NA_real_
    
    list(newdose = newdose, q_max = q_grid[idx], grid = grid, q_grid = q_grid)
  }
  
  # insert_check(): BKP + (optional) ADM-style PAVA stabilization on under/over probs
  insert_check <- function(
    j, post_alpha, post_beta, key_L, key_U, C1, C2, n_treated,
    use_pava = TRUE, return_probs = FALSE,
    LOC_BELOW_MIN = -1L, LOC_ABOVE_MAX = -2L
  ) {
    J <- length(post_alpha)
    if (length(post_beta) != J) stop("post_alpha and post_beta must have same length.")
    if (length(n_treated) != J) stop("n_treated must have same length as post_alpha.")
    if (j < 1 || j > J) stop("j must be in 1..J.")
    if (!(0 < key_L && key_L < key_U && key_U < 1)) stop("require 0 < key_L < key_U < 1.")
    if (!(0 <= C1 && C1 <= 1 && 0 <= C2 && C2 <= 1)) stop("C1,C2 must be in [0,1].")
    
    prob_under  <- stats::pbeta(key_L, post_alpha, post_beta)
    prob_over   <- 1 - stats::pbeta(key_U, post_alpha, post_beta)
    prob_target <- stats::pbeta(key_U, post_alpha, post_beta) - stats::pbeta(key_L, post_alpha, post_beta)
    
    if (use_pava) {
      prob_under_adj  <- rev(pava(rev(prob_under)))  # non-increasing
      prob_over_adj   <- pava(prob_over)             # non-decreasing
      prob_target_adj <- pmax(0, 1 - prob_under_adj - prob_over_adj)
      # (optional) renormalize to sum to 1
      s <- prob_under_adj + prob_over_adj + prob_target_adj
      s <- pmax(s, 1e-12)
      prob_under_adj  <- prob_under_adj / s
      prob_over_adj   <- prob_over_adj / s
      prob_target_adj <- prob_target_adj / s
    } else {
      prob_under_adj  <- prob_under
      prob_over_adj   <- prob_over
      prob_target_adj <- prob_target
    }
    
    insert_code <- NA_integer_
    
    # Left interval (d_{j-1}, d_j)
    if (j == 1) {
      # boundary: pi0=0 => only need right bracket evidence
      left_ok <- (n_treated[1] > 0) && (prob_over_adj[1] > C2)
      if (left_ok) insert_code <- LOC_BELOW_MIN
    } else {
      left_ok <- (n_treated[j - 1] > 0) && (n_treated[j] > 0) &&
        (prob_under_adj[j - 1] > C1) && (prob_over_adj[j] > C2)
      if (left_ok) insert_code <- as.integer(j - 1)
    }
    
    # Right interval (d_j, d_{j+1})
    if (is.na(insert_code)) {
      if (j == J) {
        # boundary: pi_{J+1}=1 => only need left bracket evidence
        right_ok <- (n_treated[J] > 0) && (prob_under_adj[J] > C1)
        if (right_ok) insert_code <- LOC_ABOVE_MAX
      } else {
        right_ok <- (n_treated[j] > 0) && (n_treated[j + 1] > 0) &&
          (prob_under_adj[j] > C1) && (prob_over_adj[j + 1] > C2)
        if (right_ok) insert_code <- as.integer(j)
      }
    }
    
    out <- list(
      need_insert = !is.na(insert_code),
      insert_code = insert_code,
      LOC_BELOW_MIN = LOC_BELOW_MIN,
      LOC_ABOVE_MAX = LOC_ABOVE_MAX
    )
    if (return_probs) {
      out$prob_under_raw  <- prob_under
      out$prob_over_raw   <- prob_over
      out$prob_target_raw <- prob_target
      out$prob_under      <- prob_under_adj
      out$prob_over       <- prob_over_adj
      out$prob_target     <- prob_target_adj
    }
    out
  }
  
  # ins_true_prob(): simulation-only true toxicity for inserted dose (ADM-style)
  ins_true_prob <- function(
    insert_code, d, p, target,
    LOC_BELOW_MIN = -1L, LOC_ABOVE_MAX = -2L, eps_prob = 1e-6
  ) {
    J <- length(p)
    if (J < 1) stop("ins_true_prob(): p must have positive length.")
    if (!(0 < target && target < 1)) stop("ins_true_prob(): target must be in (0,1).")
    if (d < 1 || d > J) stop("ins_true_prob(): d out of range.")
    
    p_new <- NA_real_
    
    if (!is.na(insert_code) && insert_code >= 1) {
      i <- as.integer(insert_code)
      if (i < 1 || i >= J) stop("ins_true_prob(): insert_code (between) out of range.")
      if (p[i] < target && target < p[i + 1]) p_new <- target else p_new <- (p[i] + p[i + 1]) / 2
    } else if (!is.na(insert_code) && insert_code == LOC_BELOW_MIN) {
      # ADM: if target < p[d] -> target else p[d]/2
      if (target < p[d]) p_new <- target else p_new <- p[d] / 2
    } else if (!is.na(insert_code) && insert_code == LOC_ABOVE_MAX) {
      # ADM: if target > p[d] -> target else p[d] + 0.1
      if (target > p[d]) p_new <- target else p_new <- p[d] + 0.1
    } else {
      stop("ins_true_prob(): unsupported insert_code.")
    }
    
    p_new <- max(min(p_new, 1 - eps_prob), eps_prob)
    p_new
  }
  
  # BKP posterior update for ALL doses on current working grid
  bkp_update_all <- function(dose_set_work, n_dlt_work, n_treated_work, pri_alpha, pri_beta, symmetric, theta) {
    n_dose_work <- length(dose_set_work)
    obs_idx <- which(n_treated_work > 0)
    
    post_alpha <- rep(pri_alpha, n_dose_work)
    post_beta  <- rep(pri_beta,  n_dose_work)
    
    if (length(obs_idx) == 0) return(list(post_alpha = post_alpha, post_beta = post_beta))
    
    dose_obs <- dose_set_work[obs_idx]
    y_obs    <- n_dlt_work[obs_idx]
    m_obs    <- n_treated_work[obs_idx]
    z_obs    <- m_obs - y_obs
    
    K <- matrix(0, nrow = n_dose_work, ncol = length(obs_idx))
    W <- matrix(0, nrow = n_dose_work, ncol = length(obs_idx))
    for (jj in 1:n_dose_work) {
      K[jj, ] <- kernel(dose_set_work[jj], dose_obs, symmetric, theta)
      km <- max(K[jj, ])
      if (!is.finite(km) || km <= 0) km <- 1
      K[jj, ] <- K[jj, ] / km
      den <- sum(K[jj, ])
      if (!is.finite(den) || den <= 1e-12) den <- 1e-12
      W[jj, ] <- K[jj, ] / den
    }
    post_alpha <- pri_alpha + as.numeric(W %*% y_obs)
    post_beta  <- pri_beta  + as.numeric(W %*% z_obs)
    list(post_alpha = post_alpha, post_beta = post_beta)
  }
  
  # Next-dose decision at current dose index d, based on BKP posterior at d
  decide_next_dose <- function(d, post_alpha, post_beta, key_L, key_U) {
    p_under  <- stats::pbeta(key_L, post_alpha[d], post_beta[d])
    p_over   <- 1 - stats::pbeta(key_U, post_alpha[d], post_beta[d])
    p_target <- max(0, 1 - p_under - p_over)
    
    if (p_under > max(p_over, p_target)) return("E")  # escalate
    if (p_over  > max(p_under, p_target)) return("D") # de-escalate
    "S"                                                # stay
  }
  
  # move to nearest non-eliminated dose in given direction
  move_to_safe <- function(d, is_eliminated, dir = c("up", "down")) {
    dir <- match.arg(dir)
    J <- length(is_eliminated)
    if (dir == "up") {
      cand <- (d + 1):J
      cand <- cand[!is_eliminated[cand]]
      if (length(cand) == 0) return(d)
      return(min(cand))
    } else {
      cand <- 1:(d - 1)
      if (length(cand) == 0) return(d)
      cand <- cand[!is_eliminated[cand]]
      if (length(cand) == 0) return(d)
      return(max(cand))
    }
  }
  
  # ============================================================
  # Argument checks
  # ============================================================
  
  # target
  if (!is.numeric(target_prob) || length(target_prob) != 1 || is.na(target_prob) ||
      target_prob <= 0 || target_prob >= 1) {
    stop("`target_prob` must be a single number in (0, 1).")
  }
  
  # margins
  if (!is.numeric(margin_left) || !is.numeric(margin_right) ||
      length(margin_left) != 1 || length(margin_right) != 1 ||
      anyNA(c(margin_left, margin_right)) || margin_left < 0 || margin_right < 0) {
    stop("`margin_left` and `margin_right` must be nonnegative scalars.")
  }
  key_L <- target_prob - margin_left
  key_U <- target_prob + margin_right
  if (key_L <= 0 || key_U >= 1 || key_L >= key_U) {
    stop("Invalid target key: require 0 < target_prob - margin_left < target_prob + margin_right < 1.")
  }
  
  # tox_prob
  if (!is.numeric(tox_prob) || anyNA(tox_prob) || any(tox_prob < 0 | tox_prob > 1)) {
    stop("`tox_prob` must be numeric in [0, 1] with no NA.")
  }
  
  # n_cohort / cohort_size / n_trial
  if (!int1(n_cohort))    stop("`n_cohort` must be a positive integer.")
  if (!int1(cohort_size)) stop("`cohort_size` must be a positive integer.")
  if (!int1(n_trial))     stop("`n_trial` must be a positive integer.")
  
  # n_earlystop: per-dose cap (keep your semantics)
  if (!int1(n_earlystop)) stop("`n_earlystop` must be a positive integer.")
  
  # cutoff / extra safety
  if (!is.numeric(cutoff_elimin) || length(cutoff_elimin) != 1 || is.na(cutoff_elimin) ||
      cutoff_elimin <= 0 || cutoff_elimin >= 1) {
    stop("`cutoff_elimin` must be a scalar in (0, 1).")
  }
  if (!is.numeric(offset) || length(offset) != 1 || is.na(offset) || offset < 0 || offset >= 0.5) {
    stop("`offset` must be a scalar in [0, 0.5).")
  }
  if (extra_safe && (cutoff_elimin - offset) <= 0) {
    stop("When `extra_safe=TRUE`, require `cutoff_elimin - offset > 0`.")
  }
  
  # dose_set vs tox_prob
  if (!is.numeric(dose_set) || anyNA(dose_set)) stop("`dose_set` must be numeric with no NA.")
  if (length(dose_set) != length(tox_prob)) {
    stop("`dose_set` must have the same length as `tox_prob`.")
  }
  if (anyDuplicated(dose_set)) stop("`dose_set` must have unique values.")
  if (is.unsorted(dose_set, strictly = TRUE)) stop("`dose_set` must be strictly increasing (sorted).")
  
  # dose_range interval
  if (!is.numeric(dose_range) || length(dose_range) != 2 || anyNA(dose_range)) {
    stop("`dose_range` must be numeric length-2: c(lower, upper).")
  }
  if (dose_range[1] >= dose_range[2]) stop("`dose_range[1]` must be < `dose_range[2]`.")
  if (min(dose_set) < dose_range[1] - 1e-12 || max(dose_set) > dose_range[2] + 1e-12) {
    stop("All `dose_set` values must lie inside `dose_range`.")
  }
  
  # start_dose index
  if (!int1(start_dose) || start_dose < 1 || start_dose > length(tox_prob)) {
    stop("`start_dose` must be an integer in 1:n_dose.")
  }
  
  # Standardize doses to [0,1]
  dose_set_std <- (dose_set - dose_range[1]) / (dose_range[2] - dose_range[1])
  
  # theta defaults
  if (is.null(theta)) {
    dose_diff_min <- min(diff(dose_set_std))
    theta <- if (symmetric) {
      -log(0.5) / dose_diff_min^2
    } else {
      c(-log(0.2), -log(0.8)) / dose_diff_min^2
    }
  }
  if (symmetric) {
    if (!is.numeric(theta) || length(theta) != 1 || is.na(theta) || theta <= 0) {
      stop("When `symmetric=TRUE`, `theta` must be a single positive number.")
    }
  } else {
    if (!is.numeric(theta) || length(theta) != 2 || anyNA(theta) || any(theta <= 0)) {
      stop("When `symmetric=FALSE`, `theta` must be length-2 positive vector c(theta1, theta2).")
    }
    if (!(theta[1] > theta[2])) {
      stop("For asymmetric kernel, require theta[1] > theta[2].")
    }
  }
  
  # ============================================================
  # Simulation containers
  # ============================================================
  n_dose     <- length(tox_prob)
  n_patients <- n_cohort * cohort_size
  
  sel_dose_idx_work <- integer(n_trial)  # index on FINAL working grid (0 if NA)
  sel_dose_std      <- rep(NA_real_, n_trial)
  sel_dose          <- rep(NA_real_, n_trial)
  
  insert_at_cohort <- integer(n_trial)   # first insertion cohort; 0 if none
  n_insertions     <- integer(n_trial)
  early_stop       <- logical(n_trial)
  
  # keep details optionally
  trial_detail <- if (!light_return) vector("list", n_trial) else NULL
  
  # Priors (you used 0.5,0.5)
  pri_alpha <- 0.5
  pri_beta  <- 0.5
  
  # Optional ADM-like cap (comment out if you want unlimited)
  max_insertions_per_trial <- 2L
  
  # ============================================================
  # Main simulation loop
  # ============================================================
  for (trial in seq_len(n_trial)) {
    
    d <- start_dose
    
    # working grid & true tox (simulation-only)
    dose_set_work <- dose_set_std
    tox_prob_work <- tox_prob
    n_dose_work   <- length(dose_set_work)
    
    # trial statistics
    n_dlt_work     <- rep(0, n_dose_work)
    n_treated_work <- rep(0, n_dose_work)
    
    is_eliminated <- rep(FALSE, n_dose_work)
    is_earlystop  <- FALSE
    
    n_insertions_trial <- 0L
    insert_at_cohort[trial] <- 0L
    
    # current posterior (init = prior)
    post_alpha_work <- rep(pri_alpha, n_dose_work)
    post_beta_work  <- rep(pri_beta,  n_dose_work)
    
    # treat sequentially by cohort until total N or early stop
    while (sum(n_treated_work) < n_patients && !is_earlystop) {
      
      # If current dose eliminated, move down to highest non-eliminated; else stop
      if (is_eliminated[d]) {
        cand <- which(!is_eliminated & (seq_len(n_dose_work) < d))
        if (length(cand) == 0) {
          is_earlystop <- TRUE
          break
        } else {
          d <- max(cand)
        }
      }
      
      # per-dose cap (your n_earlystop semantics)
      if (n_treated_work[d] >= n_earlystop) break
      
      # ---------------- Treat one cohort at current dose ----------------
      y <- stats::rbinom(1, cohort_size, tox_prob_work[d])
      n_dlt_work[d]     <- n_dlt_work[d] + y
      n_treated_work[d] <- n_treated_work[d] + cohort_size
      
      # ---------------- Safety elimination (local posterior at dose d) ----------------
      alpha_loc <- pri_alpha + n_dlt_work[d]
      beta_loc  <- pri_beta  + (n_treated_work[d] - n_dlt_work[d])
      p_tox_gt_phi <- 1 - stats::pbeta(target_prob, alpha_loc, beta_loc)
      
      cutoff_here <- cutoff_elimin
      if (extra_safe && d == 1) cutoff_here <- max(cutoff_elimin - offset, 0)
      
      if (p_tox_gt_phi > cutoff_here) {
        is_eliminated[d:n_dose_work] <- TRUE
        if (d == 1) {
          is_earlystop <- TRUE
          break
        } else {
          d <- move_to_safe(d, is_eliminated, dir = "down")
          next
        }
      }
      
      # ---------------- BKP update for ALL doses ----------------
      tmp_post <- bkp_update_all(dose_set_work, n_dlt_work, n_treated_work,
                                 pri_alpha, pri_beta, symmetric, theta)
      post_alpha_work <- tmp_post$post_alpha
      post_beta_work  <- tmp_post$post_beta
      
      # ---------------- Insertion trigger (Eq.(5) logic + optional PAVA) ----------------
      res <- insert_check(
        j = d,
        post_alpha = post_alpha_work,
        post_beta  = post_beta_work,
        key_L = key_L,
        key_U = key_U,
        C1 = C1,
        C2 = C2,
        n_treated = n_treated_work,
        use_pava = TRUE,
        return_probs = FALSE
      )
      
      need_insert   <- res$need_insert
      insert_code   <- res$insert_code
      LOC_BELOW_MIN <- res$LOC_BELOW_MIN
      LOC_ABOVE_MAX <- res$LOC_ABOVE_MAX
      
      # ---------------- Implement insertion ----------------
      if (need_insert &&
          n_insertions_trial < max_insertions_per_trial &&
          sum(n_treated_work) <= (n_patients - cohort_size)) {
        
        tol <- 1e-12
        inserted_any <- FALSE
        d_old <- d
        
        # =========================================================
        # ADM-style BELOW-MIN: insert half, treat, then choose again
        # =========================================================
        if (insert_code == LOC_BELOW_MIN) {
          
          if (abs(dose_set_work[d_old] - dose_set_work[1]) < tol &&
              sum(n_treated_work) <= (n_patients - 2 * cohort_size) &&
              n_insertions_trial <= (max_insertions_per_trial - 1L)) {
            
            halfdose <- dose_set_work[1] / 2
            
            if (is.finite(halfdose) && halfdose > tol &&
                !any(abs(dose_set_work - halfdose) < tol)) {
              
              tox_half <- ins_true_prob(
                insert_code = LOC_BELOW_MIN,
                d = d_old,
                p = tox_prob_work,
                target = target_prob,
                LOC_BELOW_MIN = LOC_BELOW_MIN,
                LOC_ABOVE_MAX = LOC_ABOVE_MAX
              )
              
              tmp <- sort(c(dose_set_work, halfdose), index.return = TRUE)
              dose_set_work <- tmp$x
              
              tox_prob_work  <- c(tox_prob_work, tox_half)[tmp$ix]
              n_dlt_work     <- c(n_dlt_work, 0)[tmp$ix]
              n_treated_work <- c(n_treated_work, 0)[tmp$ix]
              is_eliminated  <- c(is_eliminated, FALSE)[tmp$ix]
              post_alpha_work <- c(post_alpha_work, pri_alpha)[tmp$ix]
              post_beta_work  <- c(post_beta_work,  pri_beta)[tmp$ix]
              
              n_dose_work <- n_dose_work + 1L
              inserted_any <- TRUE
              
              if (insert_at_cohort[trial] == 0L) {
                insert_at_cohort[trial] <- sum(n_treated_work) / cohort_size
              }
              
              # treat one cohort at halfdose
              d_half <- which.min(abs(dose_set_work - halfdose))
              y1 <- stats::rbinom(1, cohort_size, tox_prob_work[d_half])
              n_dlt_work[d_half]     <- n_dlt_work[d_half] + y1
              n_treated_work[d_half] <- n_treated_work[d_half] + cohort_size
              
              n_insertions_trial <- n_insertions_trial + 1L
              
              # Step 2: choose again within (halfdose, original d1) -> now (dose_set_work[1], dose_set_work[2])
              if (sum(n_treated_work) <= (n_patients - cohort_size) &&
                  n_insertions_trial < max_insertions_per_trial) {
                
                sel2 <- choose_newdose(
                  dl = dose_set_work[1],
                  dr = dose_set_work[2],
                  dose_set_work = dose_set_work,
                  n_dlt_work = n_dlt_work,
                  n_treated_work = n_treated_work,
                  pri_alpha = pri_alpha,
                  pri_beta  = pri_beta,
                  key_L = key_L,
                  key_U = key_U,
                  symmetric = symmetric,
                  theta = theta,
                  M = 200
                )
                newdose2 <- sel2$newdose
                
                if (is.finite(newdose2) &&
                    (newdose2 > dose_set_work[1] + tol) && (newdose2 < dose_set_work[2] - tol) &&
                    !any(abs(dose_set_work - newdose2) < tol)) {
                  
                  # BETWEEN insertion in interval i=1 (between positions 1 and 2)
                  tox_new2 <- ins_true_prob(
                    insert_code = 1L,
                    d = 1L,
                    p = tox_prob_work,   # aligned with current working grid after half insertion
                    target = target_prob,
                    LOC_BELOW_MIN = LOC_BELOW_MIN,
                    LOC_ABOVE_MAX = LOC_ABOVE_MAX
                  )
                  
                  tmp2 <- sort(c(dose_set_work, newdose2), index.return = TRUE)
                  dose_set_work <- tmp2$x
                  
                  tox_prob_work  <- c(tox_prob_work, tox_new2)[tmp2$ix]
                  n_dlt_work     <- c(n_dlt_work, 0)[tmp2$ix]
                  n_treated_work <- c(n_treated_work, 0)[tmp2$ix]
                  is_eliminated  <- c(is_eliminated, FALSE)[tmp2$ix]
                  post_alpha_work <- c(post_alpha_work, pri_alpha)[tmp2$ix]
                  post_beta_work  <- c(post_beta_work,  pri_beta)[tmp2$ix]
                  
                  n_dose_work <- n_dose_work + 1L
                  inserted_any <- TRUE
                  
                  # treat one cohort at second inserted dose
                  d2 <- which.min(abs(dose_set_work - newdose2))
                  y2 <- stats::rbinom(1, cohort_size, tox_prob_work[d2])
                  n_dlt_work[d2]     <- n_dlt_work[d2] + y2
                  n_treated_work[d2] <- n_treated_work[d2] + cohort_size
                  
                  n_insertions_trial <- n_insertions_trial + 1L
                  
                  d <- d2
                } else {
                  d <- d_half
                }
              } else {
                d <- d_half
              }
            }
          }
          
        } else {
          
          # =========================================================
          # Metric one-shot for BETWEEN or ABOVE-MAX: choose_newdose once
          # =========================================================
          if (insert_code == LOC_ABOVE_MAX) {
            dl <- dose_set_work[n_dose_work]
            
            # ADM-style cap: dr corresponds to min(upper bound, 1.5 * current max orig)
            range_len <- (dose_range[2] - dose_range[1])
            dmax_orig <- dose_range[1] + dl * range_len
            upper_orig <- min(dose_range[2], 1.5 * dmax_orig)
            dr <- (upper_orig - dose_range[1]) / range_len
          } else {
            i <- as.integer(insert_code)
            if (i < 1 || i >= n_dose_work) {
              dl <- dr <- NA_real_
            } else {
              dl <- dose_set_work[i]
              dr <- dose_set_work[i + 1]
            }
          }
          
          if (is.finite(dl) && is.finite(dr) && (dr - dl) > 1e-12) {
            
            sel <- choose_newdose(
              dl = dl, dr = dr,
              dose_set_work = dose_set_work,
              n_dlt_work = n_dlt_work,
              n_treated_work = n_treated_work,
              pri_alpha = pri_alpha,
              pri_beta  = pri_beta,
              key_L = key_L,
              key_U = key_U,
              symmetric = symmetric,
              theta = theta,
              M = 200
            )
            newdose <- sel$newdose
            
            if (is.finite(newdose) &&
                (newdose > dl + tol) && (newdose < dr - tol) &&
                !any(abs(dose_set_work - newdose) < tol)) {
              
              tox_new <- ins_true_prob(
                insert_code = insert_code,
                d = d_old,
                p = tox_prob_work,
                target = target_prob,
                LOC_BELOW_MIN = LOC_BELOW_MIN,
                LOC_ABOVE_MAX = LOC_ABOVE_MAX
              )
              
              tmp <- sort(c(dose_set_work, newdose), index.return = TRUE)
              dose_set_work <- tmp$x
              
              tox_prob_work  <- c(tox_prob_work, tox_new)[tmp$ix]
              n_dlt_work     <- c(n_dlt_work, 0)[tmp$ix]
              n_treated_work <- c(n_treated_work, 0)[tmp$ix]
              is_eliminated  <- c(is_eliminated, FALSE)[tmp$ix]
              post_alpha_work <- c(post_alpha_work, pri_alpha)[tmp$ix]
              post_beta_work  <- c(post_beta_work,  pri_beta)[tmp$ix]
              
              n_dose_work <- n_dose_work + 1L
              inserted_any <- TRUE
              
              if (insert_at_cohort[trial] == 0L) {
                insert_at_cohort[trial] <- sum(n_treated_work) / cohort_size
              }
              
              d_new <- which.min(abs(dose_set_work - newdose))
              y_new <- stats::rbinom(1, cohort_size, tox_prob_work[d_new])
              n_dlt_work[d_new]     <- n_dlt_work[d_new] + y_new
              n_treated_work[d_new] <- n_treated_work[d_new] + cohort_size
              
              n_insertions_trial <- n_insertions_trial + 1L
              d <- d_new
            }
          }
        }
        
        # After any insertion/treatment, recompute BKP posterior for ALL doses
        if (inserted_any) {
          tmp_post2 <- bkp_update_all(dose_set_work, n_dlt_work, n_treated_work,
                                      pri_alpha, pri_beta, symmetric, theta)
          post_alpha_work <- tmp_post2$post_alpha
          post_beta_work  <- tmp_post2$post_beta
        }
      }
      
      # ---------------- Decide next dose (SKBD-style; based on BKP posterior at current d) ----------------
      if (sum(n_treated_work) >= n_patients) break
      if (is_earlystop) break
      
      action <- decide_next_dose(d, post_alpha_work, post_beta_work, key_L, key_U)
      
      if (action == "E") {
        d_next <- move_to_safe(d, is_eliminated, dir = "up")
      } else if (action == "D") {
        d_next <- move_to_safe(d, is_eliminated, dir = "down")
      } else {
        d_next <- d
      }
      
      d <- d_next
    } # end while
    
    # ---------------- Select final MTD (argmax target probability under BKP posterior) ----------------
    early_stop[trial] <- is_earlystop
    n_insertions[trial] <- n_insertions_trial
    
    if (!is_earlystop) {
      # ensure posterior is up-to-date
      tmp_post_end <- bkp_update_all(dose_set_work, n_dlt_work, n_treated_work,
                                     pri_alpha, pri_beta, symmetric, theta)
      a_end <- tmp_post_end$post_alpha
      b_end <- tmp_post_end$post_beta
      
      prob_target_end <- stats::pbeta(key_U, a_end, b_end) - stats::pbeta(key_L, a_end, b_end)
      
      admissible <- (!is_eliminated) & (n_treated_work > 0)
      if (any(admissible)) {
        idx <- which.max(ifelse(admissible, prob_target_end, -Inf))
        sel_dose_idx_work[trial] <- idx
        sel_dose_std[trial] <- dose_set_work[idx]
        sel_dose[trial] <- dose_range[1] + sel_dose_std[trial] * (dose_range[2] - dose_range[1])
      } else {
        sel_dose_idx_work[trial] <- 0L
        sel_dose_std[trial] <- NA_real_
        sel_dose[trial] <- NA_real_
      }
    } else {
      sel_dose_idx_work[trial] <- 0L
      sel_dose_std[trial] <- NA_real_
      sel_dose[trial] <- NA_real_
    }
    
    if (!light_return) {
      trial_detail[[trial]] <- list(
        dose_set_work_std = dose_set_work,
        dose_set_work     = dose_range[1] + dose_set_work * (dose_range[2] - dose_range[1]),
        tox_prob_work     = tox_prob_work,
        n_dlt_work        = n_dlt_work,
        n_treated_work    = n_treated_work,
        is_eliminated     = is_eliminated,
        post_alpha_work   = post_alpha_work,
        post_beta_work    = post_beta_work,
        early_stop        = is_earlystop,
        n_insertions      = n_insertions_trial,
        insert_at_cohort  = insert_at_cohort[trial]
      )
    }
  } # end trial loop
  
  # ============================================================
  # Summary outputs
  # ============================================================
  out <- list(
    # per-trial
    sel_dose_idx_work = sel_dose_idx_work,
    sel_dose_std      = sel_dose_std,
    sel_dose          = sel_dose,
    early_stop        = early_stop,
    n_insertions      = n_insertions,
    insert_at_cohort  = insert_at_cohort,
    
    # summary
    early_stop_rate   = mean(early_stop),
    mean_insertions   = mean(n_insertions),
    pct_any_insertion = mean(n_insertions > 0)
  )
  
  if (!light_return) out$trial_detail <- trial_detail
  
  return(out)
}
