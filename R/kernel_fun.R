#' Asymmetric/symmetric kernel with interpretable neighbor borrowing strengths
#'
#' @description
#' `kernel_fun()` computes kernel similarities between a target dose value `dose`
#' and a set of dose locations `dose_set`, typically on a standardized scale,
#' e.g., `[0, 1]`.
#'
#' Unlike the conventional parameterization that requires users to specify decay
#' parameters `theta` directly, this kernel is parameterized by interpretable
#' neighbor kernel values: `k_left` and `k_right`. These values represent the
#' desired kernel similarity at a reference distance `ref_gap` on the left and
#' right sides, respectively.
#'
#' @details
#' For interior values `0 < k < 1`, the kernel takes the Gaussian form
#' \deqn{
#'   k(d,d_s)=\exp\{-\theta(d-d_s)^2\},
#' }
#' where
#' \deqn{
#'   \theta = -\log(k) / \text{ref\_gap}^2.
#' }
#'
#' Boundary values are handled explicitly:
#' \itemize{
#'   \item `k = 0` means no borrowing from that side, except the target dose
#'         itself, whose kernel weight is always 1.
#'   \item `k = 1` means full borrowing from that side, with no distance decay.
#' }
#'
#' In asymmetric mode (`symmetric = FALSE`), `k_left` controls borrowing from
#' dose locations strictly below `dose`, whereas `k_right` controls borrowing
#' from dose locations strictly above `dose`. The target dose itself always has
#' kernel weight 1.
#'
#' In symmetric mode (`symmetric = TRUE`), a single symmetric kernel is used.
#' By convention, `k_right` is used as the reference similarity.
#'
#' The default `ref_gap` is the minimum adjacent spacing in `dose_set`.
#'
#' @param dose Numeric scalar. Target dose value, usually a standardized dose
#'   value rather than a dose index.
#' @param dose_set Numeric vector. Dose locations at which to evaluate the
#'   kernel similarity.
#' @param symmetric Logical. If `TRUE`, use a single symmetric decay. If
#'   `FALSE`, use separate left/right borrowing strengths.
#' @param k_left Numeric scalar in `[0, 1]`. Desired kernel value at distance
#'   `ref_gap` for dose locations strictly below `dose`. If `NULL`, defaults
#'   to `0.2`.
#' @param k_right Numeric scalar in `[0, 1]`. Desired kernel value at distance
#'   `ref_gap` for dose locations strictly above `dose`. If `NULL`, defaults
#'   to `0.8`.
#' @param ref_gap Numeric scalar. Reference distance used to interpret
#'   `k_left` and `k_right`. If `NULL`, defaults to the minimum adjacent
#'   spacing in `dose_set`.
#'
#' @return Numeric vector of kernel similarities with length `length(dose_set)`.
#'   Entries are in `[0, 1]`.
#'
#' @examples
#' dose_set = c(0.25, 0.50, 0.75)
#'
#' kernel_fun(
#'   dose = 0.50,
#'   dose_set = dose_set,
#'   symmetric = FALSE,
#'   k_left = 0.2,
#'   k_right = 0.8
#' )
#'
#' kernel_fun(
#'   dose = 0.50,
#'   dose_set = dose_set,
#'   symmetric = TRUE,
#'   k_right = 0.5
#' )
#'
#' kernel_fun(
#'   dose = 0.50,
#'   dose_set = dose_set,
#'   symmetric = FALSE,
#'   k_left = 0,
#'   k_right = 1
#' )
#'
#' @noRd
kernel_fun <- function(
    dose,
    dose_set,
    symmetric = FALSE,
    k_left = NULL,
    k_right = NULL,
    ref_gap = NULL
) {
  
  # ---- Validate dose ----
  if (!is.numeric(dose) || length(dose) != 1L || !is.finite(dose)) {
    stop("`dose` must be a finite numeric scalar.", call. = FALSE)
  }
  
  # ---- Validate dose_set ----
  if (!is.numeric(dose_set) ||
      length(dose_set) < 1L ||
      any(!is.finite(dose_set))) {
    stop("`dose_set` must be a finite numeric vector.", call. = FALSE)
  }
  
  dose_set <- as.numeric(dose_set)
  
  # ---- Validate symmetric ----
  if (!is.logical(symmetric) || length(symmetric) != 1L || is.na(symmetric)) {
    stop("`symmetric` must be either TRUE or FALSE.", call. = FALSE)
  }
  
  # ---- Set default neighbor kernel values ----
  if (is.null(k_left)) {
    k_left <- 0.2
  }
  
  if (is.null(k_right)) {
    k_right <- 0.8
  }
  
  # ---- Validate neighbor kernel values ----
  if (!is.numeric(k_left) ||
      length(k_left) != 1L ||
      !is.finite(k_left) ||
      k_left < 0 ||
      k_left > 1) {
    stop("`k_left` must be a scalar in [0, 1].", call. = FALSE)
  }
  
  if (!is.numeric(k_right) ||
      length(k_right) != 1L ||
      !is.finite(k_right) ||
      k_right < 0 ||
      k_right > 1) {
    stop("`k_right` must be a scalar in [0, 1].", call. = FALSE)
  }
  
  # ---- Choose reference distance ----
  if (is.null(ref_gap)) {
    
    dose_set_unique <- sort(unique(dose_set))
    
    if (length(dose_set_unique) < 2L) {
      stop(
        "Need at least two distinct dose points to infer `ref_gap`.",
        call. = FALSE
      )
    }
    
    ref_gap <- min(diff(dose_set_unique))
  }
  
  if (!is.numeric(ref_gap) ||
      length(ref_gap) != 1L ||
      !is.finite(ref_gap) ||
      ref_gap <= 0) {
    stop("`ref_gap` must be a positive finite scalar.", call. = FALSE)
  }
  
  # ---- Compute kernel values ----
  tol <- sqrt(.Machine$double.eps)
  k <- numeric(length(dose_set))
  
  if (isTRUE(symmetric)) {
    
    dist <- abs(dose - dose_set)
    same_id <- dist <= tol
    
    if (k_right == 0) {
      
      k[] <- 0
      k[same_id] <- 1
      
    } else if (k_right == 1) {
      
      k[] <- 1
      
    } else {
      
      theta <- -log(k_right) / (ref_gap^2)
      k <- exp(-theta * dist^2)
      k[same_id] <- 1
    }
    
  } else {
    
    left_id <- dose_set < dose - tol
    same_id <- abs(dose_set - dose) <= tol
    right_id <- dose_set > dose + tol
    
    # Target dose itself always has weight 1.
    k[same_id] <- 1
    
    # Left-side borrowing.
    if (any(left_id)) {
      
      dist_left <- dose - dose_set[left_id]
      
      if (k_left == 0) {
        
        k[left_id] <- 0
        
      } else if (k_left == 1) {
        
        k[left_id] <- 1
        
      } else {
        
        theta1 <- -log(k_left) / (ref_gap^2)
        k[left_id] <- exp(-theta1 * dist_left^2)
      }
    }
    
    # Right-side borrowing.
    if (any(right_id)) {
      
      dist_right <- dose_set[right_id] - dose
      
      if (k_right == 0) {
        
        k[right_id] <- 0
        
      } else if (k_right == 1) {
        
        k[right_id] <- 1
        
      } else {
        
        theta2 <- -log(k_right) / (ref_gap^2)
        k[right_id] <- exp(-theta2 * dist_right^2)
      }
    }
  }
  
  return(k)
}