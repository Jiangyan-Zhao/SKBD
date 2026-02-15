#' @title Pre-tabulated Decision Tables for the Shared Keyboard Design (SKBD)
#'
#' @description
#' \code{get_boundary_SKBD()} generates pre-tabulated dose escalation/de-escalation
#' (and overdose elimination) decision tables for the Shared Keyboard Design (SKBD)
#' at a \strong{current dose} \code{d}, given the accumulated data across all doses.
#' The output is a keyboard-like boundary table (similar to the original Keyboard design),
#' but the decisions are based on \emph{borrowed information} through a kernel-weighted
#' Beta posterior at the current dose.
#'
#' @details
#' The SKBD uses a kernel function to borrow information across dose levels. At the
#' current dose \code{d}, define kernel weights \eqn{w_s(d)} for each dose \eqn{s}.
#' Let \eqn{y_s} and \eqn{n_s} be the observed number of DLTs and the number treated
#' at dose \eqn{s}. The borrowed (effective) Beta posterior parameters are
#' \deqn{
#' \alpha_{\text{post}} = \alpha_0 + \sum_s w_s(d)\,y_s, \qquad
#' \beta_{\text{post}}  = \beta_0  + \sum_s w_s(d)\,(n_s-y_s),
#' }
#' where \eqn{(\alpha_0,\beta_0)=(1,1)} by default (non-informative prior).
#'
#' Given \eqn{(\alpha_{\text{post}},\beta_{\text{post}})}, the strongest keyboard
#' interval ("strongest key") is the interval with the largest posterior probability
#' mass under \eqn{\text{Beta}(\alpha_{\text{post}},\beta_{\text{post}})}. The action
#' is then determined by comparing the strongest key to the target key:
#' \itemize{
#'   \item \code{"E"}: escalate if the strongest key is to the left of the target key;
#'   \item \code{"S"}: stay if the strongest key is the target key;
#'   \item \code{"D"}: de-escalate if the strongest key is to the right of the target key.
#' }
#'
#' For patient safety, an overdose control rule is applied:
#' if at least 3 patients have been treated at the current dose and
#' \eqn{\Pr(p_d > \phi \mid \text{data}) > \text{cutoff\_elimin}}, the current dose is
#' eliminated and labeled as \code{"DU"} in the decision table.
#'
#' \strong{Two table modes.}
#' This function supports two ways to construct a keyboard-like table:
#' \itemize{
#'   \item \code{table_type = "baseline"}: ignore the current dose's existing history
#'         (\code{n[d]}, \code{y[d]}) when tabulating; only other doses contribute as fixed
#'         borrowed information. This is useful for creating a "baseline" conditional
#'         table given other-dose information.
#'   \item \code{table_type = "continue"}: condition on the current dose's existing
#'         history (\code{n[d]}, \code{y[d]}) and tabulate decisions \emph{from now on}.
#'         Only \eqn{n \ge n[d]} columns are relevant, and (optionally) past columns can be
#'         hidden via \code{show_past = FALSE}.
#' }
#'
#' @param target_prob Target toxicity (DLT) probability \eqn{\phi} used to define the target key.
#' @param d Integer index of the current dose level (between 1 and \code{length(y)}).
#' @param y Integer vector of length \eqn{J} giving the number of observed DLTs at each dose.
#' @param n Integer vector of length \eqn{J} giving the number of treated patients at each dose.
#' @param n_cohort Total number of cohorts in the trial (used to define the maximal table size).
#' @param cohort_size Number of patients per cohort (used for cohort-aligned \code{boundary_tab} output).
#' @param symmetric Logical; if \code{TRUE}, use a symmetric kernel; otherwise use an asymmetric kernel.
#' @param k_left Numeric scalar in `(0,1)`. Left-side neighbor borrowing strength passed to `kernel_fun()`.
#' @param k_right Numeric scalar in `(0,1)`. Right-side neighbor borrowing strength passed to `kernel_fun()`
#'        (and used for symmetric borrowing when \code{symmetric=TRUE}).
#' @param ref_gap Optional positive scalar. Reference spacing passed to `kernel_fun()`. If `NULL`,
#'        kernel defaults to the minimum adjacent spacing in `dose_set`.
#' @param dose_set Numeric vector of dose labels (default \code{1:length(y)}). Used only to compute kernel distances.
#' @param n_earlystop Maximum number of patients to display in output tables (columns are truncated at this value).
#' @param margin_left Left margin of the target key (lower bound is \code{target_prob - margin_left}).
#' @param margin_right Right margin of the target key (upper bound is \code{target_prob + margin_right}).
#' @param cutoff_elimin Overdose elimination cutoff; eliminate if \eqn{\Pr(p_d > \phi \mid \text{data}) >} this value.
#' @param extra_safe Logical; if \code{TRUE} and \code{d==1}, return an additional conservative stopping boundary
#'        (similar to the BOIN extra-safe rule).
#' @param offset Nonnegative offset in (0, 0.5); used only when \code{extra_safe=TRUE} to tighten the stopping rule.
#' @param table_type Character; either \code{"baseline"} or \code{"continue"} (see Details).
#' @param show_past Logical; if \code{FALSE} and \code{table_type="continue"}, hide columns before the current
#'        treated count \code{n[d]} (or \code{n[d]+1} depending on \code{start_from}).
#' @param start_from Character; only relevant when \code{show_past=FALSE} and \code{table_type="continue"}.
#'        \code{"current"} starts from column \code{n[d]}; \code{"next"} starts from column \code{n[d]+1}.
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item \code{boundary_tab}: a 4-row boundary table (cohort-aligned columns) containing escalation,
#'         de-escalation, and elimination boundaries.
#'   \item \code{full_boundary_tab}: the full 4-row boundary table (possibly truncated/hiding past columns).
#'   \item \code{decision_table}: the underlying decision table with entries \code{"E"}, \code{"S"}, \code{"D"}, \code{"DU"}.
#'   \item \code{weight}: the normalized kernel weights used for borrowing at dose \code{d}.
#'   \item \code{meta}: a list recording \code{table_type}, \code{show_past}, \code{start_from}, and the starting column used.
#' }
#' If \code{extra_safe=TRUE} and \code{d==1}, additional elements \code{cutoff} and \code{stop_boundary} are returned.
#'
#' @seealso
#' \code{\link{get_Key}}, \code{\link{get_strongKey}}, \code{\link{kernel_fun}}
#'
#' @examples
#' ## Example data across 5 doses:
#' y <- c(0, 1, 2, 2, 0)
#' n <- c(3, 6, 9, 3, 0)
#' d <- 3
#' target_prob <- 0.3
#'
#' ## 1) Baseline conditional table:
#' ##    Ignore (n[d], y[d]) at the current dose when tabulating
#' out_baseline <- get_boundary_SKBD(target_prob = target_prob,
#'                                  d = d, y = y, n = n,
#'                                  table_type = "baseline")
#' out_baseline$full_boundary_tab
#'
#' ## 2) Continue table:
#' ##    Condition on (n[d], y[d]) and tabulate decisions from now on
#' out_continue <- get_boundary_SKBD(target_prob = target_prob,
#'                                  d = d, y = y, n = n,
#'                                  table_type = "continue")
#' out_continue$full_boundary_tab
#'
#' ## 3) Continue table (display only future columns):
#' out_future <- get_boundary_SKBD(target_prob = target_prob,
#'                                d = d, y = y, n = n,
#'                                table_type = "continue",
#'                                show_past = FALSE,
#'                                start_from = "next")
#' out_future$full_boundary_tab
#'
#' @export

get_boundary_SKBD <- function(target_prob, d,
                              y, n,
                              n_cohort = 10, cohort_size = 3,
                              symmetric = FALSE,
                              k_left = 0.2, k_right = 0.8, ref_gap = NULL,
                              dose_set = 1:length(y),
                              n_earlystop = 1000,
                              margin_left = 0.05, margin_right = 0.05,
                              cutoff_elimin = 0.95,
                              extra_safe = FALSE, offset = 0.05,
                              table_type = c("baseline", "continue"),
                              show_past = TRUE,
                              start_from = c("current", "next")) {
  ## ---------------------------
  ## Match arguments
  ## ---------------------------
  table_type <- match.arg(table_type)
  start_from <- match.arg(start_from)
  
  ## ---------------------------
  ## Basic setup and constants
  ## ---------------------------
  n_dose    <- length(y)                        # number of dose levels
  n_patient <- n_cohort * cohort_size           # maximum total sample size in the trial
  
  ## Prior for toxicity probability at the current dose after borrowing:
  ## Beta(alpha_Pri, beta_Pri); default is non-informative Beta(1,1)
  alpha_Pri <- 1
  beta_Pri  <- 1
  
  ## ---------------------------
  ## Construct keyboard "keys"
  ## ---------------------------
  keys <- get_Key(target_prob, margin_left, margin_right)
  target_key <- which(keys == (target_prob - margin_left))
  
  ## ---------------------------
  ## Standardize doses and set kernel parameters
  ## ---------------------------
  if (max(dose_set) == min(dose_set)) stop("dose_set must have at least two distinct values.")
  dose_set_std <- (dose_set - min(dose_set)) / (max(dose_set) - min(dose_set))
  
  ker_vals <- kernel_fun(
    dose = dose_set_std[d],
    dose_set = dose_set_std,
    symmetric = symmetric,
    k_left = k_left,
    k_right = k_right,
    ref_gap = ref_gap
  )
  
  ## ---------------------------
  ## Define borrowing weights (normalize over observed doses + current dose)
  ## ---------------------------
  obs_mask <- (n > 0)
  obs_mask[d] <- TRUE
  
  w_raw <- ker_vals
  w_raw[!obs_mask] <- 0
  if (sum(w_raw) == 0) stop("All weights are zero; check kernel k_left/k_right settings.")
  weight <- w_raw / sum(w_raw)
  
  ## ---------------------------
  ## Prepare baseline information depending on table_type
  ## ---------------------------
  n0 <- n[d]
  y0 <- y[d]
  
  if (table_type == "baseline") {
    ## Baseline table: ignore current-dose history, only other doses contribute fixed parts
    y_fix <- y; n_fix <- n
    y_fix[d] <- 0
    n_fix[d] <- 0
    base_n0 <- 0
    base_y0 <- 0
  } else {
    ## Continue table: condition on current-dose history as baseline
    y_fix <- y
    n_fix <- n
    base_n0 <- n0
    base_y0 <- y0
  }
  
  fixed_y_part        <- sum(weight * y_fix)                 # sum w_s * y_s (baseline included if continue)
  fixed_nminusy_part  <- sum(weight * (n_fix - y_fix))       # sum w_s * (n_s - y_s)
  w_d <- weight[d]                                          # weight on current dose
  
  ## ---------------------------
  ## Build decision table
  ## ---------------------------
  decision_table <- matrix(NA_character_, nrow = n_patient + 1, ncol = n_patient)
  
  if (table_type == "baseline") {
    
    ## Enumerate total n_d and y_d at current dose from scratch
    for (n_d in 1:n_patient) {
      
      eliminate <- FALSE
      
      for (y_d in 0:n_d) {
        
        ## Posterior with effective counts
        alpha_Post <- alpha_Pri + fixed_y_part + w_d * y_d
        beta_Post  <- beta_Pri  + fixed_nminusy_part + w_d * (n_d - y_d)
        
        ## Elimination only if at least 3 patients at current dose
        if (n_d >= 3 && (1 - pbeta(target_prob, alpha_Post, beta_Post) > cutoff_elimin)) {
          eliminate <- TRUE
          break
        }
        
        strong_key <- get_strongKey(alpha_Post, beta_Post, keys, margin_left, margin_right)
        
        if (strong_key > target_key) {
          decision_table[y_d + 1, n_d] <- "D"
        } else if (strong_key == target_key) {
          decision_table[y_d + 1, n_d] <- "S"
        } else {
          decision_table[y_d + 1, n_d] <- "E"
        }
      }
      
      if (eliminate) {
        decision_table[(y_d + 1):(n_d + 1), n_d] <- rep("DU", n_d - y_d + 1)
      }
    }
    
  } else {
    
    ## Continue table: enumerate additional patients m and additional DLTs t
    ## Total counts: n_tot = n0 + m, y_tot = y0 + t
    for (m in 0:(n_patient - base_n0)) {
      
      n_tot <- base_n0 + m
      eliminate <- FALSE
      
      for (t in 0:m) {
        
        y_tot <- base_y0 + t
        
        alpha_Post <- alpha_Pri + fixed_y_part + w_d * 0 + 0  # kept for readability
        beta_Post  <- beta_Pri  + fixed_nminusy_part + w_d * 0 + 0
        
        ## The fixed parts already include baseline y0,n0 at dose d (through y_fix,n_fix),
        ## so we only add increments (t, m-t) at current dose:
        alpha_Post <- alpha_Post + w_d * t
        beta_Post  <- beta_Post  + w_d * (m - t)
        
        ## Elimination only if total patients at current dose >= 3
        if (n_tot >= 3 && (1 - pbeta(target_prob, alpha_Post, beta_Post) > cutoff_elimin)) {
          eliminate <- TRUE
          y_trigger <- y_tot
          break
        }
        
        strong_key <- get_strongKey(alpha_Post, beta_Post, keys, margin_left, margin_right)
        
        if (strong_key > target_key) {
          decision_table[y_tot + 1, n_tot] <- "D"
        } else if (strong_key == target_key) {
          decision_table[y_tot + 1, n_tot] <- "S"
        } else {
          decision_table[y_tot + 1, n_tot] <- "E"
        }
      }
      
      if (eliminate) {
        ## Mark DU from the first triggering y (y_trigger) upward at this n_tot
        decision_table[(y_trigger + 1):(n_tot + 1), n_tot] <- "DU"
      }
    }
  }
  
  colnames(decision_table) <- 1:n_patient
  rownames(decision_table) <- 0:n_patient
  
  ## ---------------------------
  ## Extract escalation/de-escalation/elimination boundaries
  ## ---------------------------
  boundary <- matrix(NA_real_, nrow = 4, ncol = n_patient)
  boundary[1, ] <- 1:n_patient
  
  for (i in 1:n_patient) {
    
    ## Escalate boundary: largest y that still leads to "E"
    if (length(which(decision_table[, i] == "E"))) {
      boundary[2, i] <- max(which(decision_table[, i] == "E")) - 1
    } else {
      boundary[2, i] <- NA_integer_   # no escalation possible at this n
    }
    
    ## De-escalate boundary: smallest y that leads to "D"; if none, use "DU"
    if (length(which(decision_table[, i] == "D"))) {
      boundary[3, i] <- min(which(decision_table[, i] == "D")) - 1
    } else if (length(which(decision_table[, i] == "DU"))) {
      boundary[3, i] <- min(which(decision_table[, i] == "DU")) - 1
    } else {
      boundary[3, i] <- NA
    }
    
    ## Elimination boundary: smallest y that leads to "DU"
    if (length(which(decision_table[, i] == "DU"))) {
      boundary[4, i] <- min(which(decision_table[, i] == "DU")) - 1
    } else {
      boundary[4, i] <- NA
    }
  }
  
  # At least one patient exists in the current dose level,
  # if there are patients enrolled in the higher dose level. 
  if(d < n_dose && n[d+1] > 0){
    boundary[2:4, 1] <- NA
  }
  
  colnames(boundary) <- rep("", n_patient)
  rownames(boundary) <- c("Number of patients treated",
                          "Escalate if # of DLT <=",
                          "de-escalate if # of DLT >=",
                          "Eliminate if # of DLT >=")
  
  ## ---------------------------
  ## Optionally hide past columns for display
  ## ---------------------------
  start_n <- 1
  if (!show_past && table_type == "continue" && base_n0 > 0) {
    start_n <- if (start_from == "current") base_n0 else min(n_patient, base_n0 + 1)
  }
  
  keep_cols_all <- start_n:min(n_patient, n_earlystop)
  
  ## Cohort-aligned columns (multiples of cohort_size), intersect with keep_cols_all
  cohort_cols <- (1:floor(min(n_patient, n_earlystop) / cohort_size)) * cohort_size
  cohort_cols <- cohort_cols[cohort_cols >= start_n]
  
  out <- list(
    boundary_tab      = boundary[, cohort_cols, drop = FALSE],
    full_boundary_tab = boundary[, keep_cols_all, drop = FALSE],
    decision_table    = decision_table[, keep_cols_all, drop = FALSE],
    weight            = weight,
    meta              = list(table_type = table_type,
                             show_past = show_past,
                             start_from = start_from,
                             current_n = n0,
                             current_y = y0,
                             start_n = start_n)
  )
  
  ## ---------------------------
  ## Optional: extra-safe stopping boundary (only for the lowest dose)
  ## ---------------------------
  if (extra_safe && d == 1) {
    
    stopbd <- NULL
    n_dt   <- NULL
    
    ## For extra-safe rule, it is most meaningful in "continue" mode at lowest dose.
    ## We still compute it generically using the same baseline convention as table_type.
    for (n_d in 1:n_patient) {
      
      n_dt <- c(n_dt, n_d)
      
      if (n_d < 3) {
        stopbd <- c(stopbd, NA)
      } else {
        
        stopneed <- 0
        y_hit <- NA
        
        ## In baseline mode: y_d ranges 0..n_d
        ## In continue mode at d=1: interpret n_d as total at dose 1 (includes baseline)
        if (table_type == "baseline") {
          for (y_d in 0:n_d) {
            alpha_Post <- alpha_Pri + fixed_y_part + w_d * y_d
            beta_Post  <- beta_Pri  + fixed_nminusy_part + w_d * (n_d - y_d)
            
            if (1 - pbeta(target_prob, alpha_Post, beta_Post) > cutoff_elimin - offset) {
              stopneed <- 1
              y_hit <- y_d
              break
            }
          }
        } else {
          ## continue mode: let total at dose be n_d, so increments are:
          ## m = n_d - n0, t = y_d - y0 (only feasible when n_d>=n0 and y_d>=y0)
          if (n_d < base_n0) {
            stopneed <- 0
          } else {
            m <- n_d - base_n0
            for (t in 0:m) {
              y_tot <- base_y0 + t
              
              alpha_Post <- alpha_Pri + fixed_y_part + w_d * t
              beta_Post  <- beta_Pri  + fixed_nminusy_part + w_d * (m - t)
              
              if (1 - pbeta(target_prob, alpha_Post, beta_Post) > cutoff_elimin - offset) {
                stopneed <- 1
                y_hit <- y_tot
                break
              }
            }
          }
        }
        
        stopbd <- c(stopbd, if (stopneed == 1) y_hit else NA)
      }
    }
    
    stopboundary <- rbind(n_dt, stopbd)[, keep_cols_all, drop = FALSE]
    rownames(stopboundary) <- c("The number of patients treated at the lowest dose",
                                "Stop the trial if # of DLT >= ")
    colnames(stopboundary) <- rep("", ncol(stopboundary))
    
    out <- c(out, list(target_prob = target_prob,
                       cutoff = cutoff_elimin - offset,
                       stop_boundary = stopboundary))
  }
  
  return(out)
}
