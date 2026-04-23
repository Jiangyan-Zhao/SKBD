#' @title Generate Random Monotone Toxicity Scenarios by the Pseudo-Uniform Algorithm
#'
#' @description
#' Generate random monotone dose-toxicity scenarios using a pseudo-uniform
#' algorithm (PUA) for simulation studies in phase I dose-finding.
#'
#' @details
#' This function generates \code{n_scenarios} monotone toxicity curves over a
#' prespecified dose set. For each scenario, it first randomly selects a dose
#' level to serve as the target-dose location, then samples an upper bound for
#' the toxicity range, and finally draws a sorted vector of toxicity
#' probabilities from a uniform distribution on that range.
#'
#' A generated scenario is accepted only if:
#' \itemize{
#'   \item the selected dose is the closest dose to the target toxicity rate
#'   \code{target_prob}, and
#'   \item the toxicity probability at that dose lies within the target key
#'   \eqn{(\phi - \epsilon_1, \phi + \epsilon_2)}, and
#'   \item the adjacent gaps around that dose, when applicable, are both greater
#'   than 0.05 and smaller than 0.3.
#' }
#'
#' The resulting scenarios therefore have a unique and reasonably well-separated
#' target dose, making them suitable for benchmarking dose-finding designs under
#' random monotone settings.
#'
#' @param dose_set Numeric vector of prespecified dose levels. Only its length is
#' used by the current implementation.
#' @param target_prob Scalar in \eqn{(0,1)} giving the target toxicity
#' probability \eqn{\phi}.
#' @param n_scenarios Positive integer. Number of random toxicity scenarios to
#' generate.
#' @param margin_left Nonnegative scalar. Left margin \eqn{\epsilon_1} defining
#' the target key \eqn{(\phi-\epsilon_1, \phi+\epsilon_2)}.
#' @param margin_right Nonnegative scalar. Right margin \eqn{\epsilon_2} defining
#' the target key \eqn{(\phi-\epsilon_1, \phi+\epsilon_2)}.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return
#' A numeric matrix with \code{n_scenarios} rows and \code{length(dose_set)}
#' columns. Each row is an accepted monotone toxicity scenario, with entries
#' giving the true toxicity probabilities at the corresponding dose levels.
#'
#' @references
#' Clertant, M. and O'Quigley, J. (2017).
#' Semiparametric model-assisted designs for phase I dose-finding clinical trials.
#' \emph{Statistics in Medicine}, \bold{36}(30), 4638--4649.
#'
#' Zhou, Y., Lee, J. J., and Yuan, Y. (2021).
#' A utility-based Bayesian optimal interval design for phase I/II dose-finding trials.
#' \emph{Pharmaceutical Statistics}, \bold{20}(1), 125--136.
#'
#' @examples
#' set.seed(1)
#' scen <- PUA(
#'   dose_set = 1:5,
#'   target_prob = 0.30,
#'   n_scenarios = 10
#' )
#' scen
#'
#' @importFrom stats rbeta runif
#' @export

PUA = function(dose_set, target_prob, n_scenarios = 1000, 
               margin_left = 0.05, margin_right = 0.05, seed = 6){
  set.seed(seed)
  n_dose = length(dose_set)
  tox_prob = matrix(NA, nrow = n_scenarios, ncol = n_dose)
  
  # browser()
  s = 1
  while(s <= n_scenarios){
    ## Step 1: select j uniformly
    d = sample(1:n_dose, 1)
    
    ## Step 2: sample M ~ Beta(max{J-j, 0.5}, 1)
    M = rbeta(1, max(n_dose - d, 0.5), 1)
    
    ## Step 3: rejection sampling
    B = target_prob + (1 - target_prob) * M # upper bound B
    
    n_repeat = 0
    repeat {
      n_repeat = n_repeat + 1
      if(n_repeat > 100) break
      
      tox_prob[s, ] = sort(runif(n_dose, 0, B))
      if(d != which.min(abs(tox_prob[s, ] - target_prob))) next
      
      gap_left  = if(d > 1) tox_prob[s, d] - tox_prob[s, d - 1] else NA
      gap_right = if(d < n_dose) tox_prob[s, d + 1] - tox_prob[s, d] else NA
      
      # well-definited MTD, see 21_Zhou_PS_iBOIN
      is_admissible = (
        tox_prob[s, d] >= target_prob - margin_left &&
          tox_prob[s, d] <= target_prob + margin_right &&
          (is.na(gap_left)  || (gap_left  > 0.05 && gap_left  < 0.3)) &&
          (is.na(gap_right) || (gap_right > 0.05 && gap_right < 0.3))
      )
      
      if(is_admissible) {
        s = s + 1
        break
      }
    }
    
    # cat("tox_prob[", s - 1, ", ] = (", sep = "")
    # cat(sprintf("%.3f", tox_prob[s - 1, ]), sep = ", ")
    # cat(")\n", sep = "")
  }
  return(tox_prob)
}