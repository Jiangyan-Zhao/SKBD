"_PACKAGE"

#' @name SKBD-package
#'
#' @title Shared Keyboard Designs for Phase I Dose-Finding Trials
#'
#' @description The \pkg{SKBD} package provides tools for implementing the
#'   Shared Keyboard Design (SKBD) for model-assisted phase I dose-finding
#'   trials. SKBD extends the Keyboard design by incorporating kernel-weighted
#'   information sharing across dose levels while preserving the transparent
#'   key-based dose-escalation and de-escalation structure.
#'
#'   The package supports decision-boundary construction, operating-characteristic
#'   simulation, adaptive dose insertion, and time-to-event extensions for
#'   late-onset toxicity outcomes. It also provides utilities for generating
#'   random monotone dose-toxicity scenarios and an interactive Shiny interface
#'   for trial-planning workflows.
#'
#' @section Main Functions: Core functionality is organized as follows:
#' \describe{
#'   \item{\code{\link{get_boundary_SKBD}}}{
#'     Generate pre-tabulated dose-escalation, de-escalation, and overdose
#'     elimination boundaries for the Shared Keyboard Design, conditional on
#'     accumulated trial data across dose levels.
#'   }
#'   \item{\code{\link{get_OC_SKBD}}}{
#'     Simulate operating characteristics for the standard SKBD under fixed
#'     dose-toxicity scenarios, including selection accuracy, patient allocation,
#'     overdose risk, and monotonicity diagnostics.
#'   }
#'   \item{\code{\link{get_OC_TITE_SKBD}}}{
#'     Simulate operating characteristics for the time-to-event extension of
#'     SKBD, which accommodates delayed toxicity outcomes through weighted
#'     partial follow-up information.
#'   }
#'   \item{\code{\link{get_OC_Insert_SKBD}}}{
#'     Simulate operating characteristics for the dose-insertion extension of
#'     SKBD, where the working dose grid can be adaptively refined during the
#'     trial.
#'   }
#'   \item{\code{\link{PUA}}}{
#'     Generate random monotone dose-toxicity scenarios using a pseudo-uniform
#'     algorithm for simulation studies.
#'   }
#'   \item{\code{\link{run_SKBD_shiny}}}{
#'     Launch an interactive Shiny application for generating SKBD decision
#'     tables and running operating-characteristic simulations.
#'   }
#' }
#'
#' @references
#' Zhao J, Shi X, Xu J (2026). Shared Keyboard: An improved Bayesian design 
#'   for phase I clinical trials via Beta kernel process. \emph{ArXiv}. 
#'   https://arxiv.org/abs/2605.25043
#'
#' Yan F, Mandrekar SJ, Yuan Y (2017). Keyboard: A Novel Bayesian Toxicity
#'   Probability Interval Design for Phase I Clinical Trials.
#'   \emph{Clinical Cancer Research}, 23(15), 3994--4003.
#'
#' Lin R, Yuan Y (2020). Time-to-event model-assisted designs for dose-finding
#'   trials with delayed toxicity. \emph{Biostatistics}, 21(4), 807--824.
#'
#' Chu Y, Pan H, Yuan Y (2016). Adaptive dose modification for phase I clinical
#'   trials. \emph{Statistics in Medicine}, 35(20), 3497--3508.
#'
#' Clertant M, O'Quigley J (2017). Semiparametric dose finding methods.
#'   \emph{Journal of the Royal Statistical Society: Series B}, 79(5),
#'   1487--1508.
#' 
#' @importFrom shiny runApp
#' @importFrom stats median pbeta rbeta rbinom rexp runif sd
#' @importFrom utils tail
NULL