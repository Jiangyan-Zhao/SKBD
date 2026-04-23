#' @title Launch an Interactive Shiny App for SKBD
#'
#' @description
#' `run_SKBD_shiny()` launches a Shiny interface for exploring key functions in
#' the `SKBD` package, including decision-boundary generation and operating-
#' characteristic simulation.
#'
#' @details
#' The app currently provides two tabs:
#' \itemize{
#'   \item \strong{Boundary (SKBD):} interactively compute
#'   `get_boundary_SKBD()` and inspect escalation/de-escalation/elimination
#'   boundaries.
#'   \item \strong{OC Simulation (SKBD):} interactively run `get_OC_SKBD()`
#'   and inspect summary metrics.
#' }
#'
#' @param launch.browser Logical; passed to [shiny::runApp()].
#'
#' @return
#' Invisibly returns the value from [shiny::runApp()].
#'
#' @examples
#' \dontrun{
#' run_SKBD_shiny()
#' }
#'
#' @export
run_SKBD_shiny <- function(launch.browser = TRUE) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run this app. Please install it first.")
  }

  app_dir <- system.file("app", package = "SKBD")
  if (app_dir == "") {
    stop("Cannot find Shiny app directory in installed package.")
  }

  invisible(shiny::runApp(appDir = app_dir, launch.browser = launch.browser))
}
