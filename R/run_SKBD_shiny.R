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
#'   \item \strong{Trial Setting:} interactively generate dose-escalation and
#'   de-escalation boundary tables.
#'   \item \strong{Simulation:} run batch simulations under user-specified
#'   scenarios and summarize the operating characteristics.
#' }
#'
#' By default, the app follows the standard behavior of [shiny::runApp()]. In
#' RStudio, this typically opens the app in the RStudio Viewer or Shiny window;
#' in other R sessions, it may open in the system default browser.
#'
#' @param launch.browser Logical or function; controls whether and how the app is
#'   opened after launch. The default is
#'   `getOption("shiny.launch.browser", interactive())`, which follows the
#'   standard behavior of [shiny::runApp()]. In RStudio, this typically opens the
#'   app in the RStudio Viewer or Shiny window. Set `launch.browser = TRUE` to
#'   open the app in the system default browser, or `launch.browser = FALSE` to
#'   start the app without opening a browser.
#' @param ... Additional arguments passed to [shiny::runApp()], such as `port`,
#'   `host`, or `display.mode`.
#'
#' @return
#' Invisibly returns the value from [shiny::runApp()].
#'
#' @examples
#' \dontrun{
#' run_SKBD_shiny()
#'
#' # Open in the system default browser
#' run_SKBD_shiny(launch.browser = TRUE)
#'
#' # Start the app without opening a browser
#' run_SKBD_shiny(launch.browser = FALSE)
#' }
#'
#' @export
run_SKBD_shiny <- function(
    launch.browser = getOption("shiny.launch.browser", interactive()),
    ...
) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop(
      "Package 'shiny' is required for this function. Please install it with install.packages('shiny').",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop(
      "Package 'DT' is required for this function. Please install it with install.packages('DT').",
      call. = FALSE
    )
  }
  
  app_dir <- system.file("app", package = "SKBD")
  
  if (!nzchar(app_dir)) {
    stop(
      "Cannot find the Shiny app directory in the installed SKBD package.",
      call. = FALSE
    )
  }
  
  invisible(
    shiny::runApp(
      appDir = app_dir,
      launch.browser = launch.browser,
      ...
    )
  )
}
