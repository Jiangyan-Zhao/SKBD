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

  ui <- shiny::fluidPage(
    shiny::titlePanel("SKBD Interactive Explorer"),
    shiny::tabsetPanel(
      shiny::tabPanel(
        title = "Boundary (SKBD)",
        shiny::sidebarLayout(
          shiny::sidebarPanel(
            shiny::numericInput("b_target", "Target toxicity (phi)", value = 0.30, min = 0.01, max = 0.99, step = 0.01),
            shiny::textInput("b_y", "DLTs by dose (comma separated)", value = "0,1,2,2,0"),
            shiny::textInput("b_n", "Treated by dose (comma separated)", value = "3,6,9,3,0"),
            shiny::numericInput("b_d", "Current dose index d", value = 3, min = 1, step = 1),
            shiny::selectInput("b_type", "Table type", choices = c("continue", "baseline"), selected = "continue"),
            shiny::actionButton("b_run", "Run boundary")
          ),
          shiny::mainPanel(
            shiny::verbatimTextOutput("b_msg"),
            shiny::tableOutput("b_table")
          )
        )
      ),
      shiny::tabPanel(
        title = "OC Simulation (SKBD)",
        shiny::sidebarLayout(
          shiny::sidebarPanel(
            shiny::numericInput("o_target", "Target toxicity (phi)", value = 0.30, min = 0.01, max = 0.99, step = 0.01),
            shiny::textInput("o_tox", "True toxicity probs (comma separated)", value = "0.05,0.12,0.30,0.45,0.60"),
            shiny::numericInput("o_ncohort", "Number of cohorts", value = 10, min = 1, step = 1),
            shiny::numericInput("o_csize", "Cohort size", value = 3, min = 1, step = 1),
            shiny::numericInput("o_ntrial", "Number of simulated trials", value = 500, min = 10, step = 10),
            shiny::numericInput("o_seed", "Seed", value = 1, min = 1, step = 1),
            shiny::actionButton("o_run", "Run simulation")
          ),
          shiny::mainPanel(
            shiny::verbatimTextOutput("o_msg"),
            shiny::tableOutput("o_summary")
          )
        )
      )
    )
  )

  server <- function(input, output, session) {
    parse_num_vec <- function(x) {
      vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
      as.numeric(vals)
    }

    boundary_res <- shiny::eventReactive(input$b_run, {
      y <- parse_num_vec(input$b_y)
      n <- parse_num_vec(input$b_n)

      validate_msg <- NULL
      if (any(is.na(y)) || any(is.na(n))) {
        validate_msg <- "Input y/n contains non-numeric values."
      } else if (length(y) != length(n)) {
        validate_msg <- "Vectors y and n must have the same length."
      } else if (input$b_d < 1 || input$b_d > length(y)) {
        validate_msg <- "Current dose index d is out of range."
      }

      if (!is.null(validate_msg)) {
        return(list(error = validate_msg))
      }

      out <- try(
        get_boundary_SKBD(
          target_prob = input$b_target,
          d = as.integer(input$b_d),
          y = y,
          n = n,
          table_type = input$b_type
        ),
        silent = TRUE
      )

      if (inherits(out, "try-error")) {
        return(list(error = as.character(out)))
      }

      list(result = out)
    })

    output$b_msg <- shiny::renderText({
      res <- boundary_res()
      if (!is.null(res$error)) {
        return(paste("Error:", res$error))
      }
      "Boundary table generated successfully."
    })

    output$b_table <- shiny::renderTable({
      res <- boundary_res()
      if (!is.null(res$error)) {
        return(NULL)
      }
      as.data.frame(res$result$boundary_tab)
    }, rownames = TRUE)

    oc_res <- shiny::eventReactive(input$o_run, {
      tox <- parse_num_vec(input$o_tox)
      if (any(is.na(tox))) {
        return(list(error = "tox_prob contains non-numeric values."))
      }

      out <- try(
        get_OC_SKBD(
          target_prob = input$o_target,
          tox_prob = tox,
          n_cohort = as.integer(input$o_ncohort),
          cohort_size = as.integer(input$o_csize),
          n_trial = as.integer(input$o_ntrial),
          seed = as.integer(input$o_seed)
        ),
        silent = TRUE
      )

      if (inherits(out, "try-error")) {
        return(list(error = as.character(out)))
      }

      list(result = out)
    })

    output$o_msg <- shiny::renderText({
      res <- oc_res()
      if (!is.null(res$error)) {
        return(paste("Error:", res$error))
      }
      "Simulation completed successfully."
    })

    output$o_summary <- shiny::renderTable({
      res <- oc_res()
      if (!is.null(res$error)) {
        return(NULL)
      }

      data.frame(
        metric = c("MTD", "PCS", "PCA", "ROD60", "ROD80"),
        value = c(
          res$result$MTD,
          round(res$result$PCS, 4),
          round(res$result$PCA, 4),
          round(res$result$ROD60, 4),
          round(res$result$ROD80, 4)
        )
      )
    }, rownames = FALSE)
  }

  invisible(shiny::runApp(shiny::shinyApp(ui = ui, server = server), launch.browser = launch.browser))
}
