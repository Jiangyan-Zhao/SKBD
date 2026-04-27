if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is required to run this app. Please install it first.")
}

ui <- shiny::fluidPage(
  shiny::tags$head(
    shiny::tags$style(shiny::HTML("\n      body { background-color: #efefef; font-family: 'Helvetica Neue', Arial, sans-serif; }\n      .app-title-wrap { text-align: center; margin-top: 10px; margin-bottom: 18px; }\n      .app-title { color: #7d1a8e; font-size: 42px; font-weight: 500; margin-bottom: 6px; }\n      .app-subtitle { color: #4b4b4b; font-size: 22px; margin-bottom: 4px; }\n      .app-meta { color: #8c8c8c; font-size: 16px; font-style: italic; margin-bottom: 12px; }\n      .nav-tabs > li > a { color: #1d2a7a; font-size: 24px; padding: 14px 22px; }\n      .panel-card {\n        background: #bce6e7;\n        border: 1px solid #add8da;\n        border-radius: 4px;\n        padding: 16px;\n        margin-bottom: 14px;\n      }\n      .panel-heading {\n        color: #890b97;\n        text-align: center;\n        font-size: 28px;\n        font-weight: 700;\n        margin-bottom: 12px;\n      }\n      .btn-primary {\n        background-color: #2d79be !important;\n        border-color: #236299 !important;\n      }\n      .result-box {\n        background: #ffffff;\n        border: 1px solid #dadada;\n        border-radius: 4px;\n        padding: 14px;\n        min-height: 500px;\n      }\n      .result-title { color: #ff7d5a; font-size: 30px; margin-bottom: 10px; }\n      .help-text { font-size: 18px; color: #444; margin-bottom: 12px; }\n      .small-note { color: #555; font-size: 13px; margin-top: -6px; margin-bottom: 8px; }\n    "))
  ),
  shiny::div(
    class = "app-title-wrap",
    shiny::div(class = "app-title", "SKBD: a Bayesian toxicity probability interval design explorer"),
    shiny::div(class = "app-subtitle", "Interactive app for dose-finding boundary and operating characteristics"),
    shiny::div(class = "app-meta", "Version 1.0.0 | Style inspired by trialdesign.org")
  ),
  shiny::tabsetPanel(
    shiny::tabPanel(
      title = "Boundary (SKBD)",
      shiny::fluidRow(
        shiny::column(
          width = 4,
          shiny::div(
            class = "panel-card",
            shiny::div(class = "panel-heading", "Doses"),
            shiny::textInput("b_y", "DLTs by dose (comma separated)", value = "0,1,2,2,0"),
            shiny::textInput("b_n", "Treated by dose (comma separated)", value = "3,6,9,3,0"),
            shiny::numericInput("b_d", "Current dose index d", value = 3, min = 1, step = 1),
            shiny::selectInput("b_type", "Table type", choices = c("continue", "baseline"), selected = "continue")
          ),
          shiny::div(
            class = "panel-card",
            shiny::div(class = "panel-heading", "Target Probability"),
            shiny::numericInput("b_target", "Target Toxicity Probability ϕ :", value = 0.30, min = 0.01, max = 0.99, step = 0.01),
            shiny::sliderInput(
              "b_interval",
              "Acceptable toxicity probability interval:",
              min = 0,
              max = 1,
              value = c(0.25, 0.35),
              step = 0.01
            )
          ),
          shiny::div(
            class = "panel-card",
            shiny::div(class = "panel-heading", "Sample Size"),
            shiny::numericInput("b_csize", "Cohort size", value = 3, min = 1, step = 1),
            shiny::numericInput("b_ncohort", "Number of cohorts", value = 10, min = 1, step = 1),
            shiny::numericInput("b_earlystop", "Display columns up to #patients", value = 1000, min = 1, step = 1)
          ),
          shiny::div(
            class = "panel-card",
            shiny::div(class = "panel-heading", "Overdose Control"),
            shiny::numericInput("b_cutoff", "Eliminate if Pr(p_d > phi | data) >", value = 0.95, min = 0.5, max = 0.999, step = 0.01),
            shiny::checkboxInput("b_extra_safe", "Apply extra-safe stopping rule at lowest dose", value = FALSE),
            shiny::numericInput("b_offset", "Extra-safe offset", value = 0.05, min = 0, max = 0.49, step = 0.01),
            shiny::actionButton("b_run", "Get Decision Table", class = "btn-primary btn-block")
          )
        ),
        shiny::column(
          width = 8,
          shiny::div(
            class = "result-box",
            shiny::div(class = "result-title", "Decision Table"),
            shiny::div(class = "help-text", "Dose escalation/de-escalation recommendation generated from SKBD."),
            shiny::verbatimTextOutput("b_msg"),
            shiny::tableOutput("b_table"),
            shiny::br(),
            shiny::tableOutput("b_extra_safe_table")
          )
        )
      )
    ),
    shiny::tabPanel(
      title = "OC Simulation (SKBD)",
      shiny::fluidRow(
        shiny::column(
          width = 4,
          shiny::div(
            class = "panel-card",
            shiny::div(class = "panel-heading", "Simulation"),
            shiny::numericInput("o_target", "Target toxicity (phi)", value = 0.30, min = 0.01, max = 0.99, step = 0.01),
            shiny::textInput("o_tox", "True toxicity probs (comma separated)", value = "0.05,0.12,0.30,0.45,0.60"),
            shiny::numericInput("o_ncohort", "Number of cohorts", value = 10, min = 1, step = 1),
            shiny::numericInput("o_csize", "Cohort size", value = 3, min = 1, step = 1),
            shiny::numericInput("o_ntrial", "Number of simulated trials", value = 500, min = 10, step = 10),
            shiny::numericInput("o_seed", "Seed", value = 1, min = 1, step = 1),
            shiny::actionButton("o_run", "Run Simulation", class = "btn-primary btn-block")
          )
        ),
        shiny::column(
          width = 8,
          shiny::div(
            class = "result-box",
            shiny::div(class = "result-title", "Operating Characteristics"),
            shiny::verbatimTextOutput("o_msg"),
            shiny::tableOutput("o_summary")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  parse_num_vec <- function(x) {
    parts <- strsplit(x, ",", fixed = TRUE)[[1]]
    vals <- trimws(parts)
    vals <- vals[nzchar(vals)]
    if (!length(vals)) {
      return(numeric(0))
    }
    as.numeric(vals)
  }

  boundary_res <- shiny::eventReactive(input$b_run, {
    y <- parse_num_vec(input$b_y)
    n <- parse_num_vec(input$b_n)
    n_dose <- length(y)
    interval <- input$b_interval
    margin_left <- input$b_target - interval[1]
    margin_right <- interval[2] - input$b_target

    validate_msg <- NULL
    if (!length(y) || !length(n)) {
      validate_msg <- "DLTs by dose and Treated by dose cannot be empty."
    } else if (any(is.na(y)) || any(is.na(n))) {
      validate_msg <- "Input y/n contains non-numeric values."
    } else if (length(y) != length(n)) {
      validate_msg <- "Vectors y and n must have the same length."
    } else if (any(y < 0) || any(n < 0) || any(y > n)) {
      validate_msg <- "Require 0 <= y <= n for every dose."
    } else if (input$b_d < 1 || input$b_d > n_dose) {
      validate_msg <- sprintf("Current dose index d must be between 1 and %d.", n_dose)
    } else if (length(interval) != 2 || any(is.na(interval)) || interval[1] <= 0 || interval[2] >= 1) {
      validate_msg <- "Target probability interval must stay inside (0, 1)."
    } else if (interval[1] >= input$b_target || interval[2] <= input$b_target) {
      validate_msg <- "Target toxicity probability must be inside the acceptable interval."
    }

    if (!is.null(validate_msg)) {
      return(list(error = validate_msg))
    }

    out <- try(
      SKBD::get_boundary_SKBD(
        target_prob = input$b_target,
        d = as.integer(input$b_d),
        y = y,
        n = n,
        n_cohort = as.integer(input$b_ncohort),
        cohort_size = as.integer(input$b_csize),
        n_earlystop = as.integer(input$b_earlystop),
        margin_left = margin_left,
        margin_right = margin_right,
        cutoff_elimin = input$b_cutoff,
        extra_safe = isTRUE(input$b_extra_safe),
        offset = input$b_offset,
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

  output$b_extra_safe_table <- shiny::renderTable({
    res <- boundary_res()
    if (!is.null(res$error) || is.null(res$result$stop_boundary)) {
      return(NULL)
    }
    as.data.frame(res$result$stop_boundary)
  }, rownames = TRUE)

  oc_res <- shiny::eventReactive(input$o_run, {
    tox <- parse_num_vec(input$o_tox)
    if (any(is.na(tox))) {
      return(list(error = "tox_prob contains non-numeric values."))
    }

    out <- try(
      SKBD::get_OC_SKBD(
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
      metric = c("PCS", "PCA", "Above MTD", "ROD60", "ROD80"),
      value = c(
        round(res$result$PCS, 4),
        round(res$result$PCA, 4),
        round(res$result$above_MTD, 4),
        round(res$result$ROD60, 4),
        round(res$result$ROD80, 4)
      )
    )
  }, rownames = FALSE)
}

shiny::shinyApp(ui = ui, server = server)
