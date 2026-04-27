if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is required to run this app. Please install it first.")
}

ui <- shiny::fluidPage(
  shiny::tags$head(
    shiny::tags$style(shiny::HTML("\n      body { background-color: #efefef; font-family: 'Helvetica Neue', Arial, sans-serif; }\n      .app-title-wrap { text-align: center; margin-top: 10px; margin-bottom: 18px; }\n      .app-title { color: #7d1a8e; font-size: 42px; font-weight: 500; margin-bottom: 6px; }\n      .app-subtitle { color: #4b4b4b; font-size: 22px; margin-bottom: 4px; }\n      .app-meta { color: #8c8c8c; font-size: 16px; font-style: italic; margin-bottom: 12px; }\n      .nav-tabs > li > a { color: #1d2a7a; font-size: 24px; padding: 14px 22px; }\n      .panel-card {\n        background: #bce6e7;\n        border: 1px solid #add8da;\n        border-radius: 4px;\n        padding: 16px;\n        margin-bottom: 14px;\n      }\n      .panel-heading {\n        color: #890b97;\n        text-align: center;\n        font-size: 20px;\n        font-weight: 700;\n        margin-bottom: 12px;\n      }\n      .btn-primary {\n        background-color: #2d79be !important;\n        border-color: #236299 !important;\n      }\n      .result-box {\n        background: #ffffff;\n        border: 1px solid #dadada;\n        border-radius: 4px;\n        padding: 14px;\n        min-height: 500px;\n      }\n      .result-title { color: #ff7d5a; font-size: 30px; margin-bottom: 10px; }\n      .help-text { font-size: 18px; color: #444; margin-bottom: 12px; }\n      .small-note { color: #555; font-size: 13px; margin-top: -6px; margin-bottom: 8px; }\n      .dose-grid {\n        display: grid;\n        grid-template-columns: repeat(2, minmax(0, 1fr));\n        gap: 10px 16px;\n      }\n      .dose-grid .form-group {\n        margin-bottom: 8px;\n      }\n      .sim-grid-wrap {\n        max-height: 220px;\n        overflow-y: auto;\n        border: 1px solid #cfdfe0;\n        background: #f5f8f8;\n        padding: 8px;\n        margin-top: 8px;\n      }\n      .sim-grid-table {\n        width: 100%;\n        border-collapse: collapse;\n      }\n      .sim-grid-table th,\n      .sim-grid-table td {\n        border: 1px solid #d8d8d8;\n        padding: 4px;\n        text-align: center;\n      }\n      .sim-grid-table th {\n        background: #f0f0f0;\n        font-weight: 600;\n      }\n      .sim-grid-table .form-group {\n        margin-bottom: 0;\n      }\n      .btn-sim {\n        width: 100%;\n        margin-top: 4px;\n      }\n    "))
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
            shiny::div(
              class = "dose-grid",
              shiny::textInput("b_y", "DLTs by dose (comma separated)", value = "0,1,2,2,0"),
              shiny::textInput("b_n", "Treated by dose (comma separated)", value = "3,6,9,3,0"),
              shiny::numericInput("b_d", "Current dose index d", value = 3, min = 1, step = 1),
              shiny::selectInput("b_type", "Table type", choices = c("continue", "baseline"), selected = "continue")
            )
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
              step = 0.01,
              width = "100%"
            )
          ),
          shiny::div(
            class = "panel-card",
            shiny::div(class = "panel-heading", "Sample Size"),
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::numericInput("b_csize", "Cohort size", value = 3, min = 1, step = 1)
              ),
              shiny::column(
                width = 4,
                shiny::numericInput("b_ncohort", "Number of cohorts", value = 10, min = 1, step = 1)
              ),
              shiny::column(
                width = 4,
                shiny::numericInput("b_earlystop", "Display columns up to #patients", value = 1000, min = 1, step = 1)
              )
            )
          ),
          shiny::div(
            class = "panel-card",
            shiny::div(class = "panel-heading", "Overdose Control"),
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::numericInput("b_cutoff", "Eliminate if Pr(p_d > phi | data) >", value = 0.95, min = 0.5, max = 0.999, step = 0.01)
              ),
              shiny::column(
                width = 6,
                shiny::numericInput("b_offset", "Extra-safe offset", value = 0.05, min = 0, max = 0.49, step = 0.01)
              )
            ),
            shiny::checkboxInput("b_extra_safe", "Apply extra-safe stopping rule at lowest dose", value = FALSE),
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
      title = "Simulation",
      shiny::fluidRow(
        shiny::column(
          width = 4,
          shiny::div(
            class = "panel-card",
            shiny::div(class = "panel-heading", "Simulation"),
            shiny::tags$label("Method to enter simulation scenarios:"),
            shiny::radioButtons(
              "o_input_method",
              label = NULL,
              choices = c("Type in" = "type", "Upload scenario file" = "upload"),
              selected = "type"
            ),
            shiny::conditionalPanel(
              condition = "input.o_input_method == 'upload'",
              shiny::fileInput("o_scenario_file", "Scenario CSV", accept = c(".csv", ".txt"))
            ),
            shiny::numericInput("o_target", "Target toxicity (phi)", value = 0.30, min = 0.01, max = 0.99, step = 0.01),
            shiny::numericInput("o_ncohort", "Number of cohorts", value = 10, min = 1, step = 1),
            shiny::numericInput("o_csize", "Cohort size", value = 3, min = 1, step = 1),
            shiny::fluidRow(
              shiny::column(8, shiny::numericInput("o_ntrial", "Number of Simulations", value = 1000, min = 10, step = 10)),
              shiny::column(4, shiny::numericInput("o_seed", "Set Seed", value = 6, min = 1, step = 1))
            ),
            shiny::fluidRow(
              shiny::column(4, shiny::actionButton("o_add_scn", "Add", class = "btn-sim")),
              shiny::column(4, shiny::actionButton("o_remove_scn", "Remove", class = "btn-sim")),
              shiny::column(4, shiny::actionButton("o_save_scn", "Save", class = "btn-sim"))
            ),
            shiny::div("For each scenario, enter true toxicity rate of each dose level:"),
            shiny::div(class = "sim-grid-wrap", shiny::uiOutput("o_scenario_grid")),
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

  scenario_dim <- shiny::reactiveValues(n_scenario = 4, n_dose = 5)
  scenario_defaults <- matrix(
    c(
      0.30, 0.47, 0.53, 0.58, 0.64,
      0.11, 0.30, 0.45, 0.56, 0.67,
      0.02, 0.07, 0.13, 0.30, 0.47,
      0.05, 0.08, 0.12, 0.15, 0.30
    ),
    nrow = 4,
    byrow = TRUE
  )

  shiny::observeEvent(input$o_add_scn, {
    scenario_dim$n_scenario <- scenario_dim$n_scenario + 1
  })

  shiny::observeEvent(input$o_remove_scn, {
    scenario_dim$n_scenario <- max(1L, scenario_dim$n_scenario - 1L)
  })

  shiny::observeEvent(input$o_save_scn, {
    shiny::showNotification("Scenarios saved in current session.", type = "message")
  })

  shiny::observeEvent(input$o_scenario_file, {
    shiny::req(input$o_scenario_file)
    dat <- try(utils::read.csv(input$o_scenario_file$datapath, check.names = FALSE), silent = TRUE)
    if (inherits(dat, "try-error") || !nrow(dat) || !ncol(dat)) {
      return()
    }
    mat <- as.matrix(dat)
    suppressWarnings(storage.mode(mat) <- "numeric")
    if (any(is.na(mat))) {
      return()
    }
    scenario_dim$n_scenario <- nrow(mat)
    scenario_dim$n_dose <- ncol(mat)
    session$onFlushed(function() {
      for (i in seq_len(nrow(mat))) {
        for (j in seq_len(ncol(mat))) {
          shiny::updateNumericInput(session, paste0("o_scn_", i, "_", j), value = mat[i, j])
        }
      }
    }, once = TRUE)
  })

  output$o_scenario_grid <- shiny::renderUI({
    n_scenario <- scenario_dim$n_scenario
    n_dose <- scenario_dim$n_dose
    header <- c("Scenario", paste0("D", seq_len(n_dose)))

    row_nodes <- lapply(seq_len(n_scenario), function(i) {
      shiny::tags$tr(
        shiny::tags$td(paste("Scenario", i)),
        lapply(seq_len(n_dose), function(j) {
          default_value <- if (i <= nrow(scenario_defaults) && j <= ncol(scenario_defaults)) scenario_defaults[i, j] else 0.3
          shiny::tags$td(
            shiny::numericInput(
              inputId = paste0("o_scn_", i, "_", j),
              label = NULL,
              value = default_value,
              min = 0,
              max = 1,
              step = 0.01,
              width = "70px"
            )
          )
        })
      )
    })

    shiny::tags$table(
      class = "sim-grid-table",
      shiny::tags$thead(shiny::tags$tr(lapply(header, shiny::tags$th))),
      shiny::tags$tbody(row_nodes)
    )
  })

  get_scenario_matrix <- shiny::reactive({
    n_scenario <- scenario_dim$n_scenario
    n_dose <- scenario_dim$n_dose
    mat <- matrix(NA_real_, nrow = n_scenario, ncol = n_dose)
    for (i in seq_len(n_scenario)) {
      for (j in seq_len(n_dose)) {
        mat[i, j] <- input[[paste0("o_scn_", i, "_", j)]]
      }
    }
    mat
  })

  oc_res <- shiny::eventReactive(input$o_run, {
    tox_mat <- get_scenario_matrix()
    if (!nrow(tox_mat) || !ncol(tox_mat) || any(is.na(tox_mat))) {
      return(list(error = "Please provide numeric toxicity probabilities for all scenario cells."))
    }
    if (any(tox_mat <= 0 | tox_mat >= 1)) {
      return(list(error = "All toxicity probabilities must be within (0, 1)."))
    }

    scenario_outputs <- vector("list", nrow(tox_mat))
    for (i in seq_len(nrow(tox_mat))) {
      out <- try(
        SKBD::get_OC_SKBD(
          target_prob = input$o_target,
          tox_prob = as.numeric(tox_mat[i, ]),
          n_cohort = as.integer(input$o_ncohort),
          cohort_size = as.integer(input$o_csize),
          n_trial = as.integer(input$o_ntrial),
          seed = as.integer(input$o_seed + i - 1L)
        ),
        silent = TRUE
      )
      if (inherits(out, "try-error")) {
        return(list(error = sprintf("Scenario %d failed: %s", i, as.character(out))))
      }
      scenario_outputs[[i]] <- out
    }

    list(result = scenario_outputs, tox_mat = tox_mat)
  })

  output$o_msg <- shiny::renderText({
    res <- oc_res()
    if (!is.null(res$error)) {
      return(paste("Error:", res$error))
    }
    sprintf("Simulation completed successfully for %d scenarios.", nrow(res$tox_mat))
  })

  output$o_summary <- shiny::renderTable({
    res <- oc_res()
    if (!is.null(res$error)) {
      return(NULL)
    }

    n_dose <- ncol(res$tox_mat)
    rows <- list()
    for (i in seq_along(res$result)) {
      one <- res$result[[i]]
      select_pct <- as.numeric(one$select_percent[seq_len(n_dose)])
      mean_treated <- colMeans(one$N)
      early_stop <- if ("-1" %in% names(one$select_percent)) as.numeric(one$select_percent[["-1"]]) else 0

      rows[[length(rows) + 1L]] <- c(Metric = paste0("Scenario", i), as.list(rep(NA_real_, n_dose + 2L)))
      rows[[length(rows) + 1L]] <- c(Metric = "True DLT rate", as.list(res$tox_mat[i, ]), list(NA_real_, NA_real_))
      rows[[length(rows) + 1L]] <- c(Metric = "Selection %", as.list(select_pct), list(NA_real_, early_stop))
      rows[[length(rows) + 1L]] <- c(Metric = "# Pts treated", as.list(mean_treated), list(one$n_patient_mean, NA_real_))
    }

    out_df <- do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
    colnames(out_df) <- c("Metric", paste0("Dose ", seq_len(n_dose)), "Number of Patients", "% Early Stopping")

    numeric_cols <- setdiff(names(out_df), "Metric")
    for (nm in numeric_cols) {
      out_df[[nm]] <- suppressWarnings(as.numeric(out_df[[nm]]))
      out_df[[nm]] <- ifelse(is.na(out_df[[nm]]), "", sprintf("%.2f", out_df[[nm]]))
    }
    out_df
  }, rownames = FALSE)
}

shiny::shinyApp(ui = ui, server = server)
