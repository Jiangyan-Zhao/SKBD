if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is required to run this app. Please install it first.")
}

ui <- shiny::fluidPage(
  shiny::tags$head(
    shiny::tags$style(shiny::HTML("
      body { background-color: #efefef; font-family: 'Helvetica Neue', Arial, sans-serif; }
      .app-title-wrap { text-align: center; margin-top: 10px; margin-bottom: 18px; }
      .app-title { color: #7d1a8e; font-size: 42px; font-weight: 500; margin-bottom: 6px; }
      .app-subtitle { color: #4b4b4b; font-size: 22px; margin-bottom: 4px; }
      .app-meta { color: #8c8c8c; font-size: 16px; font-style: italic; margin-bottom: 12px; }
      .nav-tabs > li > a { color: #1d2a7a; font-size: 24px; padding: 14px 22px; }
      .panel-card {
        background: #bce6e7;
        border: 1px solid #add8da;
        border-radius: 4px;
        padding: 20px;
        margin-bottom: 18px;
      }
      .panel-heading {
        color: #890b97;
        text-align: center;
        font-size: 30px;
        font-weight: 700;
        margin-bottom: 14px;
      }
      .btn-primary {
        background-color: #2d79be !important;
        border-color: #236299 !important;
      }
      .result-box {
        background: #ffffff;
        border: 1px solid #dadada;
        border-radius: 4px;
        padding: 14px;
        min-height: 500px;
      }
      .result-title { color: #ff7d5a; font-size: 30px; margin-bottom: 10px; }
      .help-text { font-size: 18px; color: #444; margin-bottom: 12px; }
    "))
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
            shiny::div(class = "panel-heading", "Trial Setting"),
            shiny::numericInput("b_target", "Target toxicity (phi)", value = 0.30, min = 0.01, max = 0.99, step = 0.01),
            shiny::textInput("b_y", "DLTs by dose (comma separated)", value = "0,1,2,2,0"),
            shiny::textInput("b_n", "Treated by dose (comma separated)", value = "3,6,9,3,0"),
            shiny::numericInput("b_d", "Current dose index d", value = 3, min = 1, step = 1),
            shiny::selectInput("b_type", "Table type", choices = c("continue", "baseline"), selected = "continue"),
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
            shiny::tableOutput("b_table")
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
    } else if (length(y) != input$n_dose) {
      validate_msg <- "Number of doses must match the length of DLTs by dose."
    } else if (length(n) != input$n_dose) {
      validate_msg <- "Number of doses must match the length of Number treated by dose."
    } else if (input$b_d < 1 || input$b_d > input$n_dose) {
      validate_msg <- "Current dose index d is out of range."
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
