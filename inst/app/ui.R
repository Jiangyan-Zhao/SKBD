ui <- shiny::fluidPage(
  shiny::tags$head(
    shiny::tags$style(shiny::HTML("\n      body { background-color: #efefef; font-family: 'Helvetica Neue', Arial, sans-serif; }\n      .app-title-wrap { text-align: center; margin-top: 10px; margin-bottom: 18px; }\n      .app-title { color: #7d1a8e; font-size: 42px; font-weight: 500; margin-bottom: 6px; }\n      .app-subtitle { color: #4b4b4b; font-size: 22px; margin-bottom: 4px; }\n      .app-meta { color: #8c8c8c; font-size: 16px; font-style: italic; margin-bottom: 12px; }\n      .nav-tabs > li > a { color: #1d2a7a; font-size: 24px; padding: 14px 22px; }\n      .panel-card {\n        background: #bce6e7;\n        border: 1px solid #add8da;\n        border-radius: 4px;\n        padding: 16px;\n        margin-bottom: 14px;\n      }\n      .panel-heading {\n        color: #890b97;\n        text-align: center;\n        font-size: 20px;\n        font-weight: 700;\n        margin-bottom: 12px;\n      }\n      .btn-primary {\n        background-color: #2d79be !important;\n        border-color: #236299 !important;\n      }\n      .result-box {\n        background: #ffffff;\n        border: 1px solid #dadada;\n        border-radius: 4px;\n        padding: 14px;\n        min-height: 500px;\n      }\n      .result-title { color: #ff7d5a; font-size: 30px; margin-bottom: 10px; }\n      .help-text { font-size: 18px; color: #444; margin-bottom: 12px; }\n      .small-note { color: #555; font-size: 13px; margin-top: -6px; margin-bottom: 8px; }\n      .dose-grid {\n        display: grid;\n        grid-template-columns: repeat(2, minmax(0, 1fr));\n        gap: 10px 16px;\n      }\n      .dose-grid .form-group {\n        margin-bottom: 8px;\n      }\n      .sim-grid-wrap {\n        max-height: 220px;\n        overflow-y: auto;\n        border: 1px solid #cfdfe0;\n        background: #f5f8f8;\n        padding: 8px;\n        margin-top: 8px;\n      }\n      .sim-grid-table {\n        width: 100%;\n        border-collapse: collapse;\n      }\n      .sim-grid-table th,\n      .sim-grid-table td {\n        border: 1px solid #d8d8d8;\n        padding: 4px;\n        text-align: center;\n      }\n      .sim-grid-table th {\n        background: #f0f0f0;\n        font-weight: 600;\n      }\n      .sim-grid-table .form-group {\n        margin-bottom: 0;\n      }\n      .btn-sim {\n        width: 100%;\n        margin-top: 4px;\n      }\n      .result-actions {\n        margin-top: 10px;\n      }\n      .dataTables_wrapper .dataTables_filter input {\n        margin-left: 0.5em;\n      }\n      table.dataTable tbody tr td {\n        vertical-align: middle;\n      }\n    "))
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
            DT::dataTableOutput("o_summary"),
            shiny::div(
              class = "result-actions",
              shiny::downloadButton("o_download_csv", "Download Summary (CSV)", class = "btn-primary")
            )
          )
        )
      )
    )
  )
)
