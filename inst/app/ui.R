ui <- shiny::fluidPage(
  shiny::tags$head(
    shiny::tags$style(shiny::HTML("\n      body { background-color: #efefef; font-family: 'Helvetica Neue', Arial, sans-serif; }\n      .app-title-wrap { text-align: center; margin-top: 10px; margin-bottom: 18px; }\n      .app-title { color: #7d1a8e; font-size: 42px; font-weight: 500; margin-bottom: 6px; }\n      .app-subtitle { color: #4b4b4b; font-size: 22px; margin-bottom: 4px; }\n      .app-meta { color: #8c8c8c; font-size: 16px; font-style: italic; margin-bottom: 12px; }\n      .nav-tabs > li > a { color: #1d2a7a; font-size: 24px; padding: 14px 22px; }\n      .panel-card {\n        background: #bce6e7;\n        border: 1px solid #add8da;\n        border-radius: 4px;\n        padding: 16px;\n        margin-bottom: 14px;\n      }\n      .panel-heading {\n        color: #890b97;\n        text-align: center;\n        font-size: 20px;\n        font-weight: 700;\n        margin-bottom: 12px;\n      }\n      .btn-primary {\n        background-color: #2d79be !important;\n        border-color: #236299 !important;\n      }\n      .result-box {\n        background: #ffffff;\n        border: 1px solid #dadada;\n        border-radius: 4px;\n        padding: 14px;\n        min-height: 500px;\n      }\n      .result-title { color: #ff7d5a; font-size: 30px; margin-bottom: 10px; }\n      .help-text { font-size: 18px; color: #444; margin-bottom: 12px; }\n      .small-note { color: #555; font-size: 13px; margin-top: 4px; margin-bottom: 8px; line-height: 1.35; }\n      .dose-grid {\n        display: grid;\n        grid-template-columns: repeat(2, minmax(0, 1fr));\n        gap: 10px 16px;\n      }\n      .dose-alert-box {\n        background: #fff8e8;\n        border: 2px solid #f0b24d;\n        border-radius: 4px;\n        padding: 10px;\n      }\n      .dose-alert-title {\n        color: #9c4e00;\n        font-weight: 700;\n        margin-bottom: 8px;\n      }\n      .dose-alert-note {\n        color: #6b4b1e;\n        font-size: 13px;\n        margin-bottom: 8px;\n        line-height: 1.35;\n      }\n      .dose-alert-box .form-group {\n        margin-bottom: 8px;\n      }\n      .dose-grid .form-group {\n        margin-bottom: 8px;\n      }\n      .sim-grid-wrap {\n        max-height: 220px;\n        overflow-y: auto;\n        border: 1px solid #cfdfe0;\n        background: #f5f8f8;\n        padding: 8px;\n        margin-top: 8px;\n      }\n      .sim-grid-table {\n        width: 100%;\n        border-collapse: collapse;\n      }\n      .sim-grid-table th,\n      .sim-grid-table td {\n        border: 1px solid #d8d8d8;\n        padding: 4px;\n        text-align: center;\n      }\n      .sim-grid-table th {\n        background: #f0f0f0;\n        font-weight: 600;\n      }\n      .sim-grid-table .form-group {\n        margin-bottom: 0;\n      }\n      .btn-sim {\n        width: 100%;\n        margin-top: 4px;\n      }\n      .result-actions {\n        margin-top: 10px;\n      }\n      .dataTables_wrapper .dataTables_filter input {\n        margin-left: 0.5em;\n      }\n      table.dataTable tbody tr td {\n        vertical-align: middle;\n      }\n      #o_summary table.dataTable th,\n      #o_summary table.dataTable td {\n        text-align: center;\n        vertical-align: middle;\n      }\n      #o_summary table.dataTable th:first-child,\n      #o_summary table.dataTable td:first-child {\n        text-align: left;\n      }\n      #o_summary table.dataTable thead th {\n        white-space: nowrap;\n      }\n    ")),
    shiny::tags$script(shiny::HTML(
      "$(document).on('shown.bs.tab', 'a[data-toggle=\"tab\"]', function () {\
         if ($.fn.dataTable && $.fn.dataTable.isDataTable('#o_summary table')) {\
           $('#o_summary table').DataTable().columns.adjust().draw(false);\
         }\
       });"
    ))
  ),
  shiny::div(
    class = "app-title-wrap",
    shiny::div(class = "app-title", "SKBD: Shared Keyboard Design"),
    shiny::div(class = "app-subtitle", "Jiangyan Zhao · Xian Shi · Jin Xu"),
    shiny::div(class = "app-meta", shiny::textOutput("app_version", inline = TRUE))
  ),
  shiny::tabsetPanel(
    shiny::tabPanel(
      title = "Trial Setting",
      shiny::fluidRow(
        shiny::column(
          width = 4,
          shiny::div(
            class = "panel-card",
            shiny::div(class = "panel-heading", "Doses"),
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::numericInput("b_n_dose", "Number of doses:", value = 5, min = 1, step = 1)
              ),
              shiny::column(
                width = 6,
                shiny::numericInput("b_start_dose", "Starting dose level:", value = 1, min = 1, step = 1)
              )
            ),
            shiny::tags$label(shiny::strong("Dose values:")),
            shiny::uiOutput("b_dose_inputs")
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
            shiny::div(class = "panel-heading", "Kernel Setting"),
            shiny::fluidRow(
              shiny::column(
                width = 12,
                shiny::checkboxInput("b_shared", "Use shared borrowing", value = TRUE)
              )
            ),
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::sliderInput("b_k_left", "Left-side borrowing strength", value = 0.2, min = 0, max = 1, step = 0.01, width = "100%")
              ),
              shiny::column(
                width = 6,
                shiny::sliderInput("b_k_right", "Right-side borrowing strength", value = 0.8, min = 0, max = 1, step = 0.01, width = "100%")
              )
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
                shiny::numericInput("b_earlystop", "Display columns up to #patients", value = 30, min = 1, step = 1)
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
            class = "dose-alert-box",
            shiny::div(class = "dose-alert-title", "⚠ Important: Current trial data"),
            shiny::div(
              class = "dose-alert-note",
              "Please keep these four fields consistent with the ongoing cohort. Enter one value per dose level (D1, D2, ...)."
            ),
            shiny::fluidRow(
              shiny::column(
                width = 12,
                shiny::tags$label(shiny::strong("DLTs by dose:")),
                shiny::uiOutput("b_y_inputs")
              )
            ),
            shiny::fluidRow(
              shiny::column(
                width = 12,
                shiny::tags$label(shiny::strong("Treated by dose:")),
                shiny::uiOutput("b_n_inputs")
              )
            ),
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::numericInput("b_d", "Current dose index d", value = 3, min = 1, step = 1)
              ),
              shiny::column(
                width = 6,
                shiny::selectInput("b_type", "Table type", choices = c("continue", "baseline"), selected = "baseline")
              )
            ),
            shiny::tags$div(
              class = "small-note",
              "Tip: continue uses accumulated trial data; baseline shows the base keyboard rule at the current dose."
            )
          ),
          shiny::div(
            class = "result-box",
            shiny::div(class = "result-title", "Decision Table"),
            shiny::div(class = "help-text", "Dose escalation/de-escalation recommendation generated from SKBD."),
            shiny::verbatimTextOutput("b_msg"),
            shiny::tags$p(shiny::strong("Table 1: Dose escalation/de-escalation rule.")),
            DT::dataTableOutput("b_table"),
            shiny::uiOutput("b_note"),
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
              shiny::fileInput("o_scenario_file", "Scenario CSV", accept = c(".csv", ".txt")),
              shiny::downloadButton("o_download_template", "Download CSV Template", class = "btn-default btn-sim"),
              shiny::div(class = "small-note", "Please download and fill in this template to avoid upload format errors.")
            ),
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
            DT::dataTableOutput("o_summary")
          )
        )
      )
    )
  )
)
