info_popover <- function(label, title, content, aria_label) {
  shiny::div(
    class = "panel-heading",
    style = "position:relative;",
    label,
    shiny::tags$button(
      type = "button",
      class = "btn btn-link",
      style = "position:absolute;right:0;top:50%;transform:translateY(-50%);padding:0;line-height:1;border:none;text-decoration:none;color:#f39c12;font-size:18px;font-weight:700;",
      `aria-label` = aria_label,
      `data-toggle` = "popover",
      `data-container` = "body",
      title = title,
      `data-content` = content,
      shiny::HTML("&#9432;")
    )
  )
}

ui <- shiny::fluidPage(
  shiny::tags$head(
    shiny::tags$style(shiny::HTML("\n      body { background-color: #efefef; font-family: 'Helvetica Neue', Arial, sans-serif; }\n      .app-title-wrap { text-align: center; margin-top: 10px; margin-bottom: 18px; }\n      .app-title { color: #7d1a8e; font-size: 42px; font-weight: 500; margin-bottom: 6px; }\n      .app-subtitle { color: #4b4b4b; font-size: 22px; margin-bottom: 4px; }\n      .app-meta { color: #8c8c8c; font-size: 16px; font-style: italic; margin-bottom: 12px; }\n      .nav-tabs > li > a { color: #1d2a7a; font-size: 24px; padding: 14px 22px; }\n      .panel-card {\n        background: #bce6e7;\n        border: 1px solid #add8da;\n        border-radius: 4px;\n        padding: 16px;\n        margin-bottom: 14px;\n      }\n      .panel-heading {\n        color: #890b97;\n        text-align: center;\n        font-size: 20px;\n        font-weight: 700;\n        margin-bottom: 12px;\n      }\n      .btn-primary {\n        background-color: #2d79be !important;\n        border-color: #236299 !important;\n      }\n      .result-box {\n        background: #ffffff;\n        border: 1px solid #dadada;\n        border-radius: 4px;\n        padding: 14px;\n        min-height: 500px;\n      }\n      .result-title { color: #ff7d5a; font-size: 30px; margin-bottom: 10px; }\n      .help-text { font-size: 18px; color: #444; margin-bottom: 12px; }\n      .small-note { color: #555; font-size: 13px; margin-top: 4px; margin-bottom: 8px; line-height: 1.35; }\n      .dose-grid {\n        display: grid;\n        grid-template-columns: repeat(2, minmax(0, 1fr));\n        gap: 10px 16px;\n      }\n      .dose-alert-box {\n        background: #fff8e8;\n        border: 2px solid #f0b24d;\n        border-radius: 4px;\n        padding: 10px;\n      }\n      .dose-alert-title {\n        color: #9c4e00;\n        font-weight: 700;\n        margin-bottom: 8px;\n      }\n      .dose-alert-note {\n        color: #6b4b1e;\n        font-size: 13px;\n        margin-bottom: 8px;\n        line-height: 1.35;\n      }\n      .dose-alert-box .form-group {\n        margin-bottom: 8px;\n      }\n      .dose-grid .form-group {\n        margin-bottom: 8px;\n      }\n      .sim-grid-wrap {\n        max-height: 220px;\n        overflow-y: auto;\n        border: 1px solid #cfdfe0;\n        background: #f5f8f8;\n        padding: 8px;\n        margin-top: 8px;\n      }\n      .input-table-scroll {\n        width: 100%;\n        max-width: 100%;\n        overflow-x: auto;\n        overflow-y: visible;\n        box-sizing: border-box;\n        padding-bottom: 4px;\n      }\n      .input-table-scroll .sim-grid-table {\n        min-width: max-content;\n      }\n      .input-table-scroll .form-group {\n        margin-bottom: 0;\n      }\n      .sim-grid-table {\n        width: 100%;\n        border-collapse: collapse;\n      }\n      .sim-grid-table th,\n      .sim-grid-table td {\n        border: 1px solid #d8d8d8;\n        padding: 4px;\n        text-align: center;\n      }\n      .sim-grid-table th {\n        background: #f0f0f0;\n        font-weight: 600;\n      }\n      .sim-grid-table .form-group {\n        margin-bottom: 0;\n      }\n      .btn-sim {\n        width: 100%;\n        margin-top: 4px;\n      }\n      .result-actions {\n        margin-top: 10px;\n      }\n      .dataTables_wrapper .dataTables_filter input {\n        margin-left: 0.5em;\n      }\n      table.dataTable tbody tr td {\n        vertical-align: middle;\n      }\n      #o_summary table.dataTable th,\n      #o_summary table.dataTable td {\n        text-align: center;\n        vertical-align: middle;\n      }\n      #o_summary table.dataTable th:first-child,\n      #o_summary table.dataTable td:first-child {\n        text-align: left;\n      }\n      #o_summary table.dataTable thead th {\n        white-space: nowrap;\n      }\n    ")),
    shiny::tags$script(shiny::HTML(
      "$(document).on('shown.bs.tab', 'a[data-toggle=\"tab\"]', function () {\
         if ($.fn.dataTable && $.fn.dataTable.isDataTable('#o_summary table')) {\
           $('#o_summary table').DataTable().columns.adjust().draw(false);\
         }\
       });\
       $(function () {\
         $('[data-toggle=\"popover\"]').popover({trigger: 'focus', html: true, placement: 'left'});\
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
            info_popover(
              label = "Doses",
              title = "Dose settings",
              content = "Specify the number of dose levels, the starting dose level, and the numerical dose values.<br/>Dose values must be strictly increasing. They are used to define distances between doses for kernel-based borrowing under SKBD.",
              aria_label = "Dose settings help"
            ),
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::numericInput("b_n_dose", "Number of doses:", value = 5, min = 1, step = 1, width = "100%")
              ),
              shiny::column(
                width = 6,
                shiny::numericInput("b_start_dose", "Starting dose level:", value = 1, min = 1, step = 1, width = "100%")
              )
            ),
            shiny::uiOutput("b_dose_values")
          ),
          shiny::div(
            class = "panel-card",
            info_popover(
              label = "Target Probability",
              title = "Target probability",
              content = "Specify the target toxicity probability and the acceptable toxicity interval.<br/>The target probability <strong>&phi;</strong> is the desired DLT rate for the MTD.<br/>The acceptable interval defines the target key used by the Keyboard/SKBD decision rule.<br/>For example, if &phi; = 0.30 and the interval is [0.25, 0.35], doses with toxicity probabilities in this interval are considered acceptable.",
              aria_label = "Target probability help"
            ),
            shiny::numericInput("b_target", "Target Toxicity Probability ϕ :", value = 0.30, min = 0.01, max = 0.99, step = 0.01, width = "100%"),
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
            info_popover(
              label = "Sample Size",
              title = "Sample size",
              content = "Specify the cohort size, number of cohorts, and early-stopping threshold.<br/>The maximum sample size is approximately cohort size multiplied by number of cohorts.<br/>The early-stopping rule stops the trial if the number of patients assigned to a single dose reaches the specified threshold and the current decision is to stay.",
              aria_label = "Sample size help"
            ),
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::numericInput("b_csize", "Cohort size", value = 3, min = 1, step = 1, width = "100%")
              ),
              shiny::column(
                width = 6,
                shiny::numericInput("b_ncohort", "Number of cohorts", value = 10, min = 1, step = 1, width = "100%")
              )
            ),
            shiny::fluidRow(
              shiny::column(
                width = 12,
                shiny::numericInput("b_earlystop", "Stop trial if the number of patients assigned to single dose reaches n and the decision is to stay, where n =", value = 30, min = 1, step = 1, width = "100%")
              )
            )
          ),
          shiny::div(
            class = "panel-card",
            info_popover(
              label = "Overdose Control",
              title = "Overdose control",
              content = "Specify safety rules for eliminating overly toxic doses.<br/>A dose is eliminated if the posterior probability of toxicity exceeding the target is larger than the elimination cutoff.<br/>The extra-safe stopping rule provides additional protection when the lowest dose appears too toxic.<br/>The offset controls how conservative this extra-safe rule is.",
              aria_label = "Overdose control help"
            ),
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::numericInput("b_cutoff", "Eliminate if Pr(p_d > phi | data) >", value = 0.95, min = 0.5, max = 0.999, step = 0.01, width = "100%")
              ),
              shiny::column(
                width = 6,
                shiny::numericInput("b_offset", "Extra-safe offset", value = 0.05, min = 0, max = 0.49, step = 0.01, width = "100%")
              )
            ),
            shiny::checkboxInput("b_extra_safe", "Apply extra-safe stopping rule at lowest dose", value = FALSE),
            shiny::actionButton("b_run", "Get Decision Table", class = "btn-primary btn-block")
          )
        ),
        shiny::column(
          width = 8,
          shiny::div(
            class = "panel-card",
            info_popover(
              label = "Kernel Setting",
              title = "Shared borrowing",
              content = "Checked: SKBD uses kernel-based borrowing across dose levels.<br/>Unchecked: the design reduces to the ordinary Keyboard design, and the left- and right-side borrowing strengths are ignored.",
              aria_label = "Shared borrowing help"
            ),
            shiny::fluidRow(
              shiny::column(
                width = 12,
                shiny::div(
                  style = "background:#fff4d6;border:2px solid #f0b24d;border-left:6px solid #d98100;border-radius:6px;padding:10px 12px 6px 12px;margin-bottom:12px;box-shadow:0 1px 3px rgba(0,0,0,0.08);",
                                    shiny::checkboxInput("b_shared", shiny::tags$strong("Use shared borrowing"), value = TRUE)
                )
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
            class = "dose-alert-box",
            shiny::div(
              class = "panel-heading dose-alert-title",
              style = "position:relative;",
              "Current Trial Data",
              shiny::tags$button(
                type = "button",
                class = "btn btn-link",
                style = "position:absolute;right:0;top:50%;transform:translateY(-50%);padding:0;line-height:1;border:none;text-decoration:none;color:#f39c12;font-size:18px;font-weight:700;",
                `aria-label` = "Current trial data help",
                `data-toggle` = "popover",
                `data-container` = "body",
                title = "Current trial data",
                `data-content` = "Enter the interim data available at the current decision point.<br/>
<strong>DLTs by dose</strong> and <strong>Treated by dose</strong> are cumulative counts at each dose level.<br/>
<strong>Current dose index d</strong> is the dose currently under evaluation for the next cohort.<br/>
Use <strong>continue</strong> to generate a decision using the accumulated trial data; use <strong>baseline</strong> to display the base Keyboard rule at the current dose.<br/>
Require 0 &le; DLTs &le; treated for every dose.",
                shiny::HTML("&#9432;")
              )
            ),
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
            info_popover(
              label = "Simulation",
              title = "Simulation scenarios",
              content = "Specify true dose-toxicity scenarios for evaluating operating characteristics.<br/>Each row represents one simulation scenario, and each cell gives the true toxicity probability at one dose level.<br/>The simulation uses the trial settings defined in the Trial Setting tab, including the target probability, sample size, borrowing mode, and overdose-control settings.",
              aria_label = "Simulation scenarios help"
            ),
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
              shiny::column(4, shiny::numericInput("o_seed", "Set Seed", value = 6, min = 1, step = 1, width = "100%"))
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
