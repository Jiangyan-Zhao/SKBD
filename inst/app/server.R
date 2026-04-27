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

  output$app_version <- shiny::renderText({
    sprintf("Version %s | Style inspired by trialdesign.org", as.character(utils::packageVersion("SKBD")))
  })

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

  output$b_table <- DT::renderDataTable({
    res <- boundary_res()
    if (!is.null(res$error)) {
      return(NULL)
    }

    boundary_df <- as.data.frame(res$result$boundary_tab, stringsAsFactors = FALSE)
    boundary_df <- cbind(Decision = rownames(boundary_df), boundary_df, row.names = NULL)
    colnames(boundary_df) <- c(" ", as.character(seq_len(ncol(boundary_df) - 1)))

    for (nm in colnames(boundary_df)[-1]) {
      vals <- suppressWarnings(as.numeric(boundary_df[[nm]]))
      boundary_df[[nm]] <- ifelse(is.na(vals), "NA", as.character(as.integer(round(vals))))
    }

    DT::datatable(
      boundary_df,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel", "print"),
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        ordering = FALSE,
        scrollX = TRUE
      )
    )
  })

  output$b_note <- shiny::renderUI({
    res <- boundary_res()
    if (!is.null(res$error)) {
      return(NULL)
    }
    shiny::tags$p(
      shiny::HTML("<strong>Note.</strong> # of DLT is the number of patients with at least 1 DLT. When none of the actions (i.e., escalate, de-escalate, or eliminate) is triggered, stay at the current dose for treating the next cohort of patients. \"NA\" means that a boundary is not available under the current setting."),
      style = "margin-top: 8px; font-weight: 600;"
    )
  })

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
      0.08, 0.30, 0.38, 0.42, 0.52,
      0.04, 0.07, 0.30, 0.35, 0.42,
      0.06, 0.07, 0.12, 0.30, 0.40,
      0.01, 0.02, 0.04, 0.06, 0.30
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
          target_prob = input$b_target,
          tox_prob = as.numeric(tox_mat[i, ]),
          n_cohort = as.integer(input$b_ncohort),
          cohort_size = as.integer(input$b_csize),
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

  oc_summary_df <- shiny::reactive({
    res <- oc_res()
    shiny::req(is.null(res$error))

    n_dose <- ncol(res$tox_mat)
    method_label <- if (identical(input$o_input_method, "upload")) {
      file_name <- input$o_scenario_file$name
      if (!is.null(file_name) && nzchar(file_name)) {
        sprintf("Upload scenario file (%s)", file_name)
      } else {
        "Upload scenario file"
      }
    } else {
      "Type in"
    }
    rows <- list()
    for (i in seq_along(res$result)) {
      one <- res$result[[i]]
      select_pct <- as.numeric(one$select_percent[seq_len(n_dose)])
      mean_treated <- colMeans(one$N)
      early_stop <- if ("-1" %in% names(one$select_percent)) as.numeric(one$select_percent[["-1"]]) else 0

      rows[[length(rows) + 1L]] <- c(
        Metric = paste0("Scenario", i),
        `Scenario Input Method` = method_label,
        as.list(rep(NA_real_, n_dose + 2L))
      )
      rows[[length(rows) + 1L]] <- c(
        Metric = "True DLT rate",
        `Scenario Input Method` = "",
        as.list(res$tox_mat[i, ]),
        list(NA_real_, NA_real_)
      )
      rows[[length(rows) + 1L]] <- c(
        Metric = "Selection %",
        `Scenario Input Method` = "",
        as.list(select_pct),
        list(NA_real_, early_stop)
      )
      rows[[length(rows) + 1L]] <- c(
        Metric = "# Pts treated",
        `Scenario Input Method` = "",
        as.list(mean_treated),
        list(one$n_patient_mean, NA_real_)
      )
    }

    out_df <- do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
    colnames(out_df) <- c(
      "Metric",
      "Scenario Input Method",
      paste0("Dose ", seq_len(n_dose)),
      "Number of Patients",
      "% Early Stopping"
    )
    out_df
  })

  output$o_summary <- DT::renderDataTable({
    out_df <- oc_summary_df()

    display_df <- out_df
    numeric_cols <- setdiff(names(display_df), c("Metric", "Scenario Input Method"))
    for (nm in numeric_cols) {
      display_df[[nm]] <- suppressWarnings(as.numeric(display_df[[nm]]))
      display_df[[nm]] <- ifelse(is.na(display_df[[nm]]), "", sprintf("%.2f", display_df[[nm]]))
    }

    DT::datatable(
      display_df,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "copy", title = "SKBD_summary_", exportOptions = list(modifier = list(page = "all"))),
          list(extend = "csv", filename = "SKBD_summary_", exportOptions = list(modifier = list(page = "all"))),
          list(extend = "excel", filename = "SKBD_summary_", exportOptions = list(modifier = list(page = "all"))),
          list(extend = "print", title = "SKBD_summary_", exportOptions = list(modifier = list(page = "all")))
        ),
        pageLength = 16,
        lengthChange = FALSE,
        ordering = FALSE,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
}
