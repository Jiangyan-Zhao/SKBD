server <- function(input, output, session) {
  output$app_version <- shiny::renderText({
    sprintf("Version %s | Style inspired by trialdesign.org", as.character(utils::packageVersion("SKBD")))
  })

  output$b_y_inputs <- shiny::renderUI({
    n_dose <- as.integer(input$b_n_dose)
    if (is.na(n_dose) || n_dose < 1) {
      n_dose <- 1L
    }

    shiny::tags$table(
      class = "sim-grid-table",
      shiny::tags$thead(
        shiny::tags$tr(lapply(seq_len(n_dose), function(i) shiny::tags$th(paste0("D", i))))
      ),
      shiny::tags$tbody(
        shiny::tags$tr(
          lapply(seq_len(n_dose), function(i) {
            shiny::tags$td(
              shiny::numericInput(
                inputId = paste0("b_y_", i),
                label = NULL,
                value = 0,
                min = 0,
                step = 1,
                width = "80px"
              )
            )
          })
        )
      )
    )
  })

  output$b_n_inputs <- shiny::renderUI({
    n_dose <- as.integer(input$b_n_dose)
    if (is.na(n_dose) || n_dose < 1) {
      n_dose <- 1L
    }

    shiny::tags$table(
      class = "sim-grid-table",
      shiny::tags$thead(
        shiny::tags$tr(lapply(seq_len(n_dose), function(i) shiny::tags$th(paste0("D", i))))
      ),
      shiny::tags$tbody(
        shiny::tags$tr(
          lapply(seq_len(n_dose), function(i) {
            shiny::tags$td(
              shiny::numericInput(
                inputId = paste0("b_n_", i),
                label = NULL,
                value = 0,
                min = 0,
                step = 1,
                width = "80px"
              )
            )
          })
        )
      )
    )
  })

  boundary_res <- shiny::eventReactive(input$b_run, {
    n_dose <- as.integer(input$b_n_dose)
    y <- vapply(seq_len(n_dose), function(i) {
      as.numeric(input[[paste0("b_y_", i)]])
    }, numeric(1))
    n <- vapply(seq_len(n_dose), function(i) {
      as.numeric(input[[paste0("b_n_", i)]])
    }, numeric(1))
    interval <- input$b_interval
    margin_left <- input$b_target - interval[1]
    margin_right <- interval[2] - input$b_target

    validate_msg <- NULL
    if (!length(y) || !length(n)) {
      validate_msg <- "DLTs by dose and Treated by dose cannot be empty."
    } else if (any(is.na(y)) || any(is.na(n))) {
      validate_msg <- "Input y/n contains non-numeric values."
    } else if (length(y) != n_dose || length(n) != n_dose) {
      validate_msg <- sprintf("Vectors y and n must each contain exactly %d values.", n_dose)
    } else if (any(y < 0) || any(n < 0) || any(y > n)) {
      validate_msg <- "Require 0 <= y <= n for every dose."
    } else if (input$b_d < 1 || input$b_d > n_dose) {
      validate_msg <- sprintf("Current dose index d must be between 1 and %d.", n_dose)
    } else if (input$b_start_dose < 1 || input$b_start_dose > n_dose) {
      validate_msg <- sprintf("Starting dose level must be between 1 and %d.", n_dose)
    } else if (length(interval) != 2 || any(is.na(interval)) || interval[1] <= 0 || interval[2] >= 1) {
      validate_msg <- "Target probability interval must stay inside (0, 1)."
    } else if (interval[1] >= input$b_target || interval[2] <= input$b_target) {
      validate_msg <- "Target toxicity probability must be inside the acceptable interval."
    } else if (is.na(input$b_k_left) || is.na(input$b_k_right) || input$b_k_left < 0 || input$b_k_left > 1 || input$b_k_right < 0 || input$b_k_right > 1) {
      validate_msg <- "Left-side and right-side borrowing strengths must both be within [0, 1]."
    }

    if (!is.null(validate_msg)) {
      return(list(error = validate_msg))
    }

    out <- try(
      SKBD::get_boundary_SKBD(
        target_prob = input$b_target,
        shared = isTRUE(input$b_shared),
        symmetric = isTRUE(input$b_symmetric),
        k_left = input$b_k_left,
        k_right = input$b_k_right,
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
  uploaded_scenario_mat <- shiny::reactiveVal(NULL)
  output$o_download_template <- shiny::downloadHandler(
    filename = function() {
      sprintf("scenario_template_%sdose.csv", scenario_dim$n_dose)
    },
    content = function(file) {
      n_dose <- scenario_dim$n_dose
      n_scenario <- max(3L, scenario_dim$n_scenario)
      template <- as.data.frame(
        matrix("", nrow = n_scenario, ncol = n_dose + 1),
        stringsAsFactors = FALSE
      )
      colnames(template) <- c("", paste0("Dose", seq_len(n_dose)))
      template[[1]] <- paste("Scenario", seq_len(n_scenario))
      utils::write.csv(template, file = file, row.names = FALSE, quote = FALSE, na = "")
    }
  )

  shiny::observe({
    n_dose_input <- input$b_n_dose
    n_dose <- if (is.null(n_dose_input) || is.na(n_dose_input)) 5L else as.integer(n_dose_input)
    n_dose <- max(1L, n_dose)
    scenario_dim$n_dose <- n_dose
    start_dose <- input$b_start_dose
    if (is.null(start_dose) || is.na(start_dose)) start_dose <- 1L
    current_d <- input$b_d
    if (is.null(current_d) || is.na(current_d)) current_d <- 1L
    shiny::updateNumericInput(session, "b_start_dose", min = 1, max = n_dose, value = min(as.integer(start_dose), n_dose))
    shiny::updateNumericInput(session, "b_d", min = 1, max = n_dose, value = min(as.integer(current_d), n_dose))
  })

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
    parsers <- list(
      function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
      function(path) utils::read.csv2(path, check.names = FALSE, stringsAsFactors = FALSE),
      function(path) utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    )
    dat_list <- lapply(parsers, function(parse_fun) {
      try(parse_fun(input$o_scenario_file$datapath), silent = TRUE)
    })
    dat_list <- Filter(function(x) !inherits(x, "try-error") && nrow(x) && ncol(x), dat_list)
    if (!length(dat_list)) {
      shiny::showNotification("Scenario file format is invalid.", type = "error")
      return()
    }
    dat <- dat_list[[which.max(vapply(dat_list, ncol, integer(1)))]]

    numeric_dat <- as.data.frame(lapply(dat, function(x) {
      suppressWarnings(as.numeric(trimws(as.character(x))))
    }))
    numeric_col <- vapply(numeric_dat, function(x) all(!is.na(x)), logical(1))

    if (!all(numeric_col)) {
      if (ncol(numeric_dat) > 1 && !numeric_col[1] && all(numeric_col[-1])) {
        numeric_dat <- numeric_dat[, -1, drop = FALSE]
      } else {
        shiny::showNotification(
          "Scenario file format is invalid. Use numeric toxicity values (optionally with a text Scenario column in the first column).",
          type = "error"
        )
        return()
      }
    }

    mat <- as.matrix(numeric_dat)
    if (any(is.na(mat)) || !nrow(mat) || !ncol(mat)) {
      shiny::showNotification("Scenario file format is invalid.", type = "error")
      return()
    }

    if (ncol(mat) != scenario_dim$n_dose) {
      shiny::showNotification(
        sprintf("Scenario file must contain exactly %d dose columns to align with Trial Setting.", scenario_dim$n_dose),
        type = "error"
      )
      return()
    }

    scenario_dim$n_scenario <- nrow(mat)
    uploaded_scenario_mat(mat)
  })

  output$o_scenario_grid <- shiny::renderUI({
    n_scenario <- scenario_dim$n_scenario
    n_dose <- scenario_dim$n_dose
    header <- c("Scenario", paste0("D", seq_len(n_dose)))
    uploaded_mat <- uploaded_scenario_mat()

    row_nodes <- lapply(seq_len(n_scenario), function(i) {
      row_values <- numeric(n_dose)
      for (j in seq_len(n_dose)) {
        base_value <- if (!is.null(uploaded_mat) && i <= nrow(uploaded_mat) && j <= ncol(uploaded_mat)) {
          uploaded_mat[i, j]
        } else if (i <= nrow(scenario_defaults) && j <= ncol(scenario_defaults)) {
          scenario_defaults[i, j]
        } else if (j == 1) {
          0.05
        } else {
          row_values[j - 1] + 0.06
        }

        row_values[j] <- base_value
        if (j > 1) {
          row_values[j] <- max(row_values[j], row_values[j - 1] + 0.06)
        }
        row_values[j] <- min(row_values[j], 0.99)
      }

      shiny::tags$tr(
        shiny::tags$td(paste("Scenario", i)),
        lapply(seq_len(n_dose), function(j) {
          shiny::tags$td(
            shiny::numericInput(
              inputId = paste0("o_scn_", i, "_", j),
              label = NULL,
              value = row_values[j],
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
          shared = isTRUE(input$b_shared),
          symmetric = isTRUE(input$b_symmetric),
          k_left = input$b_k_left,
          k_right = input$b_k_right,
          tox_prob = as.numeric(tox_mat[i, ]),
          start_dose = as.integer(input$b_start_dose),
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
    rows <- list()
    for (i in seq_along(res$result)) {
      one <- res$result[[i]]
      select_pct <- as.numeric(one$select_percent[seq_len(n_dose)])
      mean_treated <- colMeans(one$N)
      early_stop <- if (!is.null(one$dose_select)) {
        mean(as.numeric(one$dose_select) == -1, na.rm = TRUE) * 100
      } else if (length(one$select_percent) >= (n_dose + 1L)) {
        as.numeric(one$select_percent[n_dose + 1L])
      } else {
        0
      }

      rows[[length(rows) + 1L]] <- c(
        Metric = paste0("Scenario", i),
        as.list(rep(NA_real_, n_dose + 2L))
      )
      rows[[length(rows) + 1L]] <- c(
        Metric = "True DLT rate",
        as.list(res$tox_mat[i, ]),
        list(NA_real_, NA_real_)
      )
      rows[[length(rows) + 1L]] <- c(
        Metric = "Selection %",
        as.list(select_pct),
        list(NA_real_, early_stop)
      )
      rows[[length(rows) + 1L]] <- c(
        Metric = "# Pts treated",
        as.list(mean_treated),
        list(one$n_patient_mean, NA_real_)
      )
    }

    out_df <- do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
    colnames(out_df) <- c(
      "Metric",
      paste0("Dose ", seq_len(n_dose)),
      "Number of Patients",
      "% Early Stopping"
    )
    out_df
  })

  output$o_summary <- DT::renderDataTable({
    out_df <- oc_summary_df()

    display_df <- out_df
    numeric_cols <- setdiff(names(display_df), "Metric")
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
          list(extend = "copy", title = "SKBD_summary", exportOptions = list(modifier = list(page = "all"))),
          list(extend = "csv", filename = "SKBD_summary", exportOptions = list(modifier = list(page = "all"))),
          list(extend = "excel", filename = "SKBD_summary", exportOptions = list(modifier = list(page = "all"))),
          list(extend = "print", title = "SKBD_summary", exportOptions = list(modifier = list(page = "all")))
        ),
        pageLength = 16,
        lengthChange = FALSE,
        ordering = FALSE,
        scrollX = TRUE
      )
    )
  }, server = FALSE)
}
