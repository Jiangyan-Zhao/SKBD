test_that("get_boundary_SKBD returns expected structure", {
  y <- c(0, 1, 2, 2, 0)
  n <- c(3, 6, 9, 3, 0)

  out <- get_boundary_SKBD(
    target_prob = 0.3,
    d = 3,
    y = y,
    n = n,
    table_type = "baseline"
  )

  expect_type(out, "list")
  expect_true(all(c("boundary_tab", "full_boundary_tab", "decision_table", "weight", "meta") %in% names(out)))
  expect_equal(length(out$weight), length(y))
  expect_equal(sum(out$weight), 1, tolerance = 1e-8)
})

test_that("continue table can hide past columns", {
  y <- c(0, 1, 2, 2, 0)
  n <- c(3, 6, 9, 3, 0)

  out <- get_boundary_SKBD(
    target_prob = 0.3,
    d = 3,
    y = y,
    n = n,
    table_type = "continue",
    show_past = FALSE,
    start_from = "next"
  )

  expect_true(out$meta$table_type == "continue")
  expect_true(out$meta$show_past == FALSE)
  expect_true(out$meta$start_from == "next")
  expect_true(ncol(out$full_boundary_tab) > 0)
})

test_that("continue table start_from current starts at current n", {
  y <- c(0, 1, 2, 2, 0)
  n <- c(3, 6, 9, 3, 0)

  out <- get_boundary_SKBD(
    target_prob = 0.3,
    d = 3,
    y = y,
    n = n,
    table_type = "continue",
    show_past = FALSE,
    start_from = "current"
  )

  expect_equal(out$meta$start_n, n[3])
})

test_that("shared = FALSE only uses current dose in weight", {
  y <- c(0, 1, 2, 2, 0)
  n <- c(3, 6, 9, 3, 0)

  out <- get_boundary_SKBD(
    target_prob = 0.3,
    d = 3,
    y = y,
    n = n,
    shared = FALSE
  )

  expect_equal(out$weight, c(0, 0, 1, 0, 0))
})

test_that("baseline sets first-column boundaries to NA when higher dose already treated", {
  y <- c(0, 1, 0)
  n <- c(3, 3, 3)

  out <- get_boundary_SKBD(
    target_prob = 0.3,
    d = 2,
    y = y,
    n = n,
    table_type = "baseline"
  )

  expect_true(all(is.na(out$full_boundary_tab[2:4, 1])))
})

test_that("extra_safe at lowest dose returns stop boundary and cutoff", {
  y <- c(0, 1, 0)
  n <- c(3, 3, 0)

  out <- get_boundary_SKBD(
    target_prob = 0.3,
    d = 1,
    y = y,
    n = n,
    table_type = "continue",
    extra_safe = TRUE,
    cutoff_elimin = 0.95,
    offset = 0.05
  )

  expect_true(all(c("cutoff", "stop_boundary", "target_prob") %in% names(out)))
  expect_equal(out$cutoff, 0.90)
  expect_equal(nrow(out$stop_boundary), 2)
})

test_that("constant dose_set errors", {
  y <- c(0, 1, 2)
  n <- c(3, 3, 3)

  expect_error(
    get_boundary_SKBD(
      target_prob = 0.3,
      d = 2,
      y = y,
      n = n,
      dose_set = c(1, 1, 1)
    ),
    "dose_set must have at least two distinct values"
  )
})

test_that("boundary table has expected components", {
  y <- c(0, 1, 2, 2, 0)
  n <- c(3, 6, 9, 3, 0)
  out <- get_boundary_SKBD(0.30, d = 3, y = y, n = n)
  expect_true(all(c("boundary_tab", "decision_table", "weight") %in% names(out)))
})