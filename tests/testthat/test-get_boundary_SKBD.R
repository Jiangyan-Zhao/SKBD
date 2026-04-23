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
