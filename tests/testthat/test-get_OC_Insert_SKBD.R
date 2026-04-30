test_that("get_OC_Insert_SKBD returns insertion and summary outputs", {
  out <- get_OC_Insert_SKBD(
    target_prob = 0.3,
    tox_prob = c(0.05, 0.2, 0.35),
    dose_set = c(1, 2, 3),
    n_cohort = 2,
    cohort_size = 3,
    n_trial = 4,
    seed = 1,
    light_return = TRUE
  )

  expect_true(all(c("sel_pct_prespec", "pts_pct_prespec", "insertion", "simdata", "n_insertions") %in% names(out)))
  expect_equal(length(out$sel_pct_prespec), 3)
  expect_equal(length(out$pts_pct_prespec), 3)
  expect_true(is.list(out$insertion))
  expect_equal(length(out$n_insertions), 4)
  expect_true(all(c("Simulation", "Dose", "N", "X", "Selection") %in% names(out$simdata)))
})

test_that("get_OC_Insert_SKBD validates dose_set ordering", {
  expect_error(
    get_OC_Insert_SKBD(
      target_prob = 0.3,
      tox_prob = c(0.05, 0.2, 0.35),
      dose_set = c(1, 3, 2),
      n_cohort = 2,
      cohort_size = 3,
      n_trial = 2
    ),
    "strictly increasing"
  )
})
