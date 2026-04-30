test_that("get_OC_SKBD returns core summary fields", {
  out <- get_OC_SKBD(
    target_prob = 0.3,
    tox_prob = c(0.05, 0.15, 0.3),
    n_cohort = 2,
    cohort_size = 3,
    n_trial = 6,
    seed = 1,
    shared = FALSE,
    light_return = TRUE
  )

  expect_true(all(c("PCS", "PCA", "Y", "N", "dose_select", "select_percent") %in% names(out)))
  expect_equal(dim(out$Y), c(6, 3))
  expect_equal(dim(out$N), c(6, 3))
  expect_equal(length(out$dose_select), 6)
})

test_that("get_OC_SKBD returns paths when light_return = FALSE", {
  out <- get_OC_SKBD(
    target_prob = 0.25,
    tox_prob = c(0.05, 0.2, 0.35),
    n_cohort = 2,
    cohort_size = 3,
    n_trial = 4,
    seed = 7,
    shared = TRUE,
    light_return = FALSE
  )

  expect_true(all(c("dose_Paths", "DLT_Paths") %in% names(out)))
  expect_equal(dim(out$dose_Paths), c(4, 6))
  expect_equal(dim(out$DLT_Paths), c(4, 6))
})

test_that("get_OC_SKBD validates invalid start_dose", {
  expect_error(
    get_OC_SKBD(
      target_prob = 0.3,
      tox_prob = c(0.05, 0.15, 0.3),
      n_cohort = 2,
      cohort_size = 3,
      start_dose = 0,
      n_trial = 2
    ),
    "start_dose"
  )
})
