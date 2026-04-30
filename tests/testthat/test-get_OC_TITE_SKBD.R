test_that("get_OC_TITE_SKBD returns core summary outputs", {
  out <- get_OC_TITE_SKBD(
    target_prob = 0.2,
    tox_prob = c(0.05, 0.15, 0.25),
    n_cohort = 2,
    cohort_size = 3,
    n_trial = 4,
    seed = 2,
    tau = 1,
    accrual = 1,
    shared = FALSE,
    light_return = TRUE,
    dist_DLT = "uniform",
    dist_enter = "uniform"
  )

  expect_true(all(c("PCS", "PCA", "Y", "N", "dose_select", "duration_mean") %in% names(out)))
  expect_equal(dim(out$Y), c(4, 3))
  expect_equal(dim(out$N), c(4, 3))
})

test_that("get_OC_TITE_SKBD returns paths when light_return = FALSE", {
  out <- get_OC_TITE_SKBD(
    target_prob = 0.25,
    tox_prob = c(0.1, 0.2, 0.35),
    n_cohort = 2,
    cohort_size = 3,
    n_trial = 3,
    seed = 3,
    tau = 1,
    accrual = 1,
    light_return = FALSE,
    dist_DLT = "uniform",
    dist_enter = "exp"
  )

  expect_true(all(c("dose_Paths", "DLT_Paths") %in% names(out)))
  expect_equal(dim(out$dose_Paths), c(3, 6))
  expect_equal(dim(out$DLT_Paths), c(3, 6))
})

test_that("get_OC_TITE_SKBD validates prior_p length", {
  expect_error(
    get_OC_TITE_SKBD(
      target_prob = 0.2,
      tox_prob = c(0.05, 0.15, 0.25),
      n_cohort = 2,
      cohort_size = 3,
      n_trial = 2,
      tau = 1,
      accrual = 1,
      prior_p = c(0.2, 0.3)
    ),
    "prior_p"
  )
})
