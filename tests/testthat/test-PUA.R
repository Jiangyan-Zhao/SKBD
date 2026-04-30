test_that("PUA returns monotone scenarios with expected dimensions", {
  scen <- PUA(
    dose_set = 1:4,
    target_prob = 0.3,
    n_scenarios = 8,
    seed = 123
  )

  expect_equal(dim(scen), c(8, 4))
  expect_true(all(apply(scen, 1, function(x) all(diff(x) >= 0))))
  expect_true(all(scen >= 0 & scen <= 1))
})

test_that("PUA is reproducible under the same seed", {
  scen1 <- PUA(dose_set = 1:5, target_prob = 0.25, n_scenarios = 6, seed = 2024)
  scen2 <- PUA(dose_set = 1:5, target_prob = 0.25, n_scenarios = 6, seed = 2024)

  expect_equal(scen1, scen2)
})

test_that("PUA target-dose position is within target interval for every scenario", {
  target_prob <- 0.3
  margin_left <- 0.05
  margin_right <- 0.05
  scen <- PUA(
    dose_set = 1:5,
    target_prob = target_prob,
    n_scenarios = 10,
    margin_left = margin_left,
    margin_right = margin_right,
    seed = 10
  )

  in_target <- apply(scen, 1, function(x) {
    d <- which.min(abs(x - target_prob))
    x[d] >= (target_prob - margin_left) && x[d] <= (target_prob + margin_right)
  })

  expect_true(all(in_target))
})
