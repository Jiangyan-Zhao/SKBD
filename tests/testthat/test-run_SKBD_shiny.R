test_that("run_SKBD_shiny fails gracefully without required runtime dependencies", {
  expect_error(
    run_SKBD_shiny(launch.browser = FALSE),
    "Cannot find the Shiny app directory|Package 'shiny' is required|Package 'DT' is required"
  )
})

test_that("run_SKBD_shiny accepts additional args but still fails predictably in source context", {
  expect_error(
    run_SKBD_shiny(launch.browser = FALSE, port = 8181),
    "Cannot find the Shiny app directory|Package 'shiny' is required|Package 'DT' is required"
  )
})
