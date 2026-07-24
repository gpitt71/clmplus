test_that("StMoMo fitting is quiet by default and verbose on request", {
  data <- AggregateDataPP(auto_bi_fixture())

  expect_output(
    quiet_fit <- suppressWarnings(clmplus(data, hazard.model = "a"))
    , NA
  )
  expect_output(
    verbose_fit <- suppressWarnings(
      clmplus(data, hazard.model = "a", verbose = TRUE)
    ),
    "Start fitting with gnm"
  )

  expect_equal(verbose_fit$model.fit$ax, quiet_fit$model.fit$ax, tolerance = 0)
  expect_equal(
    verbose_fit$fitted_development_factors,
    quiet_fit$fitted_development_factors,
    tolerance = 0
  )
})

test_that("only Poisson non-integer response warnings are suppressed", {
  data <- AggregateDataPP(auto_bi_fixture())
  observed <- character()
  withCallingHandlers(
    clmplus(data, hazard.model = "a"),
    warning = function(w) {
      observed <<- c(observed, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("^non-integer x\\s*=", observed)))
  expect_error(clmplus(data, hazard.model = "a", verbose = NA), "verbose")
})

test_that("warnings unrelated to non-integer Poisson responses remain visible", {
  data <- AggregateDataPP(auto_bi_fixture())
  baseline <- suppressWarnings(clmplus(data, hazard.model = "a"))

  testthat::local_mocked_bindings(
    fit = function(...) {
      warning("simulated convergence warning", call. = FALSE)
      baseline$model.fit
    },
    .package = "StMoMo"
  )

  expect_warning(
    clmplus(data, hazard.model = "a"),
    "simulated convergence warning"
  )
})
