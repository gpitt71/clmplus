test_that("predictReserve returns the documented interface", {
  prediction <- predict(fit_age_fixture(), forecasting_horizon = 1)
  out <- predictReserve(prediction)
  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("AP", "DP", "CP", "IBNR"))
  expect_true(all(vapply(out[1:3], is.integer, logical(1))))
  expect_equal(nrow(out), 7)
})

test_that("plot methods execute", {
  fit <- fit_age_fixture()
  prediction <- predict(fit)
  expect_s3_class(suppressWarnings(plot(fit)), "ggplot")
  expect_error(suppressWarnings(plot(prediction)), NA)
  expect_error(suppressWarnings(plot(AggregateDataPP(auto_bi_fixture()))), NA)
})
