test_that("plot methods execute", {
  fit <- fit_age_fixture()
  prediction <- predict(fit)
  expect_no_warning(p <- plot(fit))
  expect_s3_class(p, "ggplot")
  expect_no_warning(prediction_plot <- plot(prediction))
  expect_s3_class(prediction_plot, "gtable")
  expect_no_warning(data_plot <- plot(AggregateDataPP(auto_bi_fixture())))
  expect_s3_class(data_plot, "gtable")
})

test_that("plotting branches for every supported model are warning-free", {
  data <- AggregateDataPP(auto_bi_fixture())
  for (model_name in c("a", "ac", "ap", "apc")) {
    fit <- suppressWarnings(clmplus(data, model_name))
    expect_no_warning(expect_s3_class(plot(fit), "ggplot"))
    prediction <- predict(fit, forecasting_horizon = 1)
    expect_no_warning(expect_s3_class(plot(prediction), "gtable"))
  }
})
