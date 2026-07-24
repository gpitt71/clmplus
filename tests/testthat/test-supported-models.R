test_that("every documented hazard model fits and predicts", {
  pp <- AggregateDataPP(auto_bi_fixture())
  for (model_name in c("a", "ac", "ap", "apc")) {
    set.seed(42)
    fit <- suppressWarnings(clmplus(pp, model_name))
    prediction <- suppressWarnings(predict(fit, forecasting_horizon = 1))
    expect_s3_class(fit, "clmplusmodel")
    expect_s3_class(prediction, "clmpluspredictions")
    expect_identical(dim(fit$fitted_development_factors), c(8L, 8L))
    expect_identical(dim(prediction$full_triangle), c(8L, 8L))
    expect_identical(dim(prediction$apc_output$lower_triangle_apc), c(8L, 1L))
  }
})
