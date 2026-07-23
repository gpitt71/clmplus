test_that("vignette age fit matches the pre-refactor oracle", {
  fit <- fit_age_fixture()
  expected_effect <- c(
    NA, 0.195075405763404, -0.350329447669397, -0.748932357193024,
    -1.06899155194708, -1.3366344344573, -1.5554466616873, -1.74610887223585
  )
  names(expected_effect) <- 1:8
  expect_equal(fit$model.fit$ax, expected_effect, tolerance = 1e-14)
  names(expected_effect) <- 0:7
  expect_equal(fit$fitted_effects$fitted_development_effect, expected_effect, tolerance = 1e-14)
  expect_identical(dim(fit$fitted_development_factors), c(8L, 8L))
  expect_s3_class(fit, "clmplusmodel")
})

test_that("vignette age prediction and chain ladder remain equivalent", {
  fit <- fit_age_fixture()
  prediction <- predict(fit)
  expected_reserve <- c(
    0, 11995.0209831911, 28848.2274780733, 47752.6248189638,
    67734.6836210765, 74742.2503705983, 97541.2924208589, 102444.80844246
  )
  expected_ultimate <- c(
    63918, 74756.0209831911, 89937.2274780733, 99797.6248189638,
    107290.683621077, 96776.2503705983, 109482.292420859, 105245.80844246
  )
  expect_equal(unname(prediction$reserve), expected_reserve, tolerance = 1e-10)
  expect_equal(unname(prediction$ultimate_cost), expected_ultimate, tolerance = 1e-10)
  expect_equal(sum(prediction$reserve), 431058.9081352221, tolerance = 1e-10)
  expect_identical(dim(prediction$full_triangle), c(8L, 8L))
  expect_s3_class(prediction, "clmpluspredictions")

  mack <- ChainLadder::MackChainLadder(auto_bi_fixture())
  expect_equal(unname(prediction$full_triangle), unname(mack$FullTriangle), tolerance = 1e-10)
  expected_factors <- matrix(NA_real_, 8, 8)
  expected_factors[, -1L] <- mack$FullTriangle[, -1L] / mack$FullTriangle[, -ncol(mack$FullTriangle)]
  mask <- !is.na(prediction$development_factors_predicted)
  expect_equal(prediction$development_factors_predicted[mask], expected_factors[mask], tolerance = 1e-12)
})

test_that("all forecasting horizons allocate exact output dimensions", {
  fit <- fit_age_fixture()
  for (horizon in c(1L, 3L, 7L)) {
    prediction <- predict(fit, forecasting_horizon = horizon)
    expect_identical(dim(prediction$full_triangle), c(8L, 8L))
    expect_identical(dim(prediction$apc_output$lower_triangle_apc), c(8L, horizon))
    expect_equal(sum(!is.na(prediction$lower_triangle)), sum(7L - seq_len(horizon) + 1L))
  }
  expect_error(predict(fit, forecasting_horizon = 0), "forecasting_horizon")
  expect_error(predict(fit, forecasting_horizon = 8), "forecasting_horizon")
})
