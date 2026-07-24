test_that("AggregateDataPP preserves structure and transformations", {
  x <- auto_bi_fixture()
  out <- AggregateDataPP(x)
  expect_s3_class(out, "AggregateDataPP")
  expect_identical(out$cumulative.payments.triangle, x)
  expect_equal(out$incremental.payments.triangle, ChainLadder::cum2incr(x))
  expect_identical(dim(out$occurrance), dim(x))
  expect_identical(dim(out$exposure), dim(x))
})

test_that("AggregateDataPP validates inputs", {
  expect_error(AggregateDataPP(matrix(1, 2, 3)), "square numeric matrix")
  expect_error(AggregateDataPP(matrix(c(1, -1, 2, NA), 2)), "negative")
  expect_error(AggregateDataPP(matrix(c(1, 2, 0, NA), 2)), "cumulative")
  good <- matrix(c(1, 2, 2, NA), 2)
  expect_error(AggregateDataPP(good, eta = 0), "`eta`")
  expect_error(AggregateDataPP(good, entries.weights = matrix(1, 3, 3)), "entries.weights")
})

test_that("model and prediction validation is explicit", {
  pp <- AggregateDataPP(auto_bi_fixture())
  expect_error(clmplus(pp, "unknown"), "hazard.model")
  expect_error(clmplus(1, "a"), "must inherit")
})
