test_that("triangle and calendar transformations preserve valid cells", {
  x <- matrix(c(1, 2, 3, 4, 5, NA, 6, NA, NA), 3, 3)
  calendar <- clmplus:::triangle_to_calendar(x)
  expect_identical(clmplus:::calendar_to_triangle(calendar), x)
  expect_identical(dim(calendar), c(3L, 3L))
})

test_that("full calendar transformation maps every cell", {
  x <- matrix(seq_len(9), 3, 3)
  out <- clmplus:::triangle_to_full_calendar(x)
  expect_identical(dim(out), c(3L, 6L))
  index <- which(!is.na(out), arr.ind = TRUE)
  expect_length(index[, 1L], 9L)
})

test_that("forecast placement creates only the requested lower diagonals", {
  forecast <- matrix(seq_len(12), 4, 3)
  lower <- clmplus:::place_forecast_triangle(forecast)
  expect_identical(dim(lower), c(4L, 4L))
  expect_equal(sum(!is.na(lower)), 6)
  expect_equal(lower[4, 2], forecast[2, 1])
  expect_equal(lower[2, 4], forecast[4, 1])
  observed <- matrix(1, 4, 4)
  expect_true(all(clmplus:::place_forecast_triangle(forecast, observed)[row(observed) + col(observed) <= 5] == 1))
})

test_that("small triangles and NA cells are retained", {
  x <- matrix(c(2, 3, 5, NA), 2, 2)
  expect_identical(clmplus:::calendar_to_triangle(clmplus:::triangle_to_calendar(x)), x)
})
