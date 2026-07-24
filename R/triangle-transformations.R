triangle_to_calendar <- function(x) {
  nr <- nrow(x)
  nc <- ncol(x)
  out <- matrix(NA_real_, nr, nc)
  index <- which(row(x) + col(x) <= nc + 1L, arr.ind = TRUE)
  out[cbind(index[, 2L], rowSums(index) - 1L)] <- x[index]
  out
}

calendar_to_triangle <- function(x) {
  nr <- nrow(x)
  nc <- ncol(x)
  out <- matrix(NA_real_, nr, nc)
  index <- which(row(out) + col(out) <= nc + 1L, arr.ind = TRUE)
  out[index] <- x[cbind(index[, 2L], rowSums(index) - 1L)]
  out
}

triangle_to_full_calendar <- function(x) {
  nr <- nrow(x)
  nc <- ncol(x)
  out <- matrix(NA_real_, nr, 2L * nc)
  index <- which(matrix(TRUE, nr, nc), arr.ind = TRUE)
  out[cbind(index[, 2L], rowSums(index) - 1L)] <- x[index]
  out
}

place_forecast_triangle <- function(forecast, observed = NULL) {
  size <- nrow(forecast)
  out <- if (is.null(observed)) matrix(NA_real_, size, size) else observed
  index <- which(row(out) + col(out) > size + 1L, arr.ind = TRUE)
  calendar <- rowSums(index) - size - 1L
  valid <- calendar >= 1L & calendar <= ncol(forecast)
  index <- index[valid, , drop = FALSE]
  calendar <- calendar[valid]
  out[index] <- forecast[cbind(index[, 2L], calendar)]
  out
}
