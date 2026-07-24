validate_triangle <- function(x, weights = NULL, eta = 0.5) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x) || nrow(x) < 2L) {
    stop("`cumulative.payments.triangle` must be a square numeric matrix with at least two rows.", call. = FALSE)
  }
  if (length(eta) != 1L || !is.finite(eta) || eta <= 0 || eta > 1) {
    stop("`eta` must be one finite number in (0, 1].", call. = FALSE)
  }
  if (any(x < 0, na.rm = TRUE)) {
    stop("Please provide an input triangle without negative values.", call. = FALSE)
  }
  incremental <- x
  incremental[, -1L] <- x[, -1L, drop = FALSE] - x[, -ncol(x), drop = FALSE]
  if (any(incremental < 0, na.rm = TRUE)) {
    stop("Please provide an input triangle of cumulative payments. Payments recoveries are not allowed in our model framework.", call. = FALSE)
  }
  if (!is.null(weights) && (!is.matrix(weights) || !identical(dim(weights), dim(x)) ||
                            !is.numeric(weights) || any(!is.finite(weights[!is.na(weights)])) ||
                            any(weights < 0, na.rm = TRUE))) {
    stop("`entries.weights` must be a non-negative numeric matrix with the same dimensions as the triangle.", call. = FALSE)
  }
  incremental
}

validate_horizon <- function(horizon, size) {
  if (is.null(horizon)) return(size - 1L)
  if (length(horizon) != 1L || is.na(horizon) || horizon != as.integer(horizon) ||
      horizon < 1L || horizon > size - 1L) {
    stop("`forecasting_horizon` must be an integer between 1 and J - 1.", call. = FALSE)
  }
  as.integer(horizon)
}
