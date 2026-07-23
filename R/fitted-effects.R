fitted_development_factors <- function(size, age, period, cohort, eta) {
  age <- c(NA_real_, age[!is.na(age)])
  cohort <- if (is.null(cohort)) rep(0, size) else c(cohort[!is.na(cohort)], NA_real_)
  period <- if (is.null(period)) rep(0, size) else c(NA_real_, period[!is.na(period)])
  out <- matrix(NA_real_, size, size)
  index <- which(row(out) + col(out) <= size + 1L, arr.ind = TRUE)
  linear <- age[index[, 2L]] + cohort[index[, 1L]] + period[rowSums(index) - 1L]
  out[index] <- (1 + (1 - eta) * exp(linear)) / (1 - eta * exp(linear))
  out
}

extract_fitted_effects <- function(age, period, cohort, log_scale) {
  age <- c(NA_real_, age[!is.na(age)])
  names(age) <- seq_along(age) - 1L
  if (!is.null(period)) {
    period <- c(NA_real_, period[!is.na(period)])
    names(period) <- seq_along(period) - 1L
  }
  if (!is.null(cohort)) {
    cohort <- c(cohort[!is.na(cohort)], NA_real_)
    names(cohort) <- seq_along(cohort) - 1L
  }
  out <- list(
    fitted_development_effect = age,
    fitted_calendar_effect = period,
    fitted_accident_effect = cohort
  )
  if (!log_scale) out <- lapply(out, exp)
  out
}
