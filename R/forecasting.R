forecast_hazard_rates <- function(object, model_name, gk.fc.model = "a",
                                  ckj.fc.model = "a", gk.order = c(1, 1, 0),
                                  ckj.order = c(0, 1, 0)) {
  size <- ncol(object$Dxt)
  rates <- matrix(0, size, size)
  components <- model_components[[model_name]]
  kt.f <- gc.f <- NULL

  if (isTRUE(components[["cohort"]])) {
    last <- max(which(!is.na(object$gc)))
    if (gk.fc.model == "a") {
      fitted <- forecast::Arima(object$gc[seq_len(last)], order = gk.order,
                                include.constant = TRUE)
      gc.f <- forecast::forecast(fitted, h = length(object$cohorts) - last)
    } else if (gk.fc.model == "l") {
      fitted <- stats::lm(y ~ x, data = data.frame(
        y = object$gc[seq_len(last)], x = object$cohorts[seq_len(last)]
      ))
      gc.f <- forecast::forecast(fitted, newdata = data.frame(
        x = object$cohorts[(last + 1L):length(object$cohorts)]
      ))
    } else {
      stop("`gk.fc.model` must be either \"a\" or \"l\".", call. = FALSE)
    }
    cohort <- c(object$gc[!is.na(object$gc)], unname(gc.f$mean))
    cohort_matrix <- matrix(cohort, nrow = size, ncol = size)
    cohort_matrix <- triangle_to_full_calendar(cohort_matrix)
    cohort_matrix[is.na(cohort_matrix)] <- 0
    rates <- rates + cohort_matrix[, (size + 1L):(2L * size), drop = FALSE]
  }

  if (isTRUE(components[["period"]])) {
    period <- as.vector(object$kt[1L, ])
    last <- max(which(!is.na(period)))
    if (ckj.fc.model == "a") {
      fitted <- forecast::Arima(period[seq_len(last)], order = ckj.order,
                                include.constant = TRUE)
      kt.f <- forecast::forecast(fitted, h = size)
    } else if (ckj.fc.model == "l") {
      fitted <- stats::lm(y ~ x, data = data.frame(
        y = period[seq_len(last)], x = object$years[seq_len(last)]
      ))
      kt.f <- forecast::forecast(fitted, newdata = data.frame(x = seq.int(size + 1L, 2L * size)))
    } else {
      stop("`ckj.fc.model` must be either \"a\" or \"l\".", call. = FALSE)
    }
    rates <- sweep(rates, 2L, unname(kt.f$mean), `+`)
  }

  if (isTRUE(components[["age"]])) {
    age <- unname(object$ax)
    age[is.na(age)] <- 0
    rates <- sweep(rates, 1L, age, `+`)
  }

  list(rates = exp(rates), kt.f = kt.f, gc.f = gc.f)
}
