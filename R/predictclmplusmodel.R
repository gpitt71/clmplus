#' Predict the Reserve using Chain Ladder Plus Models
#'
#' Predict the lower triangle with a \code{clmplus} model.
#' 
#' @param object \code{clmplusmodel}, Model to predict from. 
#' @param gk.fc.model \code{character}, model to forecast the cohort component for the last accident period. It can be either arima ('a') or linear model ('l'). Disregarded for models that do not have a cohort effect.
#' @param ckj.fc.model \code{character}, model to forecast the calendar period effect. It can be either arima ('a') or linear model ('l'). Disregarded for models that do not have a period effect.
#' @param gk.order \code{integer}, order of the arima model with drift for the accident year effect extrapolation. Default to (1,1,0).
#' @param ckj.order \code{integer}, order of the arima model with drift for the calendar year effect extrapolation. Default to (0,1,0).
#' @param forecasting_horizon \code{integer}, between 1 and the triangle width. Calendar periods ahead for the predictions. Default predictions are to run-off. 
#' @param constrained_development_factors \code{logical}, if \code{TRUE} the predict function will set negative development factors to 1.
#' @param ... Extra arguments to be passed to the predict function.
#' 
#' @return Returns the following output:
#'   
#'   \item{reserve}{\code{numeric} The reserve for each accident period. }
#'   
#'   \item{ultimate_cost}{\code{numeric} The ultimate cost for each accident period. }
#'   
#'   \item{full_triangle}{\code{matrix array} The complete run-off triangle of cumulative payments, it includes the (input) upper triangle and the predicted (output) lower triangle.}
#'   
#'   \item{lower_triangle}{\code{matrix array} The predicted lower triangle of cumulative payments.}
#'   
#'   \item{development_factors_predicted}{\code{matrix array} The predicted lower triangle of the extrapolated development factors.}
#'   
#'   \item{apc_output}{\code{list} The following output from the age-period-cohort representation: \code{model.fit} (\code{fitStMoMo}) age-period-cohort model fit. 
#'   \code{alphaij} (\code{matrix array}) predicted claim development. 
#'   \code{lower_triangle_apc} (\code{matrix array}) predicted lower triangle of cumulative payments in age-period-cohort form.
#'   \code{development_factors_apc} (\code{matrix array}) development factors in age-period-cohort representation.}
#' 
#' 
#' @references 
#' Pittarello, Gabriele, Munir Hiabu, and Andrés M. Villegas. "Replicating and extending chain ladder via an age-period-cohort structure on the claim development in a run-off triangle." arXiv preprint arXiv:2301.03858 (2023).
#'  
#' @export
#' @method predict clmplusmodel
predict.clmplusmodel <- function(object,
                                 gk.fc.model='a',
                                 ckj.fc.model='a',
                                 gk.order=c(1,1,0),
                                 ckj.order=c(0,1,0),
                                 forecasting_horizon=NULL,
                                 constrained_development_factors=FALSE,
                                 ...){
  if (!inherits(object, "clmplusmodel")) {
    stop("`object` must be a `clmplusmodel` object.", call. = FALSE)
  }
  size <- object$apc_input$J
  horizon <- validate_horizon(forecasting_horizon, size)
  model <- object$model.fit
  eta <- object$apc_input$eta
  model_name <- object$apc_input$hazard.model
  if (!model_name %in% names(supported_models)) {
    stop("The fitted object's hazard model is not supported.", call. = FALSE)
  }

  alphaij <- forecast_hazard_rates(
    model, model_name, gk.fc.model, ckj.fc.model, gk.order, ckj.order
  )
  factors <- (1 + (1 - eta) * alphaij$rates) / (1 - eta * alphaij$rates)
  if (isTRUE(constrained_development_factors)) factors <- pmax(1, factors)

  forecast <- matrix(0, size, horizon)
  forecast[, 1L] <- c(0, object$apc_input$diagonal[seq_len(size - 1L)]) * factors[, 1L]
  if (horizon >= 2L) {
    for (calendar in 2L:horizon) {
      forecast[, calendar] <- c(0, forecast[seq_len(size - 1L), calendar - 1L]) *
        factors[, calendar]
    }
  }
  factors <- factors[, seq_len(horizon), drop = FALSE]
  observed_calendar <- triangle_to_calendar(object$apc_input$cumulative.payments.triangle)

  if (is.null(forecasting_horizon)) {
    ultimate_cost <- rev(c(rev(forecast[size, ]), observed_calendar[size, size]))
  } else {
    ultimate_cost <- rev(forecast[, horizon])
    ultimate_cost <- c(ultimate_cost[ultimate_cost == 0], ultimate_cost[ultimate_cost != 0])
    ultimate_cost[1L] <- observed_calendar[size, size]
    if (horizon >= 2L) ultimate_cost[2L:horizon] <- forecast[size, seq_len(horizon - 1L)]
  }
  reserve <- ultimate_cost - rev(observed_calendar[, size])
  names(reserve) <- names(ultimate_cost) <- seq_along(reserve) - 1L

  lower <- place_forecast_triangle(forecast)
  out <- list(
    reserve = reserve,
    ultimate_cost = ultimate_cost,
    full_triangle = place_forecast_triangle(forecast, object$apc_input$cumulative.payments.triangle),
    lower_triangle = lower,
    development_factors_predicted = place_forecast_triangle(factors),
    apc_output = list(
      model.fit = model, hazard.model = model_name, alphaij = alphaij,
      lower_triangle_apc = forecast, development_factors_apc = factors
    )
  )
  class(out) <- "clmpluspredictions"
  out
}

