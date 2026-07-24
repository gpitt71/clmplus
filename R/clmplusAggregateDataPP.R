#' Fit a Chain Ladder Plus hazard model
#'
#' Fits one of the package's age, age-cohort, age-period, or
#' age-period-cohort claim-development models to data prepared by
#' [AggregateDataPP()]. Estimation is performed by [StMoMo::fit.StMoMo()].
#'
#' @param AggregateDataPP An object created by [AggregateDataPP()]. It contains
#'   a square cumulative paid-claims triangle and the corresponding
#'   development-calendar occurrence, exposure, and weight matrices.
#' @param hazard.model A required character scalar selecting `"a"` (age only,
#'   equivalent to chain ladder), `"ac"` (age-cohort), `"ap"` (age-period), or
#'   `"apc"` (age-period-cohort).
#' @param link,staticAgeFun,periodAgeFun,cohortAgeFun,constFun Compatibility
#'   arguments retained from the original interface. The package's four
#'   built-in model definitions determine these settings, so these arguments
#'   are currently ignored.
#' @param effect_log_scale A logical scalar. If `TRUE` (the default), fitted
#'   effects are returned on the linear-predictor/log scale; if `FALSE`, they
#'   are exponentiated.
#' @param verbose A logical scalar passed to [StMoMo::fit.StMoMo()]. The default
#'   `FALSE` hides StMoMo fitting progress. `TRUE` displays progress, including
#'   zero-weighted ages, years, and cohorts and the start/finish of the gnm fit.
#' @param ... Reserved for future extensions; no arguments are currently
#'   forwarded.
#'
#' @return A `clmplusmodel` list with:
#' \describe{
#'   \item{model.fit}{The underlying `fitStMoMo` object. Its fitted `ax`, `kt`,
#'   and `gc` fields contain the selected age, period, and cohort effects;
#'   inapplicable effects are `NULL`. Other fields are supplied by StMoMo and
#'   should be treated as implementation details.}
#'   \item{apc_input}{A list containing `J` (triangle dimension), `eta`
#'   (within-cell exposure timing), `hazard.model`, `diagonal` (latest observed
#'   cumulative payments by calendar representation), and the original
#'   `cumulative.payments.triangle`.}
#'   \item{hazard_scaled_deviance_residuals}{A `J` by `J` numeric matrix in
#'   accident-year by development-year triangle orientation. Unobserved cells
#'   are `NA`.}
#'   \item{fitted_development_factors}{A `J` by `J` numeric matrix of fitted
#'   multiplicative cumulative development factors; unavailable cells are
#'   `NA`.}
#'   \item{fitted_effects}{A list with `fitted_development_effect`,
#'   `fitted_calendar_effect`, and `fitted_accident_effect`. Components not
#'   included in the selected model are `NULL`.}
#' }
#'
#' @details Incremental payment amounts can be non-integer even though the
#'   StMoMo fit uses a Poisson quasi-likelihood. Warnings whose messages begin
#'   exactly with `non-integer x =` are therefore expected and are selectively
#'   muffled. All other warnings, including convergence and numerical warnings,
#'   remain visible.
#'
#' @examples
#' data(sifa.mtpl)
#' prepared <- AggregateDataPP(sifa.mtpl)
#' age_fit <- clmplus(prepared, hazard.model = "a", verbose = FALSE)
#' age_fit$fitted_effects
#' \donttest{
#' apc_fit <- clmplus(prepared, hazard.model = "apc", verbose = FALSE)
#' plot(apc_fit)
#' }
#'
#' @seealso [AggregateDataPP()], [predict.clmplusmodel()],
#'   [predictReserve.clmplusmodel()], [plot.clmplusmodel()]
#' @export
clmplus <- function(AggregateDataPP, hazard.model = NULL,
                    link = c("log", "logit"), staticAgeFun = TRUE,
                    periodAgeFun = "NP", cohortAgeFun = NULL,
                    effect_log_scale = TRUE, verbose = FALSE,
                    constFun = function(ax, bx, kt, b0x, gc, wxt, ages) {
                      list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
                    }, ...) {
  UseMethod("clmplus")
}

#' @rdname clmplus
#' @return The default method always raises an informative error because
#'   `AggregateDataPP` does not inherit from `"AggregateDataPP"`.
#' @export
clmplus.default <- function(AggregateDataPP, hazard.model = NULL,
                            link = c("log", "logit"), staticAgeFun = TRUE,
                            periodAgeFun = "NP", cohortAgeFun = NULL,
                            effect_log_scale = TRUE, verbose = FALSE,
                            constFun = function(ax, bx, kt, b0x, gc, wxt, ages) {
                              list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
                            }, ...) {
  stop("`AggregateDataPP` must inherit from class \"AggregateDataPP\".",
       call. = FALSE)
}

#' @rdname clmplus
#' @export
clmplus.AggregateDataPP <- function(AggregateDataPP, hazard.model = NULL,
                                    link = c("log", "logit"),
                                    staticAgeFun = TRUE,
                                    periodAgeFun = "NP",
                                    cohortAgeFun = NULL,
                                    effect_log_scale = TRUE,
                                    verbose = FALSE,
                                    constFun = function(ax, bx, kt, b0x, gc,
                                                        wxt, ages) {
                                      list(ax = ax, bx = bx, kt = kt,
                                           b0x = b0x, gc = gc)
                                    }, ...) {
  if (length(hazard.model) != 1L || !is.character(hazard.model) ||
      !hazard.model %in% names(supported_models)) {
    stop("`hazard.model` must be one of: ",
         paste(names(supported_models), collapse = ", "), ".", call. = FALSE)
  }
  if (length(verbose) != 1L || !is.logical(verbose) || is.na(verbose)) {
    stop("`verbose` must be a single non-missing logical value.", call. = FALSE)
  }

  model <- withCallingHandlers(
    StMoMo::fit(
      supported_models[[hazard.model]],
      Dxt = AggregateDataPP$occurrance,
      Ext = AggregateDataPP$exposure,
      wxt = AggregateDataPP$fit.w,
      iterMax = as.integer(1e+05),
      verbose = verbose
    ),
    warning = function(w) {
      if (
        !isTRUE(verbose) ||
        grepl("^non-integer x\\s*=", conditionMessage(w))
      ) {
        invokeRestart("muffleWarning")
      }
    },
    message = function(m) {
      if (!isTRUE(verbose)) {
        invokeRestart("muffleMessage")
      }
    }
  )

  J <- ncol(AggregateDataPP$cumulative.payments.triangle)
  fij.fit <- fitted_development_factors(
    J, age = model$ax, cohort = model$gc, period = model$kt,
    eta = AggregateDataPP$eta
  )
  residuals <- stats::residuals(model)
  residual_triangle <- calendar_to_triangle(residuals$residuals)
  effects <- extract_fitted_effects(
    age = model$ax, cohort = model$gc, period = model$kt,
    log_scale = effect_log_scale
  )

  out <- list(
    model.fit = model,
    apc_input = list(
      J = J, eta = AggregateDataPP$eta, hazard.model = hazard.model,
      diagonal = AggregateDataPP$diagonal,
      cumulative.payments.triangle =
        AggregateDataPP$cumulative.payments.triangle
    ),
    hazard_scaled_deviance_residuals = residual_triangle,
    fitted_development_factors = fij.fit,
    fitted_effects = effects
  )
  class(out) <- "clmplusmodel"
  out
}
