#' Fit Chain Ladder plus on Run-off Triangles.
#'
#' Method to Estimate Chain Ladder plus models.
#' 
#' @param AggregateDataPP \code{AggregateDataPP} object, reverse time triangle to be fitted.
#' @param hazard.model \code{character}, hazard model supported from our package. The model can be chosen from:
#' \itemize{
#' \item{'a': Age model, this is equivalent to the Mack chain-ladder.}
#' \item{'ac': Age and cohort effects.}
#' \item{'ap': Age and cohort effects.}
#' \item{'apc': Age cohort and period effects.}
#' }
#' 
#' 
#' @param link \code{character}, defines the link function and random component associated with 
#'   the mortality model. \code{"log"} would assume that deaths follow a 
#'   Poisson distribution and use a log link while \code{"logit"} would assume 
#'   that deaths follow a Binomial distribution and a logit link.
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param staticAgeFun \code{logical}, indicates if a static age function 
#'   \eqn{\alpha_x} is to be included. To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param periodAgeFun \code{list}, a list of length \eqn{N} with the definitions of the 
#'   period age modulating parameters \eqn{\beta_x^{(i)}}. Each entry can take 
#'   values: \code{"NP"} for non-parametric age terms, \code{"1"} for 
#'   \eqn{\beta_x^{(i)}=1} or a predefined parametric function of 
#'   age (see details). Set this to \code{NULL} if there are no period terms 
#'   in the model.
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param cohortAgeFun \code{character} or \code{function}, defines the cohort age modulating parameter 
#'   \eqn{\beta_x^{(0)}}. It can take values: \code{"NP"} for non-parametric 
#'   age terms, \code{"1"} for \eqn{\beta_x^{(0)}=1}, a predefined parametric 
#'   function of age (see details) or \code{NULL} if there is no cohort effect. 
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param effect_log_scale \code{logical}, whether effects should be on the logarithmic scale. By default, \code{TRUE}.
#' @param constFun \code{function}, it defines the identifiability constraints of the 
#'   model. It must be a function of the form 
#'   \code{constFun <- function(ax, bx, kt, b0x, gc, wxt, ages)} taking a set
#'   of fitted model parameters and returning a list 
#'   \code{list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)}
#'   of the model parameters with the identifiability constraints applied. If 
#'   omitted no identifiability constraints are applied to the model.
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param ... parameters to be passed to clmplus.
#' 
#' @return No return value, called to pass method \code{clmplus.AggregateDataPP}. See \code{clmplus.AggregateDataPP} documentation.
#' 
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- AggregateDataPP(cumulative.payments.triangle=sifa.mtpl)
#' hz.chl=clmplus(sifa.mtpl.rtt, 'a')
#' 
#' @references 
#' Pittarello, Gabriele, Munir Hiabu, and Andrés M. Villegas. "Replicating and extending chain ladder 
#' via an age-period-cohort structure on the claim development in a run-off triangle." arXiv preprint arXiv:2301.03858 (2023).
#' 
#' @export
clmplus <- function(AggregateDataPP,
                    hazard.model=NULL,
                    link = c("log", "logit"), 
                    staticAgeFun = TRUE, 
                    periodAgeFun = "NP",
                    cohortAgeFun = NULL, 
                    effect_log_scale=TRUE,
                    constFun = function(ax, bx, kt, b0x, gc, wxt, ages) list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc),
                    ...){
  
  UseMethod("clmplus")}

#' Fit Chain Ladder Plus to reverse time triangles.
#' 
#' Default method to fit Chain Ladder plus models.
#' 
#' @param AggregateDataPP \code{AggregateDataPP} object, reverse time triangle to be fitted.
#' @param hazard.model \code{character}, hazard model supported from our package. The model can be chosen from:
#' \itemize{
#' \item{'a': Age model, this is equivalent to the Mack chain-ladder.}
#' \item{'ac': Age and cohort effects.}
#' \item{'ap': Age and cohort effects.}
#' \item{'apc': Age cohort and period effects.}
#' }
#' 
#' @param link \code{character}, defines the link function and random component associated with 
#'   the mortality model. \code{"log"} would assume that deaths follow a 
#'   Poisson distribution and use a log link while \code{"logit"} would assume 
#'   that deaths follow a Binomial distribution and a logit link.
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param staticAgeFun \code{logical}, indicates if a static age function 
#'   \eqn{\alpha_x} is to be included. To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param periodAgeFun \code{list}, a list of length \eqn{N} with the definitions of the 
#'   period age modulating parameters \eqn{\beta_x^{(i)}}. Each entry can take 
#'   values: \code{"NP"} for non-parametric age terms, \code{"1"} for 
#'   \eqn{\beta_x^{(i)}=1} or a predefined parametric function of 
#'   age (see details). Set this to \code{NULL} if there are no period terms 
#'   in the model.
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param cohortAgeFun \code{character} or \code{function}, defines the cohort age modulating parameter 
#'   \eqn{\beta_x^{(0)}}. It can take values: \code{"NP"} for non-parametric 
#'   age terms, \code{"1"} for \eqn{\beta_x^{(0)}=1}, a predefined parametric 
#'   function of age (see details) or \code{NULL} if there is no cohort effect. 
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param effect_log_scale \code{logical}, whether effects should be on the logarithmic scale. By default, \code{TRUE}.
#' @param constFun \code{function}, it defines the identifiability constraints of the 
#'   model. It must be a function of the form 
#'   \code{constFun <- function(ax, bx, kt, b0x, gc, wxt, ages)} taking a set
#'   of fitted model parameters and returning a list 
#'   \code{list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)}
#'   of the model parameters with the identifiability constraints applied. If 
#'   omitted no identifiability constraints are applied to the model.
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param ... parameters to be passed to clmplus.
#' 
#' @return No return value, called to pass method \code{clmplus.AggregateDataPP}. See \code{clmplus.AggregateDataPP} documentation.
#' 
#' @references 
#' Pittarello, Gabriele, Munir Hiabu, and Andrés M. Villegas. "Replicating and extending chain ladder 
#' via an age-period-cohort structure on the claim development in a run-off triangle." arXiv preprint arXiv:2301.03858 (2023).
#'  
#' Hiabu, Munir. “On the relationship between classical chain ladder and granular reserving.” 
#' Scandinavian Actuarial Journal 2017 (2017): 708 - 729.
#' 
#' @export
clmplus.default <- function(AggregateDataPP,
                            hazard.model=NULL,
                            link = c("log", "logit"), 
                            staticAgeFun = TRUE, 
                            periodAgeFun = "NP",
                            cohortAgeFun = NULL, 
                            effect_log_scale=TRUE,
                            constFun = function(ax, bx, kt, b0x, gc, wxt, ages) list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc),
                            ...){stop("`AggregateDataPP` must inherit from class \"AggregateDataPP\".", call. = FALSE)}

#' Fit Chain Ladder Plus to reverse time triangles.
#'
#' Method to fit Chain Ladder plus models to \code{AggregateDataPP} objects.
#' 
#' @param AggregateDataPP \code{AggregateDataPP} object, reverse time triangle to be fitted.
#' @param hazard.model \code{character}, hazard model supported from our package. The model can be chosen from:
#' \itemize{
#' \item{'a': Age model, this is equivalent to the Mack chain-ladder.}
#' \item{'ac': Age and cohort effects.}
#' \item{'ap': Age and cohort effects.}
#' \item{'apc': Age cohort and period effects.}
#' }
#' 
#' @param link \code{character}, defines the link function and random component associated with 
#'   the mortality model. \code{"log"} would assume that deaths follow a 
#'   Poisson distribution and use a log link while \code{"logit"} would assume 
#'   that deaths follow a Binomial distribution and a logit link.
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param staticAgeFun \code{logical}, indicates if a static age function 
#'   \eqn{\alpha_x} is to be included. To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param periodAgeFun \code{list}, a list of length \eqn{N} with the definitions of the 
#'   period age modulating parameters \eqn{\beta_x^{(i)}}. Each entry can take 
#'   values: \code{"NP"} for non-parametric age terms, \code{"1"} for 
#'   \eqn{\beta_x^{(i)}=1} or a predefined parametric function of 
#'   age (see details). Set this to \code{NULL} if there are no period terms 
#'   in the model.
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param cohortAgeFun \code{character} or \code{function}, defines the cohort age modulating parameter 
#'   \eqn{\beta_x^{(0)}}. It can take values: \code{"NP"} for non-parametric 
#'   age terms, \code{"1"} for \eqn{\beta_x^{(0)}=1}, a predefined parametric 
#'   function of age (see details) or \code{NULL} if there is no cohort effect. 
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param effect_log_scale \code{logical}, whether effects should be on the logarithmic scale. By default, \code{TRUE}.
#' @param constFun \code{function}, it defines the identifiability constraints of the 
#'   model. It must be a function of the form 
#'   \code{constFun <- function(ax, bx, kt, b0x, gc, wxt, ages)} taking a set
#'   of fitted model parameters and returning a list 
#'   \code{list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)}
#'   of the model parameters with the identifiability constraints applied. If 
#'   omitted no identifiability constraints are applied to the model.
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param ... parameters to be passed to clmplus.
#' 
#' @return An object of class \code{clmplusmodel}. A list with the following elements:
#'   \item{model.fit}{\code{fitStMoMo} object, specified hazard model fit from StMoMo.}
#'   
#'   \item{apc_input}{\code{list} object. A list containing the following model inputs in age-period-cohort notation: \code{J} (\code{integer}) Run-off triangle dimension.  \code{eta} (\code{numeric}) Expected time-to-event in the cell. I.e., lost exposure.   
#'   \code{diagonal} (\code{numeric}) Cumulative payments last diagonal. \code{hazard.model} (\code{character}), hazard model specified from the user. Set to \code{user.specific} when a custom model is passed. 
#'   }   
#'   
#'   \item{hazard_scaled_deviance_residuals}{ \code{matrix array} Triangle of the scaled deviance residuals. }
#'   
#'   \item{fitted_development_factors}{ \code{matrix array} Triangle of the fitted development factors. }
#'   
#'   \item{fitted_effects}{ \code{list} List of the development-accident-calendar effects fitted. }
#'    
#'   
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- AggregateDataPP(cumulative.payments.triangle=sifa.mtpl)
#' hz.chl=clmplus(sifa.mtpl.rtt, 'a')
#' 
#' @references 
#' Pittarello, Gabriele, Munir Hiabu, and Andrés M. Villegas. "Replicating and extending chain ladder 
#' via an age-period-cohort structure on the claim development in a run-off triangle." arXiv preprint arXiv:2301.03858 (2023).
#' 
#' @export
clmplus.AggregateDataPP <- function(AggregateDataPP,
                            hazard.model=NULL,
                            link = c("log", "logit"), 
                            staticAgeFun = TRUE, 
                            periodAgeFun = "NP",
                            cohortAgeFun = NULL, 
                            effect_log_scale=TRUE,
                            constFun = function(ax, bx, kt, b0x, gc, wxt, ages) list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc),
                            ...){
  
  
  if (length(hazard.model) != 1L || !is.character(hazard.model) ||
      !hazard.model %in% names(supported_models)) {
    stop("`hazard.model` must be one of: ", paste(names(supported_models), collapse = ", "), ".", call. = FALSE)
  }
  model <- StMoMo::fit(supported_models[[hazard.model]],
                       Dxt = AggregateDataPP$occurrance, 
                       Ext = AggregateDataPP$exposure,
                       wxt=AggregateDataPP$fit.w,
                       iterMax=as.integer(1e+05))
  

  
  #forecasting horizon
  J=dim(AggregateDataPP$cumulative.payments.triangle)[2]
  
  # Find fitted development factors
  fij.fit <- fitted_development_factors(
    J, age = model$ax, cohort = model$gc, period = model$kt,
    eta = AggregateDataPP$eta
  )
  
  res.m <- stats::residuals(model)
  res.tr <- calendar_to_triangle(res.m$residuals)
  
  fitted_effects <- extract_fitted_effects(
    age = model$ax, cohort = model$gc, period = model$kt,
    log_scale = effect_log_scale
  )
  
  
  
  out <- list(model.fit=model,
              apc_input=list(J=J,
                              eta=AggregateDataPP$eta,
                              hazard.model=hazard.model,
                              diagonal=AggregateDataPP$diagonal,
                              cumulative.payments.triangle=AggregateDataPP$cumulative.payments.triangle
                              ),
              hazard_scaled_deviance_residuals=res.tr,
              fitted_development_factors = fij.fit, 
              fitted_effects=fitted_effects)
  
  class(out) <- c('clmplusmodel')
    
  out
  
}
