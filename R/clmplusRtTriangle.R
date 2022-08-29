#' Fit chain-ladder+ to reverse time triangles.
#'
#' Generic method to fit the chain ladder +.
#' 
#' @param RtTriangle RtTriangle object to be fitted.
#' @param hazard.model hazard model supported from our package, must be provided as a string. The model can be choosen from:
#' \itemize{
#' \item{'a': Age model, this is equivalent to the Mack chain-ladder.}
#' \item{'ac': Age and cohort effects.}
#' \item{'ap': Age and cohort effects.}
#' \item{'apc': Age cohort and period effects.}
#' \item{'cbd': Cairns-Blake-Dowd mortality model (CBD).}
#' \item{'m6': CBD with cohorts.}
#' \item{'m7': CBD m7 extension.}
#' \item{'m8': CBD m7 extension.}
#' }
#' @param ... arguments to be passed to or from other methods.
#' 
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- RtTriangle(cumulative.payments.triangle=sifa.mtpl)
#' hz.chl=clmplus(sifa.mtpl.rtt, 'a')
#' 
#' @references 
#' Hiabu, Munir. “On the relationship between classical chain ladder and granular reserving.” 
#' Scandinavian Actuarial Journal 2017 (2017): 708 - 729.
#' 
#' @export
clmplus <- function(RtTriangle,
                 hazard.model=NULL,
                 ...){
  
  UseMethod("clmplus")}

#' Fit chain-ladder+ to reverse time triangles.
#' 
#' This function allows to fit chain-ladder+ models to cumulative payments run-off triangles.
#' 
#' @param RtTriangle RtTriangle object to be fitted.
#' @param hazard.model hazard model supported from our package, must be provided as a string. The model can be choosen from:
#' \itemize{
#' \item{'a': Age model, this is equivalent to the Mack chain-ladder.}
#' \item{'ac': Age and cohort effects.}
#' \item{'ap': Age and cohort effects.}
#' \item{'apc': Age cohort and period effects.}
#' \item{'cbd': Cairns-Blake-Dowd mortality model (CBD).}
#' \item{'m6': CBD with cohorts.}
#' \item{'m7': CBD m7 extension.}
#' \item{'m8': CBD m7 extension.}
#' }
#' @param ... parameters to be passed to clmplus.
#' 
#' @references 
#' Hiabu, Munir. “On the relationship between classical chain ladder and granular reserving.” 
#' Scandinavian Actuarial Journal 2017 (2017): 708 - 729.
#' 
#' @export
clmplus.default <- function(RtTriangle,
                            hazard.model=NULL,
                            ...){message('The object provided must be of class RtTriangle')}

#' Fit chain-ladder+ to reverse time triangles.
#'
#' This function allows to fit chain-ladder+ models to cumulative payments run-off triangles.
#' 
#' @param RtTriangle RtTriangle object to be fitted.
#' @param hazard.model hazard model supported from our package, must be provided as a string. The model can be choosen from:
#' \itemize{
#' \item{'a': Age model, this is equivalent to the Mack chain-ladder.}
#' \item{'ac': Age and cohort effects.}
#' \item{'ap': Age and cohort effects.}
#' \item{'apc': Age cohort and period effects.}
#' \item{'cbd': Cairns-Blake-Dowd mortality model (CBD).}
#' \item{'m6': CBD with cohorts.}
#' \item{'m7': CBD m7 extension.}
#' \item{'m8': CBD m7 extension.}
#' }
#' @param xc xc constant parameter to be set for the m8 model. Default to NULL.
#' @param iter.max maximum number of iterations for the Newton-Rhapson algorithm. It will be ignored for other fitting procedures.
#' @param tolerance.max maximum tolerance of parameters difference for convergence for the Newton-Rhapson algorithm implementation.Ignored for other fitting procedures.
#' @param link defines the link function and random component associated with 
#'   the mortality model. \code{"log"} would assume that deaths follow a 
#'   Poisson distribution and use a log link while \code{"logit"} would assume 
#'   that deaths follow a Binomial distribution and a logit link.
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param staticAgeFun logical value indicating if a static age function 
#'   \eqn{\alpha_x} is to be included. To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param periodAgeFun a list of length \eqn{N} with the definitions of the 
#'   period age modulating parameters \eqn{\beta_x^{(i)}}. Each entry can take 
#'   values: \code{"NP"} for non-parametric age terms, \code{"1"} for 
#'   \eqn{\beta_x^{(i)}=1} or a predefined parametric function of 
#'   age (see details). Set this to \code{NULL} if there are no period terms 
#'   in the model.
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param cohortAgeFun defines the cohort age modulating parameter 
#'   \eqn{\beta_x^{(0)}}. It can take values: \code{"NP"} for non-parametric 
#'   age terms, \code{"1"} for \eqn{\beta_x^{(0)}=1}, a predefined parametric 
#'   function of age (see details) or \code{NULL} if there is no cohort effect. 
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param constFun function defining the identifiability constraints of the 
#'   model. It must be a function of the form 
#'   \code{constFun <- function(ax, bx, kt, b0x, gc, wxt, ages)} taking a set
#'   of fitted model parameters and returning a list 
#'   \code{list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)}
#'   of the model parameters with the identifiability constraints applied. If 
#'   omitted no identifiability constraints are applied to the model.
#'   To be disregarded unless the practitioner specifies his own hazard model in StMoMo. 
#' @param ... extra arguments to be passed from or to other methods.
#' 
#' @return An object of class \code{"clmplusmodel"} with components:
#'   \item{model.fit}{Hazard model fit from StMoMo.}
#'   
#'   \item{hazard.model}{Hazard model chosen.}
#'   
#'   \item{exposure}{Matrix that contains the exposure derived from the input triangle, under the uniform claims arrival assumption.}
#'   
#'   \item{ultimate.cost}{Ultimate costs vector.}
#'   
#'   \item{model.fcst}{Hazard forecasts.}
#'   
#'   \item{converged}{logical value. Whether the fit converged.}
#'   
#'   \item{citer}{Number of Netwon-Rhapson iterations in case a lee-carter hazard-model was chosen.}
#'   
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- RtTriangle(cumulative.payments.triangle=sifa.mtpl)
#' hz.chl=clmplus(sifa.mtpl.rtt, 'a')
#' 
#' @references 
#' Hiabu, Munir. “On the relationship between classical chain ladder and granular reserving.” 
#' Scandinavian Actuarial Journal 2017 (2017): 708 - 729.
#' 
#' @export
clmplus.RtTriangle <- function(RtTriangle,
                            hazard.model=NULL,
                            xc = NULL,
                            iter.max=1e+04,
                            tolerance.max=1e-06,
                            link = c("log", "logit"), 
                            staticAgeFun = TRUE, 
                            periodAgeFun = "NP",
                            cohortAgeFun = NULL, 
                            constFun = function(ax, bx, kt, b0x, gc, wxt, ages) list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc),
                            ...){
  
  stopifnot(is.null(hazard.model) | typeof(hazard.model)=="character")
  
  if(is.null(hazard.model)){
    
    stmomo.model = StMoMo::StMoMo(link = link, 
                          staticAgeFun = staticAgeFun, 
                          periodAgeFun = periodAgeFun,
                          cohortAgeFun = cohortAgeFun, 
                          constFun = constFun)
    
    model <- StMoMo::fit(stmomo.model, 
                         Dxt = RtTriangle$occurrance, 
                         Ext = RtTriangle$exposure,
                         iterMax=as.integer(1e+05))
    
    #forecasting horizon
    J=dim(RtTriangle$cumulative.payments.triangle)[2]
    #compute the development factors
    alphaij <- forecast::forecast(model, h = J)
    fij=(2+alphaij$rates)/(2-alphaij$rates)
    # pick the last diagonal
    d=RtTriangle$diagonal[1:(J-1)]
    # extrapolate the results
    lt=array(0.,c(J,J))
    lt[,1]=c(0.,d)*fij[,1]
    for(j in 2:J){lt[,j]=c(0.,lt[1:(J-1),(j-1)])*fij[,j]} 
    
    ot_=pkg.env$t2c(RtTriangle$cumulative.payments.triangle)
    ultimate_cost=c(rev(lt[J,1:(J-1)]),ot_[J,J])
    reserve=rev(ultimate_cost-ot_[,J])
    ultimate_cost=rev(ultimate_cost)
    converged=TRUE
    citer=NULL
    
  
   out <- list(model.fit=model,
              hazard.model='user.defined',
              ultimate.cost=ultimate_cost,
              reserve=reserve,
              model.fcst = alphaij,
              converged=converged,
              citer=citer)
  
  class(out) <- c('clmplusmodel')
  
  out}
  
  if(hazard.model=='m8'){
    pkg.env$models$m8= StMoMo::m8(link = c("log"),xc=xc)
  }
  
  if(hazard.model=='lc'){
    
    J=dim(RtTriangle$cumulative.payments.triangle)[2]
    
    model <- pkg.env$fit.lc.nr(RtTriangle,
                               iter.max=iter.max,
                               tolerance.max=tolerance.max)
    
    kt.fit=stats::arima(model$kt,c(0,1,0))
    
    kt.fcst = forecast::forecast(kt.fit,
                       J)
    
    kt.mx = matrix(rep(kt.fcst$mean,
                       J),
                   byrow = T,
                   nrow=J)
    
    bx.mx = matrix(rep(model$bx,
                       J),
                   byrow = F,
                   ncol=J)
    
    ax.mx = matrix(rep(model$ax,J),
                   byrow = F,
                   ncol=J)
    
    alphaij = exp(ax.mx+bx.mx*kt.mx)
    
    fij = (2+alphaij)/(2-alphaij)  
    
    d=RtTriangle$diagonal[1:(J-1)]
    
    # extrapolate the results
    lt=array(0.,c(J,J))
    lt[,1]=c(0.,d)*fij[,1]
    for(j in 2:J){lt[,j]=c(0.,lt[1:(J-1),(j-1)])*fij[,j]} 
    
    ot_=pkg.env$t2c(RtTriangle$cumulative.payments.triangle)
    ultimate_cost=c(rev(lt[J,1:(J-1)]),ot_[J,J])
    reserve=rev(ultimate_cost-ot_[,J])
    ultimate_cost=rev(ultimate_cost)
    
    converged=model$converged
    citer=model$citer
    
  }

  
  if(hazard.model %in% names(pkg.env$models)){
  
  model <- StMoMo::fit(pkg.env$models[[hazard.model]], 
                       Dxt = RtTriangle$occurrance, 
                       Ext = RtTriangle$exposure,
                       iterMax=as.integer(1e+05))
  
  #forecasting horizon
  J=dim(RtTriangle$cumulative.payments.triangle)[2]
  #compute the development factors
  alphaij <- forecast::forecast(model, h = J)
  fij=(2+alphaij$rates)/(2-alphaij$rates)
  # pick the last diagonal
  d=RtTriangle$diagonal[1:(J-1)]
  # extrapolate the results
  lt=array(0.,c(J,J))
  lt[,1]=c(0.,d)*fij[,1]
  for(j in 2:J){lt[,j]=c(0.,lt[1:(J-1),(j-1)])*fij[,j]} 
  
  ot_=pkg.env$t2c(RtTriangle$cumulative.payments.triangle)
  ultimate_cost=c(rev(lt[J,1:(J-1)]),ot_[J,J])
  reserve=rev(ultimate_cost-ot_[,J])
  ultimate_cost=rev(ultimate_cost)
  converged=TRUE
  citer=NULL
  }
  
  out <- list(model.fit=model,
             hazard.model=hazard.model,
             ultimate.cost=ultimate_cost,
             reserve=reserve,
             model.fcst = alphaij,
             converged=converged,
             citer=citer)
  
  class(out) <- c('clmplusmodel')
    
  out
  
}
