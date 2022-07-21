#' Fit the reserving model to the triangle in reverse time
#'
#' This function allows to fit the extended chain-ladder models to the run-off triangles.
#' 
#' @param RtTriangle RtTriangle object to be fitted
#' @param hazard.model hazard model supported from our package, must be provided as a string. The model can be choosen from:
#'\describe{
#'        \item{'a': Age model, this is equivalent to the Mack chain-ladder.}
#'        \item{'ac': Age and cohort effects.}
#'        \item{'ap': Age and cohort effects.}
#'        \item{'apc': Age cohort and period effects.}
#'        \item{'cbd': Cairns-Blake-Dowd mortality model (CBD).}
#'        \item{'m6': CBD with cohorts.}
#'        \item{'m7': CBD m7 extension.}
#'        \item{'m8': CBD m7 extension.}
#'        }
#' @param xc xc constant parameter to be set for the m8 model. Default to NULL.
#' @param iter.Max maximum number of iterations for the Newton-Rhapson algorithm. Ignored for other fitting procedures.
#' @param tolerance.max maximum tolerance of parameters difference for convergence for the Newton-Rhapson algorithm implementation.Ignored for other fitting procedures.
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- RtTriangle(cumulative.payments.triangle=sifa.mtpl)
#' hz.chl=xchl(sifa.mtpl.rtt, 'a')
#' 
#' @references 
#' 
#' Hiabu, Munir. “On the relationship between classical chain ladder and granular reserving.” 
#' Scandinavian Actuarial Journal 2017 (2017): 708 - 729.
#' 
#' @export
xchl <- function(RtTriangle,
                 hazard.model,
                 xc=NULL,
                 iter.max=1e+04,
                 tolerance.max=1e-06){
  
  UseMethod("xchl")}

xchl.default <- function(RtTriangle){message('The object provided must be of class RtTriangle')}

xchl.RtTriangle <- function(RtTriangle,
                            hazard.model,
                            xc = NULL,
                            iter.max=1e+04,
                            tolerance.max=1e-06){
  
  stopifnot(typeof(hazard.model)=="character")
  
  if(hazard.model=='m8'){
    pkg.env$models$m8= StMoMo::m8(link = c("log"),xc=xc)
  }
  
  if(hazard.model=='lc'){
    
    J=dim(RtTriangle$cumulative.payments.triangle)[2]
    
    model <- pkg.env$fit.lc.nr(RtTriangle,
                               iter.max=iter.max,
                               tolerance.max=tolerance.max)
    
    kt.fit=arima(model$kt,c(0,1,0))
    
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
  # 
  # if(hazard.model=='aacohort'){
  #   
  #   J=dim(RtTriangle$cumulative.payments.triangle)[2]
  #   
  #   ac.model <- StMoMo::fit(pkg.env$models[['ac']], 
  #                        Dxt = RtTriangle$occurrance, 
  #                        Ext = RtTriangle$exposure,
  #                        iterMax=as.integer(1e+05))
  #   
  #   model <- pkg.env$fit.aac.nr(data.T=RtTriangle,
  #                               ax.start=ac.model$ax,
  #                               gc.start=rev(ac.model$gc[(J):(2*J-1)]),
  #                               iter.max=iter.max,
  #                               tolerance.max=tolerance.max)
  #   
  #   gc.mx = pkg.env$t2c.full.square(matrix(rep(model$gc,
  #                      J),
  #                  byrow = F,
  #                  ncol=J))[,(J+1):(2*J)]
  #   
  #   bx.mx = matrix(rep(model$bx,
  #                      J),
  #                  byrow = F,
  #                  ncol=J)
  #   
  #   ax.mx = matrix(rep(model$ax,J),
  #                  byrow = F,
  #                  ncol=J)
  #   
  #   alphaij = exp(ax.mx+bx.mx*gc.mx)
  #   
  #   fij = (2+alphaij)/(2-alphaij)  
  #   
  #   d=RtTriangle$diagonal[1:(J-1)]
  #   
  #   # extrapolate the results
  #   lt=array(0.,c(J,J))
  #   lt[,1]=c(0.,d)*fij[,1]
  #   for(j in 2:J){lt[,j]=c(0.,lt[1:(J-1),(j-1)])*fij[,j]} 
  #   
  #   ot_=pkg.env$t2c(RtTriangle$cumulative.payments.triangle)
  #   ultimate_cost=c(rev(lt[J,1:(J-1)]),ot_[J,J])
  #   reserve=rev(ultimate_cost-ot_[,J])
  #   ultimate_cost=rev(ultimate_cost)
  #   
  #   converged=model$converged
  #   citer=model$citer
  #   
  # }
  
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
  
  class(out) <- c('xchlmodel')
    
  return(out)
}


setMethod("xchl", "RtTriangle", xchl.RtTriangle)

