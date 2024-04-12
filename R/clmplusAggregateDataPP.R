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
#' \item{'lc': Lee-Carter parameters: age and age-period interaction effects.}
#' \item{'cbd': Cairns-Blake-Dowd mortality model (CBD).}
#' \item{'m6': CBD with cohorts.}
#' \item{'m7': CBD m7 extension.}
#' \item{'m8': CBD m7 extension.}
#' }
#' 
#' 
#' @param xc \code{integer}, xc constant parameter to be set for the m8 model. Default to NULL.
#' @param iter.max \code{integer}, maximum number of iterations for the Newton-Rhapson algorithm. It will be ignored for other fitting procedures.
#' @param tolerance.max \code{integer}, maximum tolerance of parameters difference for convergence for the Newton-Rhapson algorithm implementation.Ignored for other fitting procedures.
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
                    xc = NULL,
                    iter.max=1e+04,
                    tolerance.max=1e-06,
                    link = c("log", "logit"), 
                    staticAgeFun = TRUE, 
                    periodAgeFun = "NP",
                    cohortAgeFun = NULL, 
                    constFun = function(ax, bx, kt, b0x, gc, wxt, ages) list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc),
                    # gk.fc.model='a',
                    # ckj.fc.model='a',
                    # gk.order=c(1,1,0),
                    # ckj.order=c(0,1,0),
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
#' \item{'lc': Lee-Carter parameters: age and age-period interaction effects.}
#' \item{'cbd': Cairns-Blake-Dowd mortality model (CBD).}
#' \item{'m6': CBD with cohorts.}
#' \item{'m7': CBD m7 extension.}
#' \item{'m8': CBD m7 extension.}
#' }
#' 
#' @param gk.fc.model \code{character}, model to forecast the cohort component for the last accident period. It can be either arima ('a') or linear model ('l'). Disregarded for models that do not have a cohort effect.
#' @param ckj.fc.model \code{character}, model to forecast the calendar period effect. It can be either arima ('a') or linear model ('l'). Disregarded for models that do not have a period effect.
#' @param gk.order \code{integer}, order of the arima model with drift for the accident year effect extrapolation. Default to (1,1,0).
#' @param ckj.order \code{integer}, order of the arima model with drift for the calendar year effect extrapolation. Default to (0,1,0).
#' 
#' 
#' @param xc \code{integer}, xc constant parameter to be set for the m8 model. Default to NULL.
#' @param iter.max \code{integer}, maximum number of iterations for the Newton-Rhapson algorithm. It will be ignored for other fitting procedures.
#' @param tolerance.max \code{integer}, maximum tolerance of parameters difference for convergence for the Newton-Rhapson algorithm implementation.Ignored for other fitting procedures.
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
                            xc = NULL,
                            iter.max=1e+04,
                            tolerance.max=1e-06,
                            link = c("log", "logit"), 
                            staticAgeFun = TRUE, 
                            periodAgeFun = "NP",
                            cohortAgeFun = NULL, 
                            constFun = function(ax, bx, kt, b0x, gc, wxt, ages) list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc),
                            # gk.fc.model='a',
                            # ckj.fc.model='a',
                            # gk.order=c(1,1,0),
                            # ckj.order=c(0,1,0),
                            ...){message('The object provided must be of class AggregateDataPP')}

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
#' \item{'lc': Lee-Carter parameters: age and age-period interaction effects.}
#' \item{'cbd': Cairns-Blake-Dowd mortality model (CBD).}
#' \item{'m6': CBD with cohorts.}
#' \item{'m7': CBD m7 extension.}
#' \item{'m8': CBD m7 extension.}
#' }
#' 
#' @param xc \code{integer}, xc constant parameter to be set for the m8 model. Default to NULL.
#' @param iter.max \code{integer}, maximum number of iterations for the Newton-Rhapson algorithm. It will be ignored for other fitting procedures.
#' @param tolerance.max \code{integer}, maximum tolerance of parameters difference for convergence for the Newton-Rhapson algorithm implementation.Ignored for other fitting procedures.
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
#'   \item{hazard.model}{\code{character}, hazard model specified from the user. Set to \code{user.specific} when a custom model is passed. }
#'      
#'   \item{ultimate.cost}{\code{numeric}, vector of predicted ultimate costs.}
#'   
#'   \item{model.fcst}{\code{list}, it contains \code{rates} (hazard rates forecasted on each cell of the triangle), \code{kt.f} (forecasted calendar period effect when present) and \code{gc.f} (forecasted cohort effect when present) }
#'   
#'   \item{converged}{\code{logical}, \code{TRUE} when the model fit converged.}
#'   
#'   \item{citer}{\code{numeric}, Number of Netwon-Rhapson iterations in case a lee-carter hazard-model was chosen.}
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
#' Hiabu, Munir. “On the relationship between classical chain ladder and granular reserving.” 
#' Scandinavian Actuarial Journal 2017 (2017): 708 - 729.
#' 
#' @export
clmplus.AggregateDataPP <- function(AggregateDataPP,
                            hazard.model=NULL,
                            xc = NULL,
                            iter.max=1e+04,
                            tolerance.max=1e-06,
                            link = c("log", "logit"), 
                            staticAgeFun = TRUE, 
                            periodAgeFun = "NP",
                            cohortAgeFun = NULL, 
                            constFun = function(ax, bx, kt, b0x, gc, wxt, ages) list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc),
                            # gk.fc.model='a',
                            # ckj.fc.model='a',
                            # gk.order=c(1,1,0),
                            # ckj.order=c(0,1,0),
                            ...){
  
  
  stopifnot(is.null(hazard.model) | typeof(hazard.model)=="character")
  
  if(is.null(hazard.model)){
    

    stmomo.model = StMoMo::StMoMo(link = link, 
                          staticAgeFun = staticAgeFun, 
                          periodAgeFun = periodAgeFun,
                          cohortAgeFun = cohortAgeFun, 
                          constFun = constFun)
    
    model <- StMoMo::fit(stmomo.model, 
                         Dxt = AggregateDataPP$occurrance, 
                         Ext = AggregateDataPP$exposure,
                         wxt = AggregateDataPP$fit.w,
                         iterMax=as.integer(1e+05))
    
    #forecasting horizon
    J=dim(AggregateDataPP$cumulative.payments.triangle)[2]
    #compute the development factors
    alphaij <- forecast::forecast(model, h = J)
    # fij=(2+alphaij$rates)/(2-alphaij$rates)
    fij=(1+(1-AggregateDataPP$k)*alphaij$rates)/(1-(AggregateDataPP$k*alphaij$rates))
    # pick the last diagonal
    d=AggregateDataPP$diagonal[1:(J-1)]
    # extrapolate the results
    lt=array(0.,c(J,J))
    lt[,1]=c(0.,d)*fij[,1]
    for(j in 2:J){lt[,j]=c(0.,lt[1:(J-1),(j-1)])*fij[,j]} 
    
    ot_=pkg.env$t2c(AggregateDataPP$cumulative.payments.triangle)
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
  
  
  
  # if(hazard.model=='m8'){
  #   pkg.env$models$m8= StMoMo::m8(link = c("log"),xc=xc)
  # }
  # 
  # if(hazard.model=='lc'){
  #   
  #   J=dim(AggregateDataPP$cumulative.payments.triangle)[2]
  #   
  #   model <- pkg.env$fit.lc.nr(AggregateDataPP,
  #                              iter.max=iter.max,
  #                              tolerance.max=tolerance.max)
  #   
  #   
  #   kt.nNA <- max(which(!is.na(model$kt)))
  #   
  #   if(ckj.fc.model=='a'){
  #     
  #     kt.model=forecast::Arima(as.vector(model$kt[1:kt.nNA]),c(0,1,0),include.constant = T)
  #     kt.f <- forecast::forecast(kt.model,h=J)
  #     
  #   }else{
  #     
  #     kt.data=data.frame(y=model$kt[1:kt.nNA],
  #                        x=1:kt.nNA)
  #     new.kt.data <- data.frame(x=seq(J+1,2*J))
  #     
  #     kt.model <- stats::lm('y~x', 
  #                           data=kt.data)
  #     
  #     kt.f <- forecast::forecast(kt.model,newdata=new.kt.data)
  #     
  #   }
  #   
  #   # model$kt.fcst=kt.f
  #   
  #   
  #   kt.mx = matrix(rep(kt.f$mean,
  #                      J),
  #                  byrow = T,
  #                  nrow=J)
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
  #   alphaij = exp(ax.mx+bx.mx*kt.mx)
  #   
  #   # fij = (2+alphaij)/(2-alphaij)  
  #   fij=(1+(1-AggregateDataPP$k)*alphaij)/(1-(AggregateDataPP$k*alphaij))
  #   
  #   d=AggregateDataPP$diagonal[1:(J-1)]
  #   
  #   # extrapolate the results
  #   lt=array(0.,c(J,J))
  #   lt[,1]=c(0.,d)*fij[,1]
  #   for(j in 2:J){lt[,j]=c(0.,lt[1:(J-1),(j-1)])*fij[,j]} 
  #   
  #   ot_=pkg.env$t2c(AggregateDataPP$cumulative.payments.triangle)
  #   ultimate_cost=c(rev(lt[J,1:(J-1)]),ot_[J,J])
  #   reserve=rev(ultimate_cost-ot_[,J])
  #   ultimate_cost=rev(ultimate_cost)
  #   
  #   converged=model$converged
  #   citer=model$citer
  #   
  #   alphaij = list(rates=alphaij,
  #                  kt.f=kt.f)
  #   
  # }
  # 
  # 
  if(hazard.model %in% names(pkg.env$models)){
    # browser()
  model <- StMoMo::fit(pkg.env$models[[hazard.model]], 
                       Dxt = AggregateDataPP$occurrance, 
                       Ext = AggregateDataPP$exposure,
                       wxt=AggregateDataPP$fit.w,
                       iterMax=as.integer(1e+05))
  
  #forecasting horizon
  J=dim(AggregateDataPP$cumulative.payments.triangle)[2]
  #compute the development factors
  # with stmomo:
 
  # alphaij <- pkg.env$fcst(model, 
  #                         hazard.model = hazard.model,
  #                         gk.fc.model=gk.fc.model,
  #                         ckj.fc.model=ckj.fc.model,
  #                         gk.order=gk.order,
  #                         ckj.order=ckj.order
  #                         )

  # fij=(1+(1-AggregateDataPP$k)*alphaij$rates)/(1-(AggregateDataPP$k*alphaij$rates))
 
   # pick the last diagonal
  # d=AggregateDataPP$diagonal[1:(J-1)]
  # extrapolate the results
  # lt=array(0.,c(J,J))
  # lt[,1]=c(0.,d)*fij[,1]
  # for(j in 2:J){lt[,j]=c(0.,lt[1:(J-1),(j-1)])*fij[,j]} 
  # 
  # ot_=pkg.env$t2c(AggregateDataPP$cumulative.payments.triangle)
  # ultimate_cost=c(rev(lt[J,1:(J-1)]),ot_[J,J])
  # reserve=rev(ultimate_cost-ot_[,J])
  # ultimate_cost=rev(ultimate_cost)
  converged=TRUE
  citer=NULL
  }
  
  out <- list(model.fit=model,
              apc_input=list(J=J,
                              eta=AggregateDataPP$k,
                              hazard.model=hazard.model,
                              diagonal=AggregateDataPP$diagonal,
                              cumulative.payments.triangle=AggregateDataPP$cumulative.payments.triangle
                              ),
                                
              # ultimate.cost=ultimate_cost,
              # reserve=reserve,
              # model.fcst = alphaij,
              converged=converged,
              citer=citer)
  
  class(out) <- c('clmplusmodel')
    
  out
  
}
