#' Model selection and scoring
#'
#' This function allows to fit the extended chain-ladder models to the run-off triangles. 
#' By doing so, different model fits are compared in terms of goodness of fit.
#' @param RtTriangle Reverse time triangle on which model performance is evaluated.
#' 
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- RtTriangle(cumulative.payments.triangle=sifa.mtpl)
#' out<-modelss(sifa.mtpl.rtt)
#' 
#' @references 
#' 
#' Hiabu, Munir. “On the relationship between classical chain ladder and granular reserving.” 
#' Scandinavian Actuarial Journal 2017 (2017): 708 - 729.
#' 
#' @export
modelss <- function(RtTriangle){
  
  UseMethod("modelss")
  
  
}

modelss.default <- function(RtTriangle){
  
  cat('You must provide an object belonging to the RtTriangle class.')
  
  
}

modelss.RtTriangle <- function(RtTriangle){
  
  leave.out=2
  
  aic = NULL
  bic = NULL
  rmse = NULL
  mae=NULL
  error.pc = NULL
  model.name = NULL
  error.incidence = NULL
  mre = NULL
  
  #pre-processing
  triangle <- RtTriangle$cumulative.payments.triangle
  J <- dim(triangle)[2]
  reduced.triangle <- pkg.env$c2t(pkg.env$t2c(triangle)[1:(J-leave.out),1:(J-leave.out)])
  newt.rtt <- RtTriangle(reduced.triangle)
  # to.project <- pkg.env$t2c(triangle)[1:(J-leave.out-1),J-leave.out]
  # true.values <- pkg.env$t2c(triangle)[2:(J-leave.out),J]
  to.project <- pkg.env$t2c(triangle)[1:(J-leave.out-1),J-leave.out]
  true.values <- pkg.env$t2c(triangle)[2:(J-leave.out),(J-leave.out+1):J]
  
  for(ix in names(pkg.env$models)){
    
    # hz.fit = clmplus(newt.rtt, ix)$model.fit
    # hz.rate = forecast::forecast(hz.fit,h=1)$rates
    
    hz.fit <- StMoMo::fit(pkg.env$models[[ix]], 
                          Dxt = newt.rtt$occurrance, 
                          Ext = newt.rtt$exposure,
                          iterMax=as.integer(1e+05))
    hz.rate = forecast::forecast(hz.fit,h=leave.out)$rates
    
    # if(leave.out==2){
      
    J.new=dim(reduced.triangle)[2]
    fij = (2+hz.rate)/(2-hz.rate)
    pred.mx = fij
    pred.mx[,1]=fij[,1]*c(NA,to.project)
    temp=c(triangle[(J-leave.out)+1,1],unname(pred.mx[1:(J.new-1),1][!is.na(pred.mx[1:(J.new-1),1])]))
    pred.mx[,2]=fij[,2]*c(rep(NA,J.new-length(temp)),temp)
    true.mx= rbind(rep(NA,2),true.values)
    
    sq.errors = (pred.mx-true.mx)^2
    abs.errors = abs(pred.mx-true.mx)
    r.errors = (pred.mx-true.mx)/true.mx
    error.inc.num = apply(pred.mx-true.mx,sum,MARGIN=2,na.rm=T)
    error.inc.den = apply(true.mx,sum,MARGIN=2,na.rm=T)
    model.name.ix = c(paste0(ix,".val"),paste0(ix,".test"))
    
    model.name = c(model.name,model.name.ix)
    rmse = c(rmse,sqrt(apply(sq.errors,MARGIN = 2,mean,na.rm=T)))
    mae = c(mae,apply(abs.errors,MARGIN = 2,mean,na.rm=T))
    aic = c(aic,rep(AIC(hz.fit),2))
    bic = c(bic,rep(BIC(hz.fit),2))
    mre = c(mre,apply(r.errors,MARGIN = 2,mean,na.rm=T))
    error.incidence = c(error.incidence,error.inc.num/error.inc.den)
    
    # }else{
    #   fij = unname((2+hz.rate)/(2-hz.rate))[2:(J-leave.out)]
    #   model.name = c(model.name,ix)
    #   rmse = c(rmse,sqrt(mean((to.project*fij-true.values)^2)))
    #   mae = c(mae,sqrt(mean(abs(to.project*fij-true.values))))
    #   aic = c(aic,AIC(hz.fit))
    #   bic = c(bic,BIC(hz.fit))
    #   predicted.inc=to.project*fij-to.project
    #   actual.inc=true.values-to.project
    #   mre = c(mre,mean((to.project*fij-true.values)/true.values))
    #   error.incidence = c(error.incidence,sum(to.project*fij-true.values)/sum(true.values))
    # }
  }
  

  out1 <- data.frame(
    model.name,
    mre,
    error.incidence,
    rmse,
    mae,
    aic,
    bic)
  
  temp.ix <- grepl(".val", model.name)
  temp.df <- out1[temp.ix,]
  
  
  out2 <- data.frame(aic=temp.df$model.name[which(temp.df$aic==min(temp.df$aic))],
                     bic=temp.df$model.name[which(temp.df$bic==min(temp.df$bic))],
                     rmse=temp.df$model.name[which(abs(temp.df$rmse)==min(abs(temp.df$rmse)))],
                     mre=temp.df$model.name[which(abs(temp.df$mre)==min(abs(temp.df$mre)))],
                     mae=temp.df$model.name[which(abs(temp.df$mae)==min(abs(temp.df$mae)))],
                     error.incidence=temp.df$model.name[which(abs(temp.df$error.incidence)==min(abs(temp.df$error.incidence)))])
  
  temp.ix <- grepl(".test", model.name)
  out3 <- out1[temp.ix,]
  
  out <- list(modelss=out1,
              val.best=out2,
              modelss.test=out3)
  
  class(out) <- c("modelss")
  
  return(out)
  
}

setMethod("modelss", "RtTriangle", modelss.RtTriangle)





