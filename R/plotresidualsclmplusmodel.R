#' Plot the hazard model residuals 
#'
#' This function allows to plot the hazard model residuals on the triangle payments.
#' 
#' @param obj clmplusmodel object to be plotted.
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- RtTriangle(cumulative.payments.triangle=sifa.mtpl)
#' hz.chl<-clmplus(sifa.mtpl.rtt, 'a')
#' plotresiduals(hz.chl)
#' 
#' @references 
#' 
#' Hiabu, Munir. “On the relationship between classical chain ladder and granular reserving.” 
#' Scandinavian Actuarial Journal 2017 (2017): 708 - 729.
#' 
#' @export
plotresiduals <- function(clmplusmodel,
                          heat.lim=c(-2.5,2.5)){
  
  UseMethod("plotresiduals")}

plotresiduals.default <- function(clmplusmodel,
                                  heat.lim=c(-2.5,2.5)){message('The object provided must be of class clmplusmodel')}


plotresiduals.clmplusmodel <- function(clmplusmodel,
                                    heat.lim=c(-2.5,2.5)){
  
  if(clmplusmodel$hazard.model %in% names(pkg.env$models)){
    
    res.m = residuals(clmplusmodel$model.fit)
    res.tr=pkg.env$c2t(res.m$residuals)
    colnames(res.tr) <- rownames(res.tr) <- c(0:(dim(res.tr)[2]-1))
    longdf.no.0 = ChainLadder::as.LongTriangle(res.tr)
    
  }
  
  if(clmplusmodel$hazard.model == 'lc'){
    
    data.O= clmplusmodel$model.fit$data.T$occurrance
    data.E= clmplusmodel$model.fit$data.T$exposure
    ind <- is.na(data.E)
    W = matrix(1,nrow=dim(data.O)[1],ncol=dim(data.E)[2])
    W[ind]=0
    
    ax.mx = matrix(rep(clmplusmodel$model.fit$ax,dim(data.O)[2]),
                   byrow = F,
                   ncol=dim(data.O)[2])
    
    bx.mx = matrix(rep(clmplusmodel$model.fit$bx,
                       dim(data.O)[2]),
                   byrow = F,
                   ncol=dim(data.O)[2])
    
    kt.mx = matrix(rep(clmplusmodel$model.fit$kt,
                       dim(data.O)[1]),
                   byrow = T,
                   nrow=dim(data.O)[1])
    
    
    mu.mx = exp(ax.mx+bx.mx*kt.mx)
    
    data.O.h = data.E*mu.mx
    
    res <- array(NA, dim(W))
    
    res[!ind] <- 2 * W[!ind] * (data.O[!ind] * log(data.O[!ind] / data.O.h[!ind]) - (data.O[!ind] - data.O.h[!ind]))
    signRes <- sign(data.O - data.O.h)
    
    phi <- sum(res[!ind]) / ((dim(data.O)[1]*(dim(data.O)[1]+1))/2 - 3*dim(data.O)[1])
    res.m <- signRes * sqrt(abs(res) / phi) 
    
    res.tr=pkg.env$c2t(res.m)
    colnames(res.tr) <- rownames(res.tr) <- c(0:(dim(res.tr)[2]-1))
    longdf.no.0 = ChainLadder::as.LongTriangle(res.tr)

  }

  
  p_hm <- ggplot2::ggplot(data=longdf.no.0, ggplot2::aes(x=as.integer(dev), y=as.integer(origin))) + 
    ggplot2::geom_tile(ggplot2::aes(fill = value))+ggplot2::scale_y_reverse()+
    ggplot2::scale_fill_gradient2(name="model residuals", 
                         low="royalblue", 
                         mid="white", 
                         high="#a71429", 
                         midpoint=0, 
                         space="Lab", 
                         na.value="grey50", 
                         limits=heat.lim,
                         guide="colourbar")+
    ggplot2::labs(x="Development year", y="Accident year")+
    ggplot2::ggtitle(clmplusmodel$modelfamily)+
    ggplot2::theme(axis.title.x = ggplot2::element_text(size=8), axis.text.x  = ggplot2::element_text(size=7))+
    ggplot2::theme(axis.title.y = ggplot2::element_text(size=8), axis.text.y  = ggplot2::element_text(size=7))+
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "grey", colour = "grey", size = 2, linetype = "solid"),
          panel.grid = ggplot2::element_line(colour="grey")) + 
    NULL
  
    p_hm
  
}


setMethod("plotresiduals", "clmplusmodel",plotresiduals.clmplusmodel)








