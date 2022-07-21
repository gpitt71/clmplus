#' Plot the hazard model fitted and forecasted parameters
#'
#' This function allows to define the behavior of the triangle payments.
#' 
#' @param obj clmplus model to be plotted.
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- RtTriangle(cumulative.payments.triangle=sifa.mtpl)
#' hz.chl<-clmplus(sifa.mtpl.rtt, 'a')
#' plot(hz.chl)
#' 
#' @references 
#' 
#' Hiabu, Munir. “On the relationship between classical chain ladder and granular reserving.” 
#' Scandinavian Actuarial Journal 2017 (2017): 708 - 729.
#' 
#' @export
plot.clmplusmodel <- function(obj){
  
  if(!is.null(obj$model.fit$ax)){
    ax=obj$model.fit$ax
    df.fitted= data.frame(dy=obj$model.fit$ages,
                          ay=obj$model.fit$cohorts[length(ax):(2*length(ax)-1)],
                          cy=obj$model.fit$years,
                          ax)
    p1 <- ggplot2::ggplot(data=df.fitted,
                          ggplot2::aes(x=dy,y=ax))+
      ggplot2::geom_line()+
      ggplot2::theme_classic()+
      ggplot2::xlab('Development year')+
      ggplot2::ylab(expression(alpha[x]))
  }else{p1<-ggplot2::ggplot()
        df.fitted= data.frame(dy=obj$model.fit$ages,
                              ay=obj$model.fit$cohorts[length(obj$model.fit$ages):(2*length(obj$model.fit$ages)-1)],
                              cy=obj$model.fit$years)}
  
  if(!is.null(obj$model.fit$kt[1,])){
    kt=as.vector(obj$model.fit$kt[1,])
    kt.f=as.vector(obj$model.fcst$kt.f$mean[1,])
    kt.u80=as.vector(obj$model.fcst$kt.f$upper[1,,1])
    kt.u95=as.vector(obj$model.fcst$kt.f$upper[1,,2])
    kt.l80=as.vector(obj$model.fcst$kt.f$lower[1,,1])
    kt.l95=as.vector(obj$model.fcst$kt.f$lower[1,,2])
    
    df.fitted= cbind(df.fitted,kt)
    df.forecasted=data.frame(ay=c(obj$model.fit$cohorts[length(kt):(2*length(kt)-1)],
                                  obj$model.fcst$gc.f$cohorts),
                             cy=c(obj$model.fit$years,
                                  obj$model.fcst$kt.f$years),
                             kt=c(kt,kt.f),
                             kt.u80=c(rep(0,length(kt)),kt.u80),
                             kt.u95=c(rep(0,length(kt)),kt.u95),
                             kt.l80=c(rep(0,length(kt)),kt.l80),
                             kt.l95=c(rep(0,length(kt)),kt.l95))
    
    p2 <- ggplot2::ggplot(data=df.fitted,
                          ggplot2::aes(x=cy,
                                       y=kt))+
      ggplot2::geom_line()+
      ggplot2::theme_classic()+
      ggplot2::xlab('Calendar year')+
      ggplot2::ylab(expression(kappa[t]))
    
    p2.f <- ggplot2::ggplot(data=df.forecasted,ggplot2::aes(x=cy,y=kt))+
      ggplot2::geom_line()+
      ggplot2::geom_ribbon(aes(ymin = kt.l80, ymax = kt.u80), alpha = 0.2)+
      ggplot2::geom_ribbon(aes(ymin = kt.l95, ymax = kt.u95), alpha = 0.1)+
      ggplot2::theme_classic()+
      ggplot2::xlab('Calendar year')+
      ggplot2::ylab(expression(kappa[t]))
    
  }else{p2<-ggplot2::ggplot()
      p2.f<-ggplot2::ggplot()}
  
  if(!is.null(obj$model.fit$kt[1,])&&dim(obj$model.fit$kt)[1]>1){
    kt2=as.vector(obj$model.fit$kt[2,])
    kt2.f=as.vector(obj$model.fcst$kt.f$mean[2,])
    kt2.u80=as.vector(obj$model.fcst$kt.f$upper[2,,1])
    kt2.u95=as.vector(obj$model.fcst$kt.f$upper[2,,2])
    kt2.l80=as.vector(obj$model.fcst$kt.f$lower[2,,1])
    kt2.l95=as.vector(obj$model.fcst$kt.f$lower[2,,2])
    
    df.fitted= cbind(df.fitted,kt2)
    df.forecasted=cbind(df.forecasted,
                        data.frame(kt2=c(kt2,kt2.f),
                                   kt2.u80=c(rep(0,length(kt2)),kt2.u80),
                                   kt2.u95=c(rep(0,length(kt2)),kt2.u95),
                                   kt2.l80=c(rep(0,length(kt2)),kt2.l80),
                                   kt2.l95=c(rep(0,length(kt2)),kt2.l95)))
    
    p22 <- ggplot2::ggplot(data=df.fitted,
                          ggplot2::aes(x=cy,
                                       y=kt2))+
      ggplot2::geom_line()+
      ggplot2::theme_classic()+
      ggplot2::xlab('Calendar year')+
      ggplot2::ylab(expression(kappa[t]^2))
    
    p22.f <- ggplot2::ggplot(data=df.forecasted,ggplot2::aes(x=cy,y=kt2))+
      ggplot2::geom_line()+
      ggplot2::geom_ribbon(aes(ymin = kt2.l80, ymax = kt2.u80), alpha = 0.2)+
      ggplot2::geom_ribbon(aes(ymin = kt2.l95, ymax = kt2.u95), alpha = 0.1)+
      ggplot2::theme_classic()+
      ggplot2::xlab('Calendar year')+
      ggplot2::ylab(expression(kappa[t]^2))
    
  }else{p22<-ggplot2::ggplot()
  p22.f<-ggplot2::ggplot()}
  
  if(!is.null(obj$model.fit$kt[1,])&&dim(obj$model.fit$kt)[1]>2){
    kt3=as.vector(obj$model.fit$kt[3,])
    kt3.f=as.vector(obj$model.fcst$kt.f$mean[3,])
    kt3.u80=as.vector(obj$model.fcst$kt.f$upper[3,,1])
    kt3.u95=as.vector(obj$model.fcst$kt.f$upper[3,,2])
    kt3.l80=as.vector(obj$model.fcst$kt.f$lower[3,,1])
    kt3.l95=as.vector(obj$model.fcst$kt.f$lower[3,,2])
    
    df.fitted= cbind(df.fitted,kt3)
    df.forecasted=cbind(df.forecasted,
                        data.frame(kt3=c(kt3,kt3.f),
                                   kt3.u80=c(rep(0,length(kt3)),kt3.u80),
                                   kt3.u95=c(rep(0,length(kt3)),kt3.u95),
                                   kt3.l80=c(rep(0,length(kt3)),kt3.l80),
                                   kt3.l95=c(rep(0,length(kt3)),kt3.l95)))
    
    p23 <- ggplot2::ggplot(data=df.fitted,
                           ggplot2::aes(x=cy,
                                        y=kt3))+
      ggplot2::geom_line()+
      ggplot2::theme_classic()+
      ggplot2::xlab('Calendar year')+
      ggplot2::ylab(expression(kappa[t]^3))
    
    p23.f <- ggplot2::ggplot(data=df.forecasted,ggplot2::aes(x=cy,y=kt3))+
      ggplot2::geom_line()+
      ggplot2::geom_ribbon(aes(ymin = kt3.l80, ymax = kt3.u80), alpha = 0.2)+
      ggplot2::geom_ribbon(aes(ymin = kt3.l95, ymax = kt3.u95), alpha = 0.1)+
      ggplot2::theme_classic()+
      ggplot2::xlab('Calendar year')+
      ggplot2::ylab(expression(kappa[t]^3))
    
  }else{p23<-ggplot2::ggplot()
  p23.f<-ggplot2::ggplot()}
  
  if(!is.null(obj$model.fit$bx[,1])){
    bx=as.vector(obj$model.fit$bx[,1])
    df.fitted=cbind(df.fitted,bx)
    
    p3 <- ggplot2::ggplot(data=df.fitted,
                          ggplot2::aes(x=dy,
                                       y=bx))+
      ggplot2::geom_line()+
      ggplot2::theme_classic()+
      ggplot2::xlab('Development year')+
      ggplot2::ylab(expression(beta[x]^1))
    
  }else{p3<-ggplot2::ggplot()}
  
  if(!is.null(obj$model.fit$bx)&&dim(obj$model.fit$bx)[2]>1){
    bx2=as.vector(obj$model.fit$bx[,2])
    df.fitted=cbind(df.fitted,bx2)
    
    p32 <- ggplot2::ggplot(data=df.fitted,
                          ggplot2::aes(x=dy,
                                       y=bx2))+
      ggplot2::geom_line()+
      ggplot2::theme_classic()+
      ggplot2::xlab('Development year')+
      ggplot2::ylab(expression(beta[x]^2))
    
  }else{p32<-ggplot2::ggplot()}
  
  if(!is.null(obj$model.fit$bx)&&dim(obj$model.fit$bx)[2]>2){
    bx3=as.vector(obj$model.fit$bx[,3])
    df.fitted=cbind(df.fitted,bx3)
    
    p33 <- ggplot2::ggplot(data=df.fitted,
                           ggplot2::aes(x=dy,
                                        y=bx3))+
      ggplot2::geom_line()+
      ggplot2::theme_classic()+
      ggplot2::xlab('Development year')+
      ggplot2::ylab(expression(beta[x]^3))
    
  }else{p33<-ggplot2::ggplot()}
  
  if(!is.null(obj$model.fit$gc)){
    gc=as.vector(obj$model.fit$gc[!is.na(obj$model.fit$gc)])
    gc.f=as.vector(obj$model.fcst$gc.f$mean)
    gc.u80=obj$model.fcst$gc.f$upper[1:length(gc)]
    gc.u95=obj$model.fcst$gc.f$upper[(length(gc)+1):(2*length(gc))]
    gc.l80=obj$model.fcst$gc.f$lower[1:length(gc)]
    gc.l95=obj$model.fcst$gc.f$lower[(length(gc)+1):(2*length(gc))]
    
    df.fitted= cbind(df.fitted,gc)
    df.forecasted=data.frame(df.forecasted,
                             gc=c(gc,gc.f),
                             gc.u80=c(rep(0,length(gc)),gc.u80),
                             gc.u95=c(rep(0,length(gc)),gc.u95),
                             gc.l80=c(rep(0,length(gc)),gc.l80),
                             gc.l95=c(rep(0,length(gc)),gc.l95))
    
    p4 <- ggplot2::ggplot(data=df.fitted,
                          ggplot2::aes(x=ay,
                                       y=gc))+
      ggplot2::geom_line()+
      ggplot2::theme_classic()+
      ggplot2::xlab('Calendar year')+
      ggplot2::ylab(expression(gamma[t-x]))
    
    p4.f <- ggplot2::ggplot(data=df.forecasted,
                            ggplot2::aes(x=ay,y=gc))+
      ggplot2::geom_line()+
      ggplot2::geom_ribbon(ggplot2::aes(ymin = gc.l80, ymax = gc.u80), alpha = 0.2)+
      ggplot2::geom_ribbon(ggplot2::aes(ymin = gc.l95, ymax = gc.u95), alpha = 0.1)+
      ggplot2::theme_classic()+
      ggplot2::xlab('Calendar year')+
      ggplot2::ylab(expression(gamma[t-x]))
  }else{p4<-ggplot2::ggplot()
        p4.f<-ggplot2::ggplot()}
  
  if(!is.null(obj$model.fit$b0x)){
    b0x=as.vector(obj$model.fit$b0x)
    df.fitted=cbind(df.fitted,b0x)
    p5 <- ggplot2::ggplot(data=df.fitted,
                          ggplot2::aes(x=ay,y=b0x))+
      ggplot2::geom_line()+
      ggplot2::theme_classic()+
      ggplot2::xlab('Development year')+
      ggplot2::ylab(expression(beta[x]^0))
  }else{p5<-ggplot2::ggplot()}
  
  lay <- rbind(c(1,NA),
               c(2,3),
               c(4,5))
  
  plt<-gridExtra::grid.arrange(p1,
                          p4,
                          p5,
                          p2,
                          p3,
                          layout_matrix = lay,
                          top = grid::textGrob("Fitted effects",gp=grid::gpar(fontsize=20)))
  
  if(!is.null(obj$model.fit$bx[,2])||!is.null(obj$model.fit$bx[,3])){
  lay.extra <- rbind(c(1,2),
                c(3,4))
  plt<-gridExtra::grid.arrange(p22,
                               p32,
                               p23,
                               p33,
                               layout_matrix = lay.extra,
                               top = grid::textGrob("Fitted effects",gp=grid::gpar(fontsize=20)))}
  
  
  lay2 <- rbind(c(1),
               c(2))
  
  plt2<-gridExtra::grid.arrange(p2.f,p4.f,layout_matrix = lay2,
                          top = grid::textGrob("Forecasted effects",gp=grid::gpar(fontsize=20)))
  
  if(!is.null(obj$model.fit$bx[,2])||!is.null(obj$model.fit$bx[,3])){
    
    plt2.extra<-gridExtra::grid.arrange(p22.f,p23.f,layout_matrix = lay2,
                                  top = grid::textGrob("Forecasted effects",gp=grid::gpar(fontsize=20)))
    
    
  }
}






