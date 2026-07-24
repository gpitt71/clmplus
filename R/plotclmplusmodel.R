#' Plot the hazard model residuals 
#'
#' This function allows to plot the hazard model residuals on the triangle payments.
#' 
#' @param x A fitted `clmplusmodel` object.
#' @param heat.lim A length-two numeric vector giving the lower and upper fill
#'   scale limits for scaled deviance residuals.
#' @param ... Reserved; currently ignored.
#' 
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- AggregateDataPP(cumulative.payments.triangle=sifa.mtpl)
#' clm.fit<-clmplus(sifa.mtpl.rtt, 'a')
#' plot(clm.fit)
#' 
#' @return A `ggplot` object showing scaled deviance residuals in
#'   accident-year by development-year triangle form.
#' 
#' @references 
#' 
#' Pittarello, Gabriele, Munir Hiabu, and Andrés M. Villegas. "Replicating and extending chain ladder 
#' via an age-period-cohort structure on the claim development in a run-off triangle." arXiv preprint arXiv:2301.03858 (2023).
#' 
#' @export
plot.clmplusmodel <- function(x,
                              heat.lim=c(-2.5,2.5),
                              ...){
  
  if (x$apc_input$hazard.model %in% names(supported_models)) {
    res.tr <- x$hazard_scaled_deviance_residuals
    colnames(res.tr) <- rownames(res.tr) <- c(0:(dim(res.tr)[2]-1))
    longdf.no.0 = ChainLadder::as.LongTriangle(res.tr)
    
  }
  
  p_hm <- ggplot2::ggplot(data=longdf.no.0, ggplot2::aes(x=as.integer(dev)-1, y=as.integer(origin)-1)) + 
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
    ggplot2::ggtitle(x$modelfamily)+
    ggplot2::theme(axis.title.x = ggplot2::element_text(size=8), axis.text.x  = ggplot2::element_text(size=7))+
    ggplot2::theme(axis.title.y = ggplot2::element_text(size=8), axis.text.y  = ggplot2::element_text(size=7))+
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "grey", colour = "grey", linewidth = 2, linetype = "solid"),
          panel.grid = ggplot2::element_line(colour="grey")) + 
    NULL
  
    p_hm
  
}








