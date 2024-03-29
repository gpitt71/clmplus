#' Plot the payments behavior
#'
#' This function allows to define the behavior of the triangle payments.
#' 
#' @param x RtTriangle to be plotted.
#' @param ... Arguments to be passed to plot.
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- RtTriangle(cumulative.payments.triangle=sifa.mtpl)
#' plot(sifa.mtpl.rtt)
#' 
#' @return No return value, plots the run-off triangle cumulative payments and incremental payments.
#' 
#' @references 
#' Pittarello, Gabriele, Munir Hiabu, and Andrés M. Villegas. "Replicating and extending chain ladder 
#' via an age-period-cohort structure on the claim development in a run-off triangle." arXiv preprint arXiv:2301.03858 (2023).
#'  
#' Hiabu, Munir. “On the relationship between classical chain ladder and granular reserving.” 
#' Scandinavian Actuarial Journal 2017 (2017): 708 - 729.
#' 
#' @export
plot.RtTriangle <- function(x, ...){
  
  temp=x$incremental.payments.triangle
  colnames(temp)<-seq(0,dim(temp)[1]-1)
  rownames(temp)<-seq(1,dim(temp)[1])
  temp.long <- reshape2::melt(temp,
                        value.name=c('incrementals'))
  temp.long=cbind(temp.long,reshape2::melt(x$cumulative.payments.triangle,
                                             value.name=c('cumulatives'))[,'cumulatives'])
  colnames(temp.long) <-c('ay','dy','incrementals','cumulatives')
  temp.long$ay = as.factor(temp.long$ay)
  
  p1 <- ggplot2::ggplot(data=temp.long,
                        ggplot2::aes(x=dy,
                                     y=incrementals,
                                     by=ay,
                                     colour=ay))+
    ggplot2::geom_line()+
    ggplot2::ggtitle("Incremental payments")+
    ggplot2::theme_classic()
  
  p2 <- ggplot2::ggplot(data=temp.long,
                        ggplot2::aes(x=dy,
                                     y=cumulatives,
                                     by=ay,
                                     colour=ay))+
    ggplot2::geom_line()+
    ggplot2::ggtitle("Cumulative payments")+
    ggplot2::theme_classic()
  
  gridExtra::grid.arrange(p1,
                          p2,
                          ncol=2)
  
}


