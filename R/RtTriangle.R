#' Reverse time triangles
#'
#' This function allows to define the class of triangles for reverse time models.
#' @param cumulative.payments.triangle \code{triangle matrix} or \code{matrix array} object, input triangle of cumulative payments.
#' @param k \code{numeric}, individual claims exposure in the cell, also known as lost exposure. It must be in the interval (0,1].
#' 
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- RtTriangle(cumulative.payments.triangle=sifa.mtpl)
#' 
#' @return An object of class \code{RtTriangle}. Lists the following elements:
#'   \item{cumulative.payments.triangle}{\code{triangle matrix} object, input triangle of cumulative payments.}
#'   
#'   \item{occurrance}{\code{matrix array} object, the occurrance derived from the input triangle.}
#'   
#'   \item{exposure}{\code{matrix array} object, the exposure derived from the input triangle, under the \code{k} claims arrival assumption.}
#'   
#'   \item{incremental.payments.triangle}{\code{triangle matrix} object, incremental payments derived from the input.}
#'   
#'   \item{J}{\code{integer}, Run-off triangle dimension.}
#'   
#'   \item{diagonal}{ \code{numeric}, cumulative payments last diagonal.}
#'  
#'   
#' @references 
#' Hiabu, Munir. “On the relationship between classical chain ladder and granular reserving.” 
#' Scandinavian Actuarial Journal 2017 (2017): 708 - 729.
#' 
#' @export
RtTriangle <- function(cumulative.payments.triangle, k=1/2)
{
  
  rtt.input.env$properties.cpt(cumulative.payments.triangle)
  
  incrementals = ChainLadder::cum2incr(cumulative.payments.triangle)
  J=dim(cumulative.payments.triangle)[2]
  
  # find out occurrance and exposure
  occurrance=pkg.env$t2c(incrementals)
  drop=1-k
  exposure=pkg.env$t2c(cumulative.payments.triangle-drop*incrementals)
  
  # occurrance[is.na(occurrance)]=c(0.)
  # exposure[is.na(occurrance)]=c(0.)
  
  # find out the weights
  fit.w <- matrix(1,nrow=J,ncol = J) 
  fit.w[,1]=0
  fit.w=pkg.env$t2c(fit.w)
  fit.w[is.na(fit.w)]=0
  
  tr <- list(
    cumulative.payments.triangle = cumulative.payments.triangle,
    occurrance = occurrance,
    exposure = exposure,
    fit.w=fit.w,
    incremental.payments.triangle = incrementals,
    J=J,
    diagonal=pkg.env$t2c(cumulative.payments.triangle)[,J],
    k=k
  )
  
  ## Set the name for the class
  class(tr) <- "RtTriangle"
  tr
  
}


