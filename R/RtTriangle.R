#' Reverse time triangles
#'
#' This function allows to define the class of triangles for reverse time models.
#' @param cumulative.payments.triangle Input triangle of cumulative payments.
#' 
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- RtTriangle(cumulative.payments.triangle=sifa.mtpl)
#' 
#' @return An object of class \code{"RtTriangle"}. Lists the following elements:
#'   \item{cumulative.payments.triangle}{Input triangle of cumulative payments.}
#'   
#'   \item{occurrance}{Matrix that contains the occurrance derived from the input triangle.}
#'   
#'   \item{exposure}{Matrix that contains the exposure derived from the input triangle, under the uniform claims arrival assumption.}
#'   
#'   \item{incremental.payments.triangle}{Triangle of incremental payments derived from the input.}
#'   
#'   \item{J}{Run-off triangle dimension.}
#'   
#'   \item{diagonal}{Cumulatives payments last diagonal.}
#'  
#'   
#' @references 
#' Hiabu, Munir. “On the relationship between classical chain ladder and granular reserving.” 
#' Scandinavian Actuarial Journal 2017 (2017): 708 - 729.
#' 
#' @export
RtTriangle <- function(cumulative.payments.triangle)
{
  
  rtt.input.env$properties.cpt(cumulative.payments.triangle)
  
  incrementals = ChainLadder::cum2incr(cumulative.payments.triangle)
  J=dim(cumulative.payments.triangle)[2]
  
  # find out occurrance and exposure
  occurrance=pkg.env$t2c(incrementals)
  exposure=pkg.env$t2c(cumulative.payments.triangle-incrementals/2)
  
  # occurrance[is.na(occurrance)]=c(0.)
  # exposure[is.na(occurrance)]=c(0.)
  
  tr <- list(
    cumulative.payments.triangle = cumulative.payments.triangle,
    occurrance = occurrance,
    exposure = exposure,
    incremental.payments.triangle = incrementals,
    J=J,
    diagonal=pkg.env$t2c(cumulative.payments.triangle)[,J]
  )
  
  ## Set the name for the class
  class(tr) <- "RtTriangle"
  tr
  
}


