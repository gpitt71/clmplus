#' Pre-process Run-Off Triangles
#'
#' Pre-process Run-Off Triangles.
#' @param cumulative.payments.triangle A square numeric matrix with at least
#'   two rows. Rows are accident periods, columns are development periods, and
#'   observed upper-triangle cells satisfy `row + column <= J + 1`. Values are
#'   non-negative cumulative paid amounts in the source data's monetary units,
#'   non-decreasing across each row; unavailable cells may be `NA`. Recoveries
#'   (negative incremental payments) are not supported.
#' @param entries.weights Optional non-negative numeric `J` by `J` matrix of
#'   fitting weights in the same accident/development layout. `NULL` gives
#'   observed cells weight one. The first development period and missing cells
#'   are always zero-weighted after conversion to calendar coordinates.
#' @param eta One finite numeric value in `(0, 1]`, default `0.5`, describing
#'   expected within-cell payment timing (lost exposure). It is used to derive
#'   exposure and to convert fitted hazards to development factors.
#' 
#' @examples
#' data(sifa.mtpl)
#' sifa.mtpl.rtt <- AggregateDataPP(cumulative.payments.triangle=sifa.mtpl)
#' 
#' @return An `AggregateDataPP` list with:
#'   \item{cumulative.payments.triangle}{The input `J` by `J` cumulative paid
#'   triangle, unchanged.}
#'   
#'   \item{occurrance}{A `J` by `J` matrix of incremental paid amounts in
#'   development-period by calendar-period coordinates. The misspelling is
#'   retained as a stable public field name.}
#'   
#'   \item{exposure}{A `J` by `J` numeric matrix in the same calendar
#'   coordinates, calculated as cumulative payments minus
#'   `(1 - eta) * occurrence`.}
#'   
#'   \item{incremental.payments.triangle}{A `J` by `J`
#'   accident/development matrix of incremental paid amounts.}
#'   
#'   \item{fit.w}{The `J` by `J` fitting-weight matrix in
#'   development/calendar coordinates.}
#'   
#'   \item{J}{The integer triangle dimension.}
#'   
#'   \item{diagonal}{A length-`J` numeric vector containing the latest observed
#'   cumulative diagonal in calendar representation.}
#'   
#'   \item{eta}{The supplied within-cell timing scalar.}
#'  
#'   
#' @references 
#' Pittarello, G., Hiabu, M., & Villegas, A. M. (2023). Replicating and extending chain-ladder via an age-period-cohort structure on the claim development in a run-off triangle. arXiv preprint arXiv:2301.03858.
#' 
#' @export
AggregateDataPP <- function(cumulative.payments.triangle, entries.weights=NULL, eta=1/2)
{
  incrementals <- validate_triangle(cumulative.payments.triangle, entries.weights, eta)
  J <- ncol(cumulative.payments.triangle)
  
  # find out occurrance and exposure
  occurrance <- triangle_to_calendar(incrementals)
  calendar_cumulative <- triangle_to_calendar(cumulative.payments.triangle)
  exposure <- calendar_cumulative - (1 - eta) * occurrance
  
  
  # find out the weights
  if(is.null(entries.weights)){
    
    fit.w <- matrix(1,nrow=J,ncol = J) 

    
    }
  else{
    
    fit.w <- entries.weights
    
  }
  
  fit.w[,1]=0
  fit.w=triangle_to_calendar(fit.w)
  fit.w[is.na(fit.w)]=0
  
  
  tr <- list(
    cumulative.payments.triangle = cumulative.payments.triangle,
    occurrance = occurrance,
    exposure = exposure,
    fit.w=fit.w,
    incremental.payments.triangle = incrementals,
    J=J,
    diagonal=calendar_cumulative[,J],
    eta=eta
  )
  
  ## Set the name for the class
  class(tr) <- "AggregateDataPP"
  tr
  
}


