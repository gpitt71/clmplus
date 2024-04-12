#' Predict the Reserve using Chain Ladder Plus Models
#'
#' This function allows to define the behavior of the triangle payments.
#' 
#' @param model clmplus model to be plotted. 
#' @param gk.fc.model \code{character}, model to forecast the cohort component for the last accident period. It can be either arima ('a') or linear model ('l'). Disregarded for models that do not have a cohort effect.
#' @param ckj.fc.model \code{character}, model to forecast the calendar period effect. It can be either arima ('a') or linear model ('l'). Disregarded for models that do not have a period effect.
#' @param gk.order \code{integer}, order of the arima model with drift for the accident year effect extrapolation. Default to (1,1,0).
#' @param ckj.order \code{integer}, order of the arima model with drift for the calendar year effect extrapolation. Default to (0,1,0).
#' @param forecasting_horizon \code{integer}, between 1 and the triangle width. Calendar periods ahead for the predictions. Default predictions are to run-off. 
#' 
#' 
#' @return predictions.
#' 
#' @references 
#' Pittarello, Gabriele, Munir Hiabu, and Andrés M. Villegas. "Replicating and extending chain ladder 
#' via an age-period-cohort structure on the claim development in a run-off triangle." arXiv preprint arXiv:2301.03858 (2023).
#'  
#' @export
predict.clmplusmodel <- function(clmplusmodel,
                                 gk.fc.model='a',
                                 ckj.fc.model='a',
                                 gk.order=c(1,1,0),
                                 ckj.order=c(0,1,0),
                                 forecasting_horizon=NULL,
                                 ...){
  
  
  # forecasting horizon
  J <- clmplusmodel$apc_input$J
  # fitted model
  model <- clmplusmodel$model.fit
  # occurrences distribution
  eta <- clmplusmodel$apc_input$eta
  # hazard model
  hazard.model <- clmplusmodel$apc_input$hazard.model
  # diagonal
  d <- clmplusmodel$apc_input$diagonal[1:J-1]
  # cumulative payments 
  cumulative.payments.triangle <- clmplusmodel$apc_input$cumulative.payments.triangle
  
  if(hazard.model %in% names(pkg.env$models)){
    
  alphaij <- pkg.env$fcst(model, 
                          hazard.model = hazard.model,
                          gk.fc.model=gk.fc.model,
                          ckj.fc.model=ckj.fc.model,
                          gk.order=gk.order,
                          ckj.order=ckj.order
  )
  
  fij=(1+(1-eta)*alphaij$rates)/(1-(eta*alphaij$rates))
  
  # d=AggregateDataPP$diagonal[1:(J-1)]
  if(is.null(forecasting_horizon)){
    J.stop=J
  }else{
    J.stop=forecasting_horizon
    }
  
  
  lt=array(0.,c(J,J.stop))
  
  if(!is.null(forecasting_horizon)){
    # plus one is added artificially to make indexing consistent
    J.stop=J.stop+1
  }
  
  lt[,1]=c(0.,d)*fij[,1]
  
  if(J.stop>2){
    for(j in 2:J.stop){lt[,j]=c(0.,lt[1:(J-1),(j-1)])*fij[,j]} 
  }
  
  # ot_=pkg.env$t2c(AggregateDataPP$cumulative.payments.triangle)

  ot_=pkg.env$t2c(cumulative.payments.triangle)
  
  if(is.null(forecasting_horizon)){
    
    
    ultimate_cost=c(rev(lt[J,1:(J-1)]),ot_[J,J])
    reserve=rev(ultimate_cost-ot_[,J])
    ultimate_cost=rev(ultimate_cost)
    #remove last column (not necessary)
    lt <- lt[,1:(dim(lt)[2]-1)]
    fij <- fij[,1:(dim(fij)[2]-1)]
    
  }else{
      
      ultimate_cost=rev(lt[,J.stop-1])
      ultimate_cost<-c(ultimate_cost[ultimate_cost==0],ultimate_cost[ultimate_cost!=0])
      ultimate_cost[1]=c(ot_[J,J])
      
      fij <- fij[1:dim(lt)[1],1:dim(lt)[2],drop=FALSE]
      
      if(J.stop>2){
        
        for(cohort in 2:(J.stop-1)){
          
          ultimate_cost[cohort] <- lt[J,cohort-1]
          
        }
        
      }
      
      reserve=ultimate_cost-rev(ot_[,J])
      
    }
  
  
  out <- list(reserve=reserve,
              ultimate_cost=ultimate_cost,
              full_triangle= pkg.env$create_full_triangle(cumulative.payments.triangle=cumulative.payments.triangle,
                                                          lt),
              lower_triangle = pkg.env$create_lower_triangle(lt),
              development_factors_predicted = pkg.env$create_lower_triangle(fij),
              apc_output=list(alphaij=alphaij,
                              lower_triangle_apc=lt,
                              development_factors_apc=fij))
  
  return(out)
  
  
  
  }
  
  
  
  
  
  
  
  
}