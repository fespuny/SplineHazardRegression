#' Objective function for b-spline hazard regression
#'
#' @description Calculates the objective function for b-spline hazard regression using as input the needed pre-calculated matrices
#'
#' @param par.alpha - coefficients of a b-spline in a b-pline basis B (not provided)
#' @param yd - time-to-event data observations matrix, with each row being a pair (time,status)
#' @param Wik - matrix that multiplied by par.alpha gives the hazard function
#' @param Zik - matrix that multiplied by par.alpha gives the cumulative hazard function
#'
#' @return l - objective function (-2 *log likelihood function)
#' @export
#'
srllikb_fun = function(par.alpha, yd, Wik, Zik){

  alphaT = as.numeric(par.alpha) #arrays work as column vectors in matrix multiplication

  aux = Wik%*%alphaT

  logWikalpha = log( aux )

  if( sum( aux<=0 ) > 0 ){
    logWikalpha = Inf*matrix( 0, nrow(aux), ncol(aux) )
  }

  l = -2*sum( yd[, 2] * logWikalpha - Zik%*%alphaT )

  if( is.na(l)) l = Inf

  ## L-BFGS-B needs finite values
  if( l == Inf ) l = 0.99 * .Machine$double.xmax
  if( l == -Inf ) l = - 0.99 * .Machine$double.xmax

  return( l )
}
