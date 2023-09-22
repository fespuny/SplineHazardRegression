#' Gradient of the objective function for b-spline hazard regression
#'
#' @description Calculates the gradient of the objective function for b-spline hazard regression using as input the needed pre-calculated matrices
#'
#' @param par.alpha - coefficients of a b-spline in a b-pline basis B (not provided)
#' @param yd - time-to-event data observations matrix, with each row being a pair (time,status)
#' @param Wik - matrix that multiplied by par.alpha gives the hazard function
#' @param Zik - matrix that multiplied by par.alpha gives the cumulative hazard function
#'
#' @return dl - gradient of the objective function
#' @export
#'
srllikb_grad = function(par.alpha, yd, Wik, Zik){

  alphaT = as.numeric(par.alpha) #arrays work as column vectors in matrix multiplication
  np = length( alphaT )

  h = Wik %*% alphaT
  t1 = (yd[,2]/h) %*% matrix(1, 1, np)
  t1 = Wik*t1

  score = colSums(t1 - Zik, na.rm = T )

  dl = as.numeric( -2*score )

  return( dl )
}
