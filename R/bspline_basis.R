#' Generate B-Spline basis
#'
#' @description Wraper function to a B-spline basis calculation function in R (currently from package bspline2)
#'
#' @param time - x coordinates for the B-spline functions
#' @param Interior.knots - interior knots of the B-spline basis
#' @param Boundary.knots - boundary knots of the B-spline basis
#' @param ORDER - 1 step, 2 linear, 3 quadratic, 4 cubic
#' @param integral - TRUE only if you want the integrals of a B-spline basis
#'
#' @details The current bSpline() function from the splines2 package extends the bs() function in the splines package for B-spline basis by allowing piecewise constant (left-closed and right-open except on the right boundary) spline basis of degree zero.
#'
#' @return B - matrix of basis of b-pline functions of the specified degree, evaluated at the "time" points
#'         B has dimensions length(time) x (number of interior points + degree + 1)
#' @importFrom splines2 bSpline
#' @importFrom graphics matplot
#' @export
#'
#' @details See the function splines2::bspline for further details
#'
#' @examples # A basis of cubic B-splines with no interior points
#' B = generate_bspline_basis( time = 0:10, Interior.knots=c(), Boundary.knots=c(0,10) )
#' B
#' # Note that B has 4=3+1 functions (as many as the order of cubic B-splines)
#' # we can plot the basis using matplot
#' matplot( 0:10, B, type="l")
#'
generate_bspline_basis <- function( time, Interior.knots, Boundary.knots, ORDER=4, integral=FALSE ){

  B = bSpline( x=time,
               knots=Interior.knots,
               Boundary.knots=Boundary.knots,
               degree=ORDER-1,
               intercept=TRUE,
               integral = integral)
  return( B )
}
