#' Generate B-Spline basis
#'
#' @description Wraper function to a B-spline basis calculation function in R (currently from package bspline2)
#'
#' @param time - x coordinates for the B-spline functions
#' @param Interior.knots - interior knots of the B-spline basis
#' @param Boundary.knots - bounddary knots of the B-spline basis
#' @param degree - degree of the B-spline functions
#'
#' @details The current bSpline() function from the splines2 package extends the bs() function in the splines package for B-spline basis by allowing piecewise constant (left-closed and right-open except on the right boundary) spline basis of degree zero.
#'
#' @return B - matrix of basis of b-pline functions of the specified degree, evaluated at the "time" points
#'         B has dimensions length(time) x (number of interior points + degree + 1)
#' @importFrom splines2 bSpline
#' @importFrom graphics matplot
#' @export
#'
#' @examples # A basis of cubic B-splines with no interior points will have 4=3+1 functions (as many as the order of cubic B-splines)
#' B = generate_bspline_basis( time = 0:10, Interior.knots=c(), Boundary.knots=c(0,10) )
#' B
#' # we can plot the basis using matplot
#' matplot( 0:10, B, type="l")
#'
generate_bspline_basis <- function( time, Interior.knots, Boundary.knots, degree=3 ){

  B = bSpline( x=time,
               Interior.knots=Interior.knots,
               Boundary.knots=Boundary.knots,
               degree=degree, intercept=TRUE)

  return( B )
}
