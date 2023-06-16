#' Generate hazard and censoring distributions
#'
#' @description Generates hazard and censoring distributions and respective survival functions [IS IT OK TO IGNORE LATE ENTRY AT THIS STAGE?]
#'
#' @param Hazard - hazard function for outcome variable.
#'        'exp' for exponential, 'weib' for Weibull, 'pe' for piecewise-exponential, or 'spline' for B-Spline
#' @param HParam - matrix of parameter values for outcome variable \cr
#' @param Censor - censoring hazard ('', 'exp', 'weib', 'pe' or 'bspline')
#' @param CParam - matrix of parameter values for censoring variable
#' @param Tmax - scalar, maximum follow-up time, default 10
#' @param SampleSize - scalar, sample size
#'
#' @details For 'exp', HParm contains a scalar;
#'     For 'pe', HParm(:, 1) is left limit, HParm(:, 2) is right limit, and HParm(:, 3) is hazard over the interval
#'     For 'bspline', HParm(:, 1) lists knots, HParm(:, 2) lists spline coefficients.
#'     The number of rows of Hparm has to be equal to the number of degrees of freedom: number of interior knots plus order (=degree plus one) of the b-spline
#'
#' @return Returns a list of outputs (t, basis, h, hCensor, Sh, Scensor)
#'
#'     t = array of times
#'     basis = basis of b-spline cubic functions [WE NEED TO ADD ORDER/DEGREE AS PARAMETER IF WE WANT TO ALLOW DIFFERENT B-SPLINES]
#'     h = simulated hazard distribution
#'     hCensor = simulated censoring distribution
#'     Sh = cumulative hazard survival function
#'     SCensor = cumulative censoring survival funtion
#'
#' @export
#'
#' @examples
#'     #Generate input for a b-spline hazard simulation with piecewise exponential survival censoring
#'     knots = c(0, 1, 3, 6, 10, NaN, NaN )
#'     betac = 1 * c(0.05, 0.05, 0.05, 0.05, 0.40, 0.1, 0.05)
#'     HParm = data.frame(knots, betac) # 'A Simple B-Spline'
#'     cll = c(0, 5)
#'     cup = c(5, 10)
#'     cih = c(0.0125, 0.025)
#'     CParm = data.frame(cll, cup, cih) # 'Light Censoring'
#'     INPUTS = etsim_inputs( HParam=HParm, CParam=CParm )
#'
etsim_inputs <- function( Hazard="spline", HParam, Censor="pe", CParam, Tmax=10, SampleSize=101 ){

  ##we ignore duplicated knots and assign boundary and interior knot variables
  knots = unique( setdiff( as.numeric(HParam[,1]), c(NA, NaN) ) )
  Boundary.knots = c(min(knots),max(knots))
  Interior.knots = setdiff( knots, Boundary.knots )

  betac = as.numeric(HParam[,2])
  cll = as.numeric(CParam[,1])
  cup = as.numeric(CParam[,2])
  cih = as.numeric(CParam[,3])

  t = seq(from=Boundary.knots[1],to=Boundary.knots[2],length.out=SampleSize)

  B = generate_bspline_basis( time=t, Interior.knots=Interior.knots, Boundary.knots=Boundary.knots, ORDER=4 )

  h = B %*% betac

  # iB = ibs(x=t, knots=Interior.knots, Boundary.knots=Boundary.knots, degree=3, intercept=TRUE)
  # Sh = exp( - iB %*% betac )

  DELTA = (Boundary.knots[2]-Boundary.knots[1])/(SampleSize-1)

  Sh = exp( - cumsum(h*DELTA) )

  hCensor = rep( 0, SampleSize )
  for( i in 1:length(cll) ){
    hCensor[ t > cll[i] & t<= cup[i] ] = cih[i]
  }

  Scensor = exp( - cumsum( hCensor*DELTA) )

  return(
    list(t=t,
         basis=B,
         b=betac,
         h = h,
         hCensor=hCensor,
         Sh = Sh,
         Scensor=Scensor)
  )
}
