#' Core function used to simulate hazard and censoring distributions
#'
#' @description Generates hazard and censoring distributions and respective survival functions [IS IT OK TO IGNORE LATE ENTRY AT THIS STAGE?]
#'
#' @param Hazard_type - hazard function for outcome variable.
#'        'exp' for exponential, 'weib' for Weibull, 'pe' for piecewise-exponential, or 'spline' for B-Spline
#' @param Param - matrix of parameter values for outcome variable
#' @param Tmin - scalar, minimum follow-up time, default 0
#' @param Tmax - scalar, maximum follow-up time, default 10
#' @param SampleSize - scalar, sample size, default 101
#'
#' @details For 'exp', Parm contains a scalar;
#'     For 'pe', Parm(:, 1) is left limit, Parm(:, 2) is right limit, and Parm(:, 3) is hazard over the interval
#'     For 'bspline', Parm(:, 1) lists knots, Parm(:, 2) lists b-spline coefficients.
#'     The number of rows of Parm has to be equal to the number of degrees of freedom: number of interior knots plus order (=degree plus one) of the b-spline
#'
#' @return Returns a list of outputs (t, h, Sh )
#'     t = array of times
#'     h = simulated hazard distribution
#'     Sh = cumulative hazard survival function
#'
#' @export
#'
#' @examples
#'     #PENDING
#'     #Generate input for a cubic b-spline hazard simulation
#'     knots = c(0, 1, 3, 6, 10, NaN, NaN )
#'     betac = 1 * c(0.05, 0.05, 0.05, 0.05, 0.40, 0.1, 0.05)
#'     HParm = data.frame(knots, betac) # 'A Simple B-Spline'
#'     INPUTS_hazard = etsim_inputs_core( Param=HParm )
#'
etsim_inputs_core <- function( Hazard_type="spline", Param, Tmin=0, Tmax=10, SampleSize=101 ){

  t = seq(from=Tmin,to=Tmax,length.out=SampleSize)

  DELTA = (Tmax-Tmin)/(SampleSize-1)

  if( nchar(Hazard_type)==0 ){
    h = rep( 0, SampleSize )
  }

  ## EXPONENTIAL SURVIVAL
  if( Hazard_type=="exp" ){
    ih = as.numeric(Param)
    h = rep( ih, SampleSize )
  }

  ## PIECEWISE EXPONENTIAL SURVIVAL
  if( Hazard_type=="pe" ){
    ll = as.numeric(Param[,1])
    up = as.numeric(Param[,2])
    ih = as.numeric(Param[,3])

    h = rep( 0, SampleSize )
    for( i in 1:length(ll) ){
      h[ t > ll[i] & t<= up[i] ] = ih[i]
    }
  }

  ## BSPLINE HAZARD
  if( Hazard_type == "spline" ){

    ##we ignore duplicated knots and assign boundary and interior knot variables
    knots = unique( setdiff( as.numeric(Param[,1]), c(NA, NaN) ) )
    Boundary.knots = c(min(knots),max(knots))
    Interior.knots = setdiff( knots, Boundary.knots )

    betac = as.numeric(Param[,2])

    ORDER = length(betac) - length(Interior.knots)

    B = generate_bspline_basis( time=t, Interior.knots=Interior.knots, Boundary.knots=Boundary.knots, ORDER=ORDER )

    h = B %*% betac
  }

  Sh = exp( - cumsum(h*DELTA) )

  return(
    list(t=t,
         h = h,
         Sh = Sh)
  )
}
