#' Generate hazard and censoring distributions
#'
#' @description Generates hazard and censoring distributions and respective survival functions [IS IT OK TO IGNORE LATE ENTRY AT THIS STAGE?]
#'
#' @param fullname string [NOT USED]
#' @param Hazard hazard function for outcome variable (REQUIRED)
#'        'exp' for exponential, 'weib' for Weibull, 'pe' for piecewise-exponential, or 'spline' for B-Spline
#' @param HParm matrix of parameter values for outcome variable
#' @param Censor censoring hazard, '', 'exp', 'weib', 'pe' or 'bspline'
#' @param Cparm matrix of parameter values for censoring variable
#' @param Tmax scalar, maximum follow-up time, default 10
#' @param Mesh scalar, discrete time-increment [NOT USED]
#' @param SampleSize scalar, sample size
#' @param B scalar, number of simulations used in bootstrapping
#' @param OutputFormat string, 'double', 'tbl', or 'data_set' [NOT USED]
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
#'     INPUTS = etsim_inputs( HParm=HParm, Cparm=Cparm)
#'
etsim_inputs = function( fullname="",
                         Hazard="spline",
                         HParm,
                         Censor="pe",
                         Cparm,
                         Tmax=10,
                         Mesh=7/365,
                         SampleSize=101,
                         B=10,
                         OutputFormat="data_set"){
  # fullname: string
  # Hazard: hazard function for outcome variable (REQUIRED),
  #  'exp' for exponential, 'weib' for Weibull, 'pe' for
  #  piecewise-exponential, or 'spline' for B-Spline
  # HParm: tbl, parameter values for outcome variable
  # Censor: censoring hazard, '', 'exp', 'weib', 'pe' or 'bspline'
  # CParm: tbl, parameter values for censoring variable
  #  TMax: scalar, maximum follow-up time
  #  Mesh: scalar, discrete time-increment
  # SampleSize: scalar, sample size
  #  B: scalar, number of simulations
  # OutputFormat: string, 'double', 'tbl', or 'data_set'
  #
  # * For 'exp', HParm contains a scalar
  # * For 'pe', HParm(:, 1) is left limit, HParm(:, 2) is right limit,
  # HParm(:, 3) is hazard over the interval
  # * For 'bspline', HParm(:, 1) lists knots, HParm(:, 2) lists spline coefficients.
  #
  # Example:
  # I1 = etsim_inputs('Hazard', 'exp', 'HParm', tbl(0.1, 'Flat Hazard 10%/y'));

  knots = unique( setdiff( as.numeric(HParm$knots), c(NA, NaN) ) )   ##we ignore duplicated knots
  Boundary.knots = c(min(knots),max(knots))   ##boundary knots
  Interior.knots = knots[2:(length(knots)-1)]

  betac = as.numeric(HParm$betac)[ 1:(length(knots)+2) ] #the tail of betac is ignored

  t = seq(from=Boundary.knots[1],to=Boundary.knots[2],length.out=SampleSize)

  B = bSpline( x=t, knots=Interior.knots, Boundary.knots=Boundary.knots, degree=3, intercept=TRUE)

  h = B %*% betac

  # iB = ibs(x=t, knots=Interior.knots, Boundary.knots=Boundary.knots, degree=3, intercept=TRUE)
  # Sh = exp( - iB %*% betac )

  DELTA = (Boundary.knots[2]-Boundary.knots[1])/(SampleSize-1)

  Sh = exp( - cumsum(h*DELTA) )

  hCensor = rep(0, SampleSize)
  for( i in 1:length(CParm$cll) ){
    hCensor[ t > CParm$cll[i] & t<= CParm$cup[i] ] = CParm$cih[i]
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
