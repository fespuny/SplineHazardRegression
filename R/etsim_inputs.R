#' Generate hazard and censoring distributions
#'
#' @description Generates hazard and censoring distributions and respective survival functions [IS IT OK TO IGNORE LATE ENTRY AT THIS STAGE?]
#'
#' @param Hazard - hazard function for outcome variable.
#'        'exp' for exponential, 'weib' for Weibull, 'pe' for piecewise-exponential, or 'spline' for B-Spline
#' @param HParam - matrix of parameter values for outcome variable \cr
#' @param Censor - censoring hazard ('', 'exp', 'weib', 'pe' or 'spline')
#' @param CParam - matrix of parameter values for censoring variable
#' @param Tmin - scalar, minimum follow-up time, default 0
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
etsim_inputs <- function( Hazard="spline", HParam, Censor="pe", CParam=NULL, Tmin=0, Tmax=10, SampleSize=101 ){

  if( sum( Hazard %in% c("exp","spline","pe" ), na.rm=TRUE ) ){

    Hazard_etsim = etsim_inputs_core( Hazard, HParam, Tmin, Tmax, SampleSize )

  } else {
    stop("Argmument 'Hazard' can only take values 'exp', 'spline', or 'pe'. See help(etsim_inputs) for details.")
  }

  # Missing HParam throws an error
  # Wrong HParam makes program crash

  #Note: censoring definition requires that time field "t" has been already defined
  if( nchar(Censor)==0 | sum( Censor %in% c("exp","spline","pe" ), na.rm=TRUE ) ){

    Censor_etsim = etsim_inputs_core( Censor, CParam, Tmin, Tmax, SampleSize )

  } else {
    stop("Argmument 'Censor' can only take values '', 'exp', 'spline', or 'pe'. See help(etsim_inputs) for details.")
  }

  warning_message = "The provided knots for the Spline hazard do not include Tmin and Tmax"
  if( Hazard=="spline" & ! (sum(Tmin %in% HParam[,1],na.rm=T)>0 & sum(Tmax %in% HParam[,1],na.rm=T)>0 ) ){
    warning( warning_message, immediate. = TRUE )
  }
  if( Censor=="spline" & ! (sum(Tmin %in% CParam[,1],na.rm=T)>0 & sum(Tmax %in% CParam[,1],na.rm=T)>0 ) ){
    warning( warning_message, immediate. = TRUE )
  }

  return(
    list(
      t=Hazard_etsim$t,
      h = Hazard_etsim$h,
      hCensor=Censor_etsim$h,
      Sh = Hazard_etsim$Sh,
      Scensor=Censor_etsim$Sh
    )
  )

}
