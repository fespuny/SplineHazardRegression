#' Generate the basis functions needed for fitting the hazard B-spline model
#'
#' @description This function calculates all matrices needed for the estimation of a hazard function with censored time-to-event data using B-splines
#'
#' @param yd - matrix of time-to-event data; pneumonic y is event time, d is for delta, the status indicator
#' @param entry - late entry / left truncation times
#' @param ORDER - 1 step, 2 linear, 3 quadratic, 4 cubic
#' @param knots - sequence of knot locations
#' @param t - vector of evaluation times
#'
#' @return list(Wik, Zik, XH, Xh ) - see details for definitions
#'
#' @details When fitting a B-spline function \eqn{h(\alpha)(t) := B(t) \alpha} to time-to-event-data,
#' where \eqn{B(t)} is a basis of B-splines and \eqn{\alpha} are the coefficients to be estimated,
#' we optimise a function of the likelihood \deqn{L(h(\theta)) =  }.
#'
#' @importFrom pracma zeros ones tril kron
#' @export
#'
bspline_regression_basis_functions <- function(yd, entry, ORDER, knots, t ){

    ## Calculate B-Spline Basis Functions (B-values) in terms of 'x'
    # We cannot pass any NaN values, these must be set aside
    x = as.numeric(yd[,1])
    sx = sort(x)
    si = order(x)
    inc = which(!is.na(sx))
    n = length(x)

    knots = unique( as.numeric(knots) )
    Boundary.knots = c(min(knots),max(knots))
    Interior.knots = unique(knots)[2:(length(knots)-1)]

    # Basis functions for h
    Wik = NaN*zeros(n, length(knots)-ORDER)
    if( length(inc) ){
      Wik[si[inc], ] = #spcol(knots, order, sx[inc], 'noderiv')
        bSpline( x=sx[inc], knots=Interior.knots, Boundary.knots=Boundary.knots, degree=(ORDER-1), intercept=TRUE)
    }

    # Basis functions for the definite integral of h
    # * Lower limit of integration is the left endpoint of the basic interval
    # * This code is based on spcol and De Boor, PGS, p. 150-151.
    nk = length(knots)
    # Number of actual parameters
    p = nk - ORDER
    scl = (knots[(ORDER+1):(ORDER+p)] - knots[1:p])/ORDER
    z = (tril(kron(1:p*ones(1,p), ones(p,1))) > 0) * (ones(p,1)%*%scl)
    Zik = NaN*zeros(n,p)
    if( length(inc) ){
      Zik[ si[inc], ] = bSpline( x=sx[inc], knots=Interior.knots, Boundary.knots=Boundary.knots, degree=ORDER, intercept=TRUE)[,-1]
    }
    Zik = Zik%*%z

    # Basis functions for our vector of evaluation points t - hazard
    x = as.numeric( t )
    sx = sort(x)
    si = order(x)
    inc = which(!is.na(sx))
    n = length(x)
    Xh = NaN*zeros(n,length(knots) - ORDER)

    if( !length(inc) ){
      print( "ERROR : You have no data to analyze" )
    }
    Xh[ si[inc], ] = #spcol(knots, order, sx(inc), 'noderiv')
      bSpline( x=sx[inc], knots=Interior.knots, Boundary.knots=Boundary.knots, degree=(ORDER-1), intercept=TRUE)

    # Basis functions t - define integral
    XH = NaN*zeros(n,p);
    XH[si[inc],] = #spcol([knots, knots(end)], order+1, sx(inc), 'noderiv');
      bSpline( x=sx[inc], knots=Interior.knots, Boundary.knots=Boundary.knots, degree=ORDER, intercept=TRUE)[,-1]
    XH = XH%*%z;

    return( list(
      Wik = Wik,
      Zik = Zik,
      XH = XH,
      Xh = Xh
    ))
  }
