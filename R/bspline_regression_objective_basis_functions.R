#' Generate the basis functions needed for fitting the hazard B-spline model
#'
#' @description This function calculates all matrices needed for the estimation of a hazard function with censored time-to-event data using B-splines
#'
#' @param tin - event/censoring times in data
#' @param ORDER - 1 step, 2 linear, 3 quadratic, 4 cubic
#' @param knots - sequence of knot locations
#' @param t - vector of evaluation times
#'
#' @return list(Wik, Zik, Xh, XH ) - see details for definitions
#'
#' @details When fitting a B-spline function \eqn{h(\alpha)(t) := B(t) \alpha} to time-to-event-data,
#' where \eqn{B(t)} is a basis of B-splines and \eqn{\alpha} are the coefficients to be estimated,
#' we optimise a function of the likelihood \deqn{L(h(\theta)) =  }.
#'
#' @importFrom pracma zeros ones tril kron
#' @export
#'
bspline_regression_basis_functions <- function(tin, ORDER, knots, t ){

    ## Calculate B-Spline Basis Functions (B-values) in terms of 'x'
    # We cannot pass any NaN values, these must be set aside
    x = as.numeric(tin)
    sx = sort(x)
    si = order(x)
    inc = which(!is.na(sx))
    n = length(x)

    knots = unique( as.numeric(knots) )
    Boundary.knots = c(min(knots),max(knots))
    Interior.knots = setdiff(knots,Boundary.knots)
    p = length(Interior.knots)+ORDER     # Number of actual parameters (DOF)

    # Basis functions for h
    Wik = NaN*zeros(n, p )
    if( length(inc) ){
      Wik[si[inc], ] =
#        bSpline( x=sx[inc], knots=Interior.knots, Boundary.knots=Boundary.knots, degree=(ORDER-1), intercept=TRUE)
        generate_bspline_basis( time=sx[inc], Interior.knots, Boundary.knots, ORDER )
    }

    # Basis functions for the definite integral of h
    # * Lower limit of integration is the left endpoint of the basic interval
    # * This code is based on spcol and De Boor, PGS, p. 150-151.
    knots_with_multiplicity = c( rep(Boundary.knots[1],ORDER), Interior.knots, rep(Boundary.knots[2],ORDER) )
    scl = (knots_with_multiplicity[(ORDER+1):(ORDER+p)] - knots_with_multiplicity[1:p])/ORDER
    z = (tril(kron(1:p*ones(1,p), ones(p,1))) > 0) * (ones(p,1)%*%scl)
    Zik = NaN*zeros(n,p)
    if( length(inc) ){
      Zik[ si[inc], ] =
        # bSpline( x=sx[inc], knots=Interior.knots, Boundary.knots=Boundary.knots, degree=ORDER, intercept=TRUE)[,-1]
        generate_bspline_basis( time=sx[inc], Interior.knots, Boundary.knots, (ORDER+1) )[,-1]
    }
    Zik = Zik%*%z

    # Basis functions for our vector of evaluation points t - hazard
    x = as.numeric( t )
    sx = sort(x)
    si = order(x)
    inc = which(!is.na(sx))
    n = length(x)
    Xh = NaN*zeros(n,p)

    if( !length(inc) ){
      print( "ERROR : You have no data to analyze" ) ###//HERE HERE HERE
    }
    Xh[ si[inc], ] =
#     bSpline( x=sx[inc], knots=Interior.knots, Boundary.knots=Boundary.knots, degree=(ORDER-1), intercept=TRUE)
    generate_bspline_basis( time=sx[inc], Interior.knots, Boundary.knots, ORDER )

    # Basis functions t - define integral
    XH = NaN*zeros(n,p);
    XH[si[inc],] =
      # bSpline( x=sx[inc], knots=Interior.knots, Boundary.knots=Boundary.knots, degree=ORDER, intercept=TRUE)[,-1]
    generate_bspline_basis( time=sx[inc], Interior.knots, Boundary.knots, (ORDER+1) )[,-1]
    XH = XH%*%z;

    return( list(
      Wik = Wik,
      Zik = Zik,
      Xh = Xh,
      XH = XH
    ))
  }
