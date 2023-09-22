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
    Wik = NaN*matrix(0,n, p )
    if( length(inc) ){
      Wik[si[inc], ] =
#        bSpline( x=sx[inc], knots=Interior.knots, Boundary.knots=Boundary.knots, degree=(ORDER-1), intercept=TRUE)
        generate_bspline_basis( time=sx[inc], Interior.knots, Boundary.knots, ORDER )
    }

    # Basis functions for the definite integral of h
    Zik = NaN*matrix(0,n,p)
    if( length(inc) ) {
      Zik[ si[inc], ] = generate_bspline_basis( time=sx[inc], Interior.knots, Boundary.knots, ORDER, integral=TRUE )
    }
    # if( 0 ){ ##comparison to the De Boor method (this required the pracma package)
    #   # * Lower limit of integration is the left endpoint of the basic interval
    #   # * This code is based on spcol and De Boor, PGS, p. 150-151.
    #   knots_with_multiplicity = c( rep(Boundary.knots[1],ORDER), Interior.knots, rep(Boundary.knots[2],ORDER) )
    #   scl = (knots_with_multiplicity[(ORDER+1):(ORDER+p)] - knots_with_multiplicity[1:p])/ORDER
    #   z = (pracma::tril(pracma::kron(1:p*matrix(1,1,p), matrix(1,p,1))) > 0) * (matrix(1,p,1)%*%scl)
    #   Zik_DB = NaN*matrix(0,n,p)
    #   if( length(inc) ){
    #     Zik_DB[ si[inc], ] =
    #       # bSpline( x=sx[inc], knots=Interior.knots, Boundary.knots=Boundary.knots, degree=ORDER, intercept=TRUE)[,-1]
    #       generate_bspline_basis( time=sx[inc], Interior.knots, Boundary.knots, (ORDER+1) )[,-1]
    #   }
    #   Zik_DB = Zik_DB%*%z
    #   print( paste0( "Relative difference between integral basis methods: ", sum( abs(Zik - Zik_DB) / sum(abs(Zik_DB)) ) ) )
    # }

    # Basis functions for our vector of evaluation points t - hazard
    x = as.numeric( t )
    sx = sort(x)
    si = order(x)
    inc = which(!is.na(sx))
    n = length(x)
    Xh = NaN*matrix(0,n,p)

    if( !length(inc) ){
      print( "ERROR : You have no data to analyze" ) ###//HERE HERE HERE
    }
    Xh[ si[inc], ] =
#     bSpline( x=sx[inc], knots=Interior.knots, Boundary.knots=Boundary.knots, degree=(ORDER-1), intercept=TRUE)
    generate_bspline_basis( time=sx[inc], Interior.knots, Boundary.knots, ORDER )

    # Basis functions t - define integral
    XH = NaN*matrix(0,n,p);
    XH[si[inc],] =
      generate_bspline_basis( time=sx[inc], Interior.knots, Boundary.knots, ORDER, integral=TRUE )

    # if( 0 ){ ##comparison to the De Boor method (this required the pracma package)
    #   # # * This code is based on spcol and De Boor, PGS, p. 150-151.
    #   knots_with_multiplicity = c( rep(Boundary.knots[1],ORDER), Interior.knots, rep(Boundary.knots[2],ORDER) )
    #   scl = (knots_with_multiplicity[(ORDER+1):(ORDER+p)] - knots_with_multiplicity[1:p])/ORDER
    #   z = (tril(kron(1:p*matrix(1,1,p), matrix(1,p,1))) > 0) * (matrix(1,p,1)%*%scl)
    #   XH_DB = NaN*matrix(0,n,p);
    #   XH_DB[si[inc],] =
    #     generate_bspline_basis( time=sx[inc], Interior.knots, Boundary.knots, (ORDER+1) )[,-1]
    #   XH_DB = XH_DB%*%z;
    #
    #   print( paste0( "Relative difference between integral basis methods: ", sum( abs(XH - XH_DB) ) / sum(abs(XH_DB)) ) )
    # }

    return( list(
      Wik = Wik,
      Zik = Zik,
      Xh = Xh,
      XH = XH
    ))
  }
