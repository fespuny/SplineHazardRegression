#' Core function for B-spline hazard regression
#'
#' @description This function allows to perform B-spline hazard regression given right-censored time-to-event data.
#'              If no b-spline knots and/or b-spline order are provided, automatic selection of those parameters is performed (see Details).
#'
#' @param yd - matrix of time and status data; y is event time, d is for delta, the status indicator
#' @param ORDER - 1 step, 2 linear, 3 quadratic, 4 cubic (i.e. ORDER=degree-1)
#' @param Exterior.knots - start and end knot locations (without multiplicity)
#' @param Interior.knots - optional sequence of interior knot locations excluding endpoints (without multiplicity)
#' @param SelectBestKnots - Boolean; if TRUE an exploration for the best number of interior knots is performed
#' @param time - vector of evaluation times
#' @param Bootstrap - number of bootstrap samples (set to zero if none)
#' @param alphalevel - alpha level for confidence intervals if bootstrap is being used
#' @param verbose - boolean; if TRUE, the convergence of the two used optimisation routines is reported during bootstrap estimation of variance
#'
#' @return list(alpha1,t,h,m2loglik,hb) - the optimal coefficients alpha1, the time variable, the optimal hazard, the objective function value and gradient, and the bootstrap hazards if any
#' @importFrom pracma zeros ceil
#' @importFrom stats runif quantile
#'
#' @details See \doi{<doi.org/10.2307/2532989>} the original paper by Philip S. Rosenberg.
#'
#' @export
#'
hspcore <- function(yd, ORDER=4, Exterior.knots, Interior.knots=NULL, SelectBestKnots, time, Bootstrap=0, alphalevel=0.05, verbose=FALSE){

  Exterior.knots = unique( as.numeric(Exterior.knots) )
  if( length( Exterior.knots ) != 2 ){
    print( "ERROR: Exterior.knots should contain 2 numeric values")
    return(-1)
  }

  ## Rescale the time axis to avoid numerical problems
  TMAX = max(yd[,1]);
  yd[,1] = yd[,1]/TMAX;
  Exterior.knots = unique(Exterior.knots)/TMAX;
  Interior.knots = unique(Interior.knots)/TMAX;
  t = time/TMAX;

  if( !(ORDER %in% 1:4) ){
    print( "Error: ORDER needs to take a value between 1 and 4; see help('hspcore')" )
    return( -1 )
  }

  ## AUTOMATIC SEARCH FOR KNOTS NEEDED?
  if( SelectBestKnots==TRUE ) {
    print( "Automatic search for K the number of interior knots of the B-spline hazard function" )

    eventtimes = sort( yd[ yd[,2]==1, 1] )
    nevents = length( eventtimes )

    bestAICc = NULL
    bestK    = NULL

    for( K in 1:8 ){

      DOF = K+ORDER

      ## interior knots using K percentiles for the event times
      Interior.knots = eventtimes[ round( seq( nevents/(K+1), nevents/(K+1) * K , nevents/(K+1) ) ) ]

      ## candidate knots
      knots = c( Exterior.knots[1], Interior.knots, Exterior.knots[2] )

      ## Calculate all needed auxiliary matrices and functions
      basis_functions = bspline_regression_basis_functions(yd[,1], ORDER, knots, t )

      Wik = basis_functions$Wik
      Zik = basis_functions$Zik
      Xh  = basis_functions$Xh
      XH  = basis_functions$XH

      ## Initial guess for the coefficients
      alpha0 = (sum(yd[,2])/sum(yd[,1]))*ones(1, ncol(Wik))

      ## Maximize the likelihood function by minimising the objective = -2*logLikekihood
      maxL = hazl_ker(yd, alpha0, Wik, Zik, Xh, XH )

      ## SOLUTION GIVEN THE KNOTS
      alpha1= maxL$alpha1
      h     = maxL$h
      S     = maxL$S
      m2loglik = maxL$m2loglik #this contains

      AICc = m2loglik[1] + 2 * DOF + 2 * DOF * (DOF+1) / (nrow(yd)-DOF-1)

      print( paste0( "K= ", K, " AICc=", AICc, " knots= ", paste0( round( knots, 1), collapse = " "  ) ) )

      if( K==1 ){
        bestAICc = AICc #we initialise bestAICc
        bestK    = K
      }

      if( AICc < bestAICc ){
        bestAICc = AICc #we minimise AICc
        bestK    = K
      }
    }
    ## we use the best K found:
    print( paste0( "SEARCH RESULT: We use ", bestK , " interior B-spline knots"))

    Interior.knots <- eventtimes[ round( seq( nevents/(bestK+1), nevents/(bestK+1) * bestK , nevents/(bestK+1) ) ) ]
  }

  knots = c( Exterior.knots[1], Interior.knots, Exterior.knots[2] )

  ## REGRESSION WITH KNOWN KNOTS
  print( paste0( "K= ", length(knots)-2, " DOF= ", length(knots)-2+ORDER, " knots= ", paste0( round( knots, 1), collapse = " "  ) ) )

    ## Calculate all needed auxiliary matrices and functions
    basis_functions = bspline_regression_basis_functions(yd[,1], ORDER, knots, t )

    Wik = basis_functions$Wik
    Zik = basis_functions$Zik
    Xh  = basis_functions$Xh
    XH  = basis_functions$XH

    ## Initial guess for the coefficients
    alpha0 = (sum(yd[,2])/sum(yd[,1]))*ones(1, ncol(Wik))

    ## Maximize the likelihood function by minimising the objective = -2*logLikekihood
    maxL = hazl_ker(yd, alpha0, Wik, Zik, Xh, XH )

    ## SOLUTION GIVEN THE KNOTS
    alpha1= maxL$alpha1
    h     = maxL$h
    S     = maxL$S
    m2loglik = maxL$m2loglik #this contains

  ## BOOTSTRAP NEEDED?
  if( Bootstrap > 0 ){
  # Bootstrap - keep track of h(t), S(t), and alpha
  # We are just sampling with replacement
  print( "Variance estimation using bootstrap" )

    hb = zeros(length(t),  Bootstrap);
    Sb = zeros(length(t), Bootstrap);
    alphab = zeros(Bootstrap, length(alpha1));
    m2loglikb = zeros(Bootstrap, 1);
    convergence = zeros(Bootstrap, 3);

    for( i in 1:Bootstrap ){

      #i_b = unique( ceil( runif( n=length(t), min=0, max=nrow(yd) ) ) )
      i_b = ceil( runif( n=length(t), min=0, max=nrow(yd) ) )
      yd0 = yd[i_b,]
      Wik0 = Wik[i_b, ]
      Zik0 = Zik[i_b, ]
      alpha0i = (sum(yd0[,2])/sum(yd0[,1]))*ones(1, ncol(Wik0))
      maxLbi = hazl_ker( yd0, alpha0i, Wik0, Zik0, Xh, XH, verbose );

      alphab[i,] = maxLbi$alpha1
      hb[,i] = maxLbi$h
      Sb[,i] = maxLbi$S
      m2loglikb[i] = maxLbi$m2loglik[1]
      convergence[i,] = as.numeric( maxLbi$convergence )
      rm( maxLbi )
    }
    gc()

    if( verbose ){
      print( "1. Convergence L-BFGS-B (0=yes):" )
      print( table( convergence[,1] ) )
      print( "2. Convergence PORT (0=yes):" )
      print( table( convergence[,2] ) )
      print( "3. Winning methods:" )
      print( table( convergence[,3] ) )
    }

    ###
    # Package the results / pointwise confidence limits
    # We are saving trios of columns / MLE, lower limit, upper limit
    ###

    # library( matrixStats )
    # h = cbind( h, rowQuantiles( hb, probs=c(alphalevel/2, 1-alphalevel/2), na.rm=T ) )
    # S = cbind( S, rowQuantiles( Sb, probs=c(alphalevel/2, 1-alphalevel/2), na.rm=T ) )
    hci = hb
    Sci = Sb
    for( j in 1:nrow(hb) ){
      hci[j,] = sort( hb[j,] )
      Sci[j,] = sort( Sb[j,] )
    }
    ialphaci = round( quantile(1:Bootstrap, c( alphalevel/2, 1-alphalevel/2) ) )
    h = cbind(h, hci[,ialphaci])
    S = cbind(S, Sci[,ialphaci])

  } else { # no bootstrapping
    alphab = c();
    hb = c();
    Sb=c()
  }

  # rescale the hazard values and coefficients
  h = h/TMAX;
  alpha1 = alpha1/TMAX;


  if( Bootstrap > 0 ){

    hb = hb/TMAX;
    alphab = alphab/TMAX;
  }

  return( list(alpha1=alpha1,
               t=t*TMAX,
               h=h,
               S=S,
               m2loglik=m2loglik,
               hb=hb,
               convergenceb = convergence ) )

}
