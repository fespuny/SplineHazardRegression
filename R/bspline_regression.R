#' Core function for B-spline hazard regression
#'
#' @description This function allows to perform B-spline hazard regression given right-censored time-to-event data.
#'              If no b-spline knots and/or b-spline order are provided, automatic selection of those parameters is performed (see Details).
#'
#' @param yd - matrix of time and status data; y is event time, d is for delta, the status indicator
#' @param ORDER - 1 step, 2 linear, 3 quadratic, 4 cubic (i.e. ORDER=degree-1)
#' @param knots - optional sequence of knot locations including endpoints (without multiplicity)
#' @param time - vector of evaluation times
#' @param Bootstrap - number of bootstrap samples (set to zero if none)
#' @param alphalevel - alpha level for confidence intervals if bootrap is being used
#'
#' @return list(alpha1,t,h,m2loglik,hb) - the optimal coefficients alpha1, the time variable, the optimal hazard, the objective function value and gradient, and the bootstrap hazards if any
#' @importFrom matrixStats rowQuantiles
#' @importFrom pracma size
#'
#' @details See \doi{<doi.org/10.2307/2532989>} the original paper by Philip S. Rosenberg.
#'
#' @export
#'
hspcore <- function(yd, ORDER=4, knots, time, Bootstrap=0, alphalevel=0.95){

  knots = unique( as.numeric(knots) )

  ## Rescale the time axis to avoid numerical problems
  TMAX = max(yd[,1]);
  yd[,1] = yd[,1]/TMAX;
  knots = unique(knots)/TMAX;
  t = time/TMAX;

  if( !(ORDER %in% 1:4) ){
    print( "Error: ORDER needs to take a value between 1 and 4; see help('hspcore')" )
    return( -1 )
  }

  ## AUTOMATIC SEARCH FOR KNOTS NEEDED?
  if( length(knots)<2 ) {
    print( "Automatic search for K the number of interior knots of the B-spline hazard function" )

    Exterior.knots = c( min(yd[,1]), max(yd[,1]) )
    eventtimes = sort( yd[ yd[,2]==1, 1] )
    nevents = length( eventtimes )

    bestAICc = NULL
    bestK    = NULL

    for( K in 1:8 ){

      DOF = K+ORDER

      ## interior knots using K percentiles for the event times
      Interior.knots = eventtimes[ round( seq( nevents/(K+1), nevents/(K+1) * K , nevents/(K+1) ) ) ]

      ## candidate knots
      knots = c( min(eventtimes), Interior.knots, max(eventtimes) )
      print( paste0( "K= ", K, " knots= ", paste0( round( knots, 1), collapse = " "  ) ) )

      ## Calculate all needed auxiliary matrices and functions
      basis_functions = bspline_regression_basis_functions(yd[,1], ORDER, knots, t )

      Wik = basis_functions$Wik
      Zik = basis_functions$Zik
      Xh  = basis_functions$Xh
      XH  = basis_functions$XH
      print( paste0( "K= ", K, " knots= ", paste0( round( knots, 1), collapse = " "  ) ) )

      ## Initial guess for the coefficients
      alpha0 = (sum(yd[,2])/sum(yd[,1]))*ones(1, size(Wik,2))

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
    print( paste0( "We use ", bestK , " interior B-spline knots"))

    Interior.knots = eventtimes[ round( seq( nevents/(bestK+1), nevents/(bestK+1) * bestK , nevents/(bestK+1) ) ) ]
    knots = c( min(eventtimes), Interior.knots, max(eventtimes) )
  }

  ## REGRESSION WITH KNOWN KNOTS
  print( paste0( "K= ", length(knots)-2, " DOF= ", length(knots)-2+ORDER, " knots= ", paste0( round( knots, 1), collapse = " "  ) ) )

    ## Calculate all needed auxiliary matrices and functions
    basis_functions = bspline_regression_basis_functions(yd[,1], ORDER, knots, t )

    Wik = basis_functions$Wik
    Zik = basis_functions$Zik
    Xh  = basis_functions$Xh
    XH  = basis_functions$XH

    ## Initial guess for the coefficients
    alpha0 = (sum(yd[,2])/sum(yd[,1]))*ones(1, size(Wik,2))

    ## Maximize the likelihood function by minimising the objective = -2*logLikekihood
    maxL = hazl_ker(yd, alpha0, Wik, Zik, Xh, XH )

    ## SOLUTION GIVEN THE KNOTS
    alpha1= maxL$alpha1
    h     = maxL$h
    m2loglik = maxL$m2loglik #this contains

  ## BOOTSTRAP NEEDED?
  if( Bootrsrap > 0 ){
  # Bootstrap - keep track of h(t), S(t), and alpha
  # We are just sampling with replacement
  # Use the MLE as the starting value?

    hb = zeros(length(t),  B);
    Sb = zeros(length(t), B);
    alphab = zeros(B, length(alpha1));
    m2loglikb = zeros(B, 1);

    for( i in 1:B ){

      iteration = (1 - i) %/% 5 #5-fold cross validation for each iteration

      if( (i-1)%%5 ==0 ) shuffle = sample( x=1:size(yd,1) ) # index permutation

      i_b = shuffle[ which( (0:(size(yd, 1)-1))%%5 != (i-1)%%5) ] #only 80% of the indices are selected

      # i_b = ceil(size(yd, 1)*rand(size(yd, 1), 1));
      yd0 = yd[i_b,]
      Wik0 = Wik[i_b, ]
      Zik0 = Zik[i_b, ]
      # Eik0 = Eik[i_b, ]
      Eik0 = NULL
      maxLbi = hazl_ker( yd0, alpha0, Wik0, Zik0, Eik0, Xh, XH, smooth );

      alphab[i,] = maxLbi$alpha1
      hb[,i] = maxLbi$h
      Sb[,i] = maxLbi$S
      m2loglikb[i] = maxLbi$m2loglik
      rm( maxLbi )
    }
    gc()

    ###
    # Package the results / pointwise confidence limits
    # We are saving trios of columns / MLE, lower limit, upper limit
    ###

    # library( matrixStats )
    h = cbind( h, rowQuantiles( hb, probs=c(alphalevel/2, 1-alphalevel/2), na.rm = T ) )

  } else { # no bootstrapping
    alphab = c();
    hb = c();
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
               m2loglik=m2loglik,
               hb=hb) )

}
