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
#' @return
#' @importFrom matrixStats rowQuantiles
#'
#' @details See the \href{https://doi.org/10.2307/2532989}{Original paper by Philip S. Rosenberg}.
#'
#' @export
#'
#' @examples
hspcore <- function(yd, ORDER, knots, time, Bootstrap=0, alphalevel=0.95){

  knots = unique( as.numeric(knots) )

  ## Rescale the time axis to avoid numerical problems
  TMAX = max(yd[,1]);
  yd[,1] = yd[,1]/TMAX;
  entry = entry/TMAX;
  knots = knots/TMAX;
  t = t/TMAX;

  ## Calculate all needed auxiliary matrices and functions
  basis_functions = bspline_regression_basis_functions(yd, entry, ORDER, knots, t )
  Wik = basis_functions$Wik
  Zik = basis_functions$Zik
  XH  = basis_functions$XH
  Xh  = basis_functions$Xh

  ## Initial guess for the coefficients
  alpha0 = (sum(yd[,2])/sum(yd[,1]))*ones(1, size(Wik,2))

  ## Maximize the likelihood function by minimising ( - 2 log-Likekihood )
  maxL = hazl_ker(yd, alpha0, Wik, Zik, Xh, XH )

  alpha1= maxL$alpha1
  h     = maxL$h
  S     = maxL$S
  m2loglik = maxL$m2loglik
  m2loglik_nosmooth = srllikb_fun(alpha1, yd, Wik, Zik)


  ####
  # Bootstrap - keep track of h(t), S(t), and alpha
  # We are just sampling with replacement
  # Use the MLE as the starting value
  ####
  if( B>0 ){

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
    S = cbind( S, rowQuantiles( Sb, probs=c(alphalevel/2, 1-alphalevel/2), na.rm=T ) )

  } else { # no bootstrapping
    alphab = c();
    hb = c();
    Sb = c();
  }

  # rescale the hazard values and coefficients
  h = h/TMAX;
  alpha1 = alpha1/TMAX;

  if( B > 0 ){

    hb = hb/TMAX;
    alphab = alphab/TMAX;
  }

  return( list(alpha1=alpha1,
               t=t*TMAX,
               h=h,
               S=S,
               m2loglik=m2loglik,
               m2loglik_nosmooth=m2loglik_nosmooth,
               hb=hb) )

}
