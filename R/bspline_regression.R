#' Title
#'
#' @param yd
#' @param ORDER
#' @param knots
#' @param time
#' @param Bootstrap
#' @param alphalevel
#'
#' @return
#' @importFrom matrixStats rowQuantiles
#' @export
#'
#' @examples
hspcore <- function(yd, ORDER, knots, time, Bootstrap, alphalevel){
  # % HSPCORE core function for hazard function estimation using B-Splines.
  # %
  # % [alpha1, t, h, S, m2loglik] = hspcore(yd, entry, order, knots, t, B, alphalevel)
  #yd pneumonic y is event time, d is for delta, the status indicator
  # yd    - matrix of time and status data
  # order - 1 step, 2 linear, 3 quadratic, 4 cubic
  # knots - sequence of knot locations including enpoint multiplicities
  # t     - vector of evaluation times
  # B     - number of bootstrap samples
  # alpha - alpha level

  ## Rescale the time axis to avoid numerical problems
  TMAX = max(yd[,1]);
  yd[,1] = yd[,1]/TMAX;
  entry = entry/TMAX;
  knots = knots/TMAX;
  t = t/TMAX;

  basis_functions = bspline_regression_basis_functions(yd, entry, ORDER, knots, t )

  # An initial guess for the coefficients
  #MATLAB alpha0 = (sum(yd[,2])/sum(yd[,1]-entry))*ones(1, size(Wik,2));
  alpha0 = (sum(yd[,2])/sum(yd[,1]))*ones(1, size(Wik,2))

  # Maximize the likelihood function
  maxL = hazl_ker(yd, alpha0, Wik, Zik, Eik, Xh, XH, smooth);
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
