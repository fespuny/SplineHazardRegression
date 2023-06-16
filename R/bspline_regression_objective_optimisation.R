#' Optimisation of the objective function
#'
#' @param yd - data
#' @param alpha0 - initial solution value for the B-spline coefficients
#' @param Wik - matrix as calculated by bspline_regression_basis_functions()
#' @param Zik - matrix as calculated by bspline_regression_basis_functions()
#' @param Eik - matrix as calculated by bspline_regression_basis_functions()
#' @param Xh - matrix as calculated by bspline_regression_basis_functions()
#' @param XH - matrix as calculated by bspline_regression_basis_functions()
#' @param smooth
#'
#' @return list(alpha1,hout,Sout,m2loglik) - alpha1 is the new solution;
#' using alpha1, the function also calculates the hazard (hout), survival (Sout) and objective function and gradient (m2loglik) ;
#' @importFrom pracma zeros
#' @importFrom stats optim
#' @export
#'
hazl_ker = function(yd, alpha0, Wik, Zik, Eik, Xh, XH, smooth=FALSE){
  # wrapper function for fmincon's access to srllikb
  # MATLAB hoptions = optimset('TolFun', 1E-6, 'TolX', 1E-6, 'Display', 'none', ...
  #     'algorithm', 'interior-point', 'GradObj', 'on');
  p = size(Wik, 2)
  lwr = zeros(1, p);

  #[WE NEED TO TEST OTHER OPTIMISATION ROUTINES]
  output_optim = optim( as.numeric(alpha0), srllikb_fun, srllikb_deriv, yd, Wik, Zik, method = "L-BFGS-B", lower=0)

  alpha1 = as.numeric(output_optim$par) #array (interpreted as column vector in matrix operations)

  hout = Xh %*% alpha1;
  Sout = exp(-XH %*% alpha1 );

  if( smooth == FALSE){ #default

    m2loglik = c(srllikb_fun(alpha1, yd, Wik, Zik),
                 srllikb_deriv(alpha1, yd, Wik, Zik))
  } else {
    m2loglik = c(srllikb_smooth_fun(alpha1, yd, Wik, Zik),
                 srllikb_smooth_deriv(alpha1, yd, Wik, Zik))
  }


  return( list(alpha1=alpha1,
               h=hout,
               S=Sout,
               m2loglik=m2loglik) )
}
