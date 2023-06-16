#' Optimisation of the objective function
#'
#' @param yd - data
#' @param alpha0 - initial solution value for the B-spline coefficients
#' @param Wik - matrix as calculated by bspline_regression_basis_functions()
#' @param Zik - matrix as calculated by bspline_regression_basis_functions()
#' @param Xh - matrix as calculated by bspline_regression_basis_functions()
#' @param XH - matrix as calculated by bspline_regression_basis_functions()
#'
#' @return list(alpha1,hout,Sout,m2loglik) - alpha1 is the new solution;
#' using alpha1, the function also calculates the hazard (hout), survival (Sout) and objective function and gradient (m2loglik) ;
#' @importFrom pracma zeros
#' @importFrom stats optim nlminb
#' @export
#'
hazl_ker = function(yd, alpha0, Wik, Zik, Xh, XH ){
  # wrapper function for fmincon's access to srllikb
  # MATLAB hoptions = optimset('TolFun', 1E-6, 'TolX', 1E-6, 'Display', 'none', ...
  #     'algorithm', 'interior-point', 'GradObj', 'on');
  p = size(Wik, 2)
  lwr = zeros(1, p);

  #OPTIMISATION WITH CONSTRAINTS (we compare two different R functions based on Fortran code)
  output_optim = optim( as.numeric(alpha0), srllikb_fun, srllikb_grad, yd, Wik, Zik, method = "L-BFGS-B", lower=0)
  alpha1 = as.numeric(output_optim$par) #array (interpreted as column vector in matrix operations)

    output_optim2 = nlminb ( start=as.numeric(alpha0), objective=srllikb_fun, gradient=srllikb_grad, yd, Wik, Zik, lower=rep(0,length(alpha0)))
    alpha2 = as.numeric(output_optim2$par) #array (interpreted as column vector in matrix operations)

    if( output_optim2$objective < output_optim$objective ){
      print( "Optimisation note: The PORT optimization routine was superior to the L-BFGS-B")
      alpha1 <- alpha2
    }

  hout = Xh %*% alpha1;
  Sout = exp(-XH %*% alpha1 );

  m2loglik = c(srllikb_fun(alpha1, yd, Wik, Zik),
               srllikb_grad(alpha1, yd, Wik, Zik))

  return( list(alpha1=alpha1,
               h=hout,
               S=Sout,
               m2loglik=m2loglik) )
}
