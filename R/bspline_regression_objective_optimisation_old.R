#' Optimisation of the objective function
#'
#' @param yd - data
#' @param alpha0 - initial solution value for the B-spline coefficients
#' @param Wik - matrix as calculated by bspline_regression_basis_functions()
#' @param Zik - matrix as calculated by bspline_regression_basis_functions()
#' @param Xh - matrix as calculated by bspline_regression_basis_functions()
#' @param XH - matrix as calculated by bspline_regression_basis_functions()
#' @param verbose - Boolean; if TRUE, the convergence of the two used optimisation routines is reported.
#'
#' @return list(alpha1,h,S,m2loglik) - alpha1 is the new solution;
#' using alpha1, the function also calculates the hazard (h), survival (S) and objective function and gradient (m2loglik) ;
#' @importFrom stats optim nlminb
#' @export
#'
hazl_ker_old = function(yd, alpha0, Wik, Zik, Xh, XH, verbose=FALSE ){
  # wrapper function for fmincon's access to srllikb
  # MATLAB hoptions = optimset('TolFun', 1E-6, 'TolX', 1E-6, 'Display', 'none', ...
  #     'algorithm', 'interior-point', 'GradObj', 'on');
  p = ncol(Wik)
  lwr = rep(0, p);

  #OPTIMISATION WITH CONSTRAINTS (we compare two different R functions based on Fortran code)
  ## L-BFGS-B
  output_optim = list( par = alpha0, convergence=-1, value = srllikb_fun(alpha0, yd, Wik, Zik), message="" )   #initial values
  tryCatch({
    output_optim = optim( as.numeric(alpha0), srllikb_fun, srllikb_grad, yd, Wik, Zik, method = "L-BFGS-B", lower=0)
    }, error= function(e){ output_optim$message <<- conditionMessage(e) } )

  # print( paste0( "L-BFGS-B Convergence (0=yes,1=maxiter)? ", output_optim$convergence ))
  # if( output_optim$convergence > 50 ){ print( paste0("L-BFGS-B Warning: ", output_optim$message))}
  alpha1 = as.numeric(output_optim$par) #array (interpreted as column vector in matrix operations)

    output_optim2 = nlminb( start=as.numeric(alpha0), objective=srllikb_fun, gradient=srllikb_grad, hessian=NULL,
                            yd, Wik, Zik, lower=lwr)
    alpha2 = as.numeric(output_optim2$par) #array (interpreted as column vector in matrix operations)
    # print( paste0( "PORT Convergence (0=yes)? ", output_optim2$convergence ))
    # if( output_optim2$convergence != 0 ){ print( paste0("PORT Warning: ", output_optim2$message))}
    convergence = c( output_optim$convergence,
                     output_optim2$convergence,
                     1 )

    if( output_optim2$convergence==0 ){
      if( output_optim2$objective < output_optim$value ){
        alpha1 <- alpha2
        convergence[3] = 2
      }
    }

  if( convergence[1] * convergence[2] != 0 ){

    print( "*ERROR* none of the tested optimization methods did converge ")
    convergence[3] = 3
    alpha1 <<- alpha0


    ### we test R functions for non-linear constrained optimization
    output_optim3 <- list( par = alpha0, convergence=-1, value = srllikb_fun(alpha0, yd, Wik, Zik), message="" )
    output_optim4 <- list( par = alpha0, convergence=-1, value = srllikb_fun(alpha0, yd, Wik, Zik), message="" )
    applyDefaults <- function(fn, ...) { function(x) fn(x, ...) }
    Lalpha1 = srllikb_fun(alpha1, yd, Wik, Zik)

    ## Adaptative barrier algorithm using BFGS as optimization method (Nelder-Mead did not find a solution to trouble problems)
    tryCatch({
    output_optim3 =stats::constrOptim( theta = as.numeric(alpha0), f = applyDefaults(srllikb_fun,yd=yd,Wik=Wik,Zik=Zik), grad=applyDefaults(srllikb_grad,yd=yd,Wik=Wik,Zik=Zik),
                                       ui=pracma::eye(p), ci=rep(0,p), control = list(fnscale = 1) )
    if( output_optim3$convergence == 0 ){ print( "Adaptative barrier with BFGS converged!") } else { print( "Adaptative barrier with BFGS failed") }
    }, error= function(e){ output_optim3$message <<- conditionMessage(e) } )
    if( output_optim3$convergence == 0 ){ #the previous algorithms did not converge, so this is our best estimate
      if( sum(output_optim3$par>=0) == length( output_optim3$par ) ){ ## alpha >= 0 ?
        alpha1  <- output_optim3$par
        Lalpha1 <- srllikb_fun(alpha1, yd, Wik, Zik)
        print( paste0( "Best B-spline ", paste( round(alpha1,3), collapse = ","), " with objective ", Lalpha1 ) )
      }
    }

    ## BOBYQA, Bound Optimization by Quadatric Approximation
    tryCatch({
      output_optim4 =nloptr::bobyqa( x0 = as.numeric(alpha0), fn = applyDefaults(srllikb_fun,yd=yd,Wik=Wik,Zik=Zik), lower=rep(0,p) )
      if( output_optim4$convergence > 0 ){ print( "BOBYQA, Bound Optimization by Quadatric Approximation converged!") } else { print( "BOBYQA, Bound Optimization by Quadatric Approximation failed") }
    }, error= function(e){ output_optim4$message <<- conditionMessage(e) } )
    if( output_optim4$convergence == 0 ){ #the previous algorithms did possibly converge, so is this a better estimate?
      if( sum(output_optim4$par>=0) == length( output_optim4$par ) ){ ## alpha >= 0 ?
        if( output_optim4$value < Lalpha1 ){
          alpha1  <- output_optim4$par
          Lalpha1 <- srllikb_fun(alpha1, yd, Wik, Zik)
          printt( paste0( "Best B-spline ", paste( round(alpha1,3), collapse = ","), " with objective ", Lalpha1 ) )
        }
      }
    }
  }

  ##report convergence properties
  if( verbose ){
    if( output_optim$convergence == 0 & output_optim2$convergence== 0 ){
      print( paste0( "Objective L-BFGS-B: ", output_optim$value, " Objective PORT: ", output_optim2$objective) )
    } else {
      if( output_optim$convergence == 0 ) print( paste0( "Objective L-BFGS-B: ", output_optim$value, " PORT error ", output_optim2$message ) )
      if( output_optim2$convergence == 0 ) print( paste0( "Objective PORT: ", output_optim$value ) )
      if( output_optim$convergence != 0 & output_optim2$convergence != 0 ) {
        print( paste0( "L-BFGS-B error ", output_optim$message, ". PORT error ", output_optim2$message ) )
      }
    }
  }
  print( paste0( "Final B-spline ", paste( round(alpha1,3), collapse = ","), " with objective ", srllikb_fun(alpha1, yd, Wik, Zik) ) )

  hout = Xh %*% alpha1;
  Sout = exp(-XH %*% alpha1 );

  m2loglik = c(srllikb_fun(alpha1, yd, Wik, Zik),
               srllikb_grad(alpha1, yd, Wik, Zik))

  return( list(alpha1=alpha1,
               h=hout,
               S=Sout,
               m2loglik=m2loglik,
               convergence=convergence ) )
}
