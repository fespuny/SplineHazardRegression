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
#' @importFrom stats optim nlminb constrOptim
#' @importFrom nloptr bobyqa
#' @export
#'
hazl_ker = function(yd, alpha0, Wik, Zik, Xh, XH, verbose=FALSE ){
  # wrapper function for fmincon's access to srllikb
  # MATLAB hoptions = optimset('TolFun', 1E-6, 'TolX', 1E-6, 'Display', 'none', ...
  #     'algorithm', 'interior-point', 'GradObj', 'on');
  p = ncol(Wik)
  lwr = rep(0, p);

  #OPTIMISATION WITH CONSTRAINTS
  ## WE COMPARE FOUR ALGORITHMS

  alpha1 = alpha0
  best_objective = srllikb_fun(alpha1, yd, Wik, Zik)
  convergence = c( -1,-1,-1,-1, 0 ) # 4 first will be 0 if convergence, last indicates the winning method
  applyDefaults <- function(fn, ...) { function(x) fn(x, ...) }

  ## L-BFGS-B
  output_optim = list( par = alpha0, convergence=-1, value = srllikb_fun(alpha0, yd, Wik, Zik), message="" )   #initial values
  tryCatch({
    output_optim = optim( as.numeric(alpha0), srllikb_fun, srllikb_grad, yd, Wik, Zik, method = "L-BFGS-B", lower=0)
    }, error= function(e){ output_optim$message <<- conditionMessage(e) } )
  convergence[1] = output_optim$convergence
  if( output_optim$convergence == 0 ){
    ## double-check feasibility
    if( sum(output_optim$par>=0)==p ){
      alpha1 = as.numeric(output_optim$par) #array (interpreted as column vector in matrix operations)
      best_objective = output_optim$value
      convergence[5] = 1 #first tested method did give the solution
    }
  }

  ## PORT
  output_optim2 = nlminb( start=as.numeric(alpha0), objective=srllikb_fun, gradient=srllikb_grad, hessian=NULL,yd, Wik, Zik, lower=lwr)
  convergence[2] = output_optim2$convergence
  if( output_optim2$convergence == 0 ){
    ## double-check feasibility
    if( sum(output_optim2$par>=0)==p ){
      if( output_optim2$objective < best_objective ){
        alpha1 = as.numeric(output_optim2$par) #array (interpreted as column vector in matrix operations)
        best_objective = output_optim2$objective
        convergence[5] = 2 #second tested method did give the solution
      }
    }
  }

  ## Adaptative barrier algorithm using BFGS as optimization method (Nelder-Mead did not find a solution to trouble problems)
  output_optim3 <- list( par=alpha0, convergence=-1, value = srllikb_fun(alpha0, yd, Wik, Zik), message="" )
  tryCatch({
    output_optim3 =constrOptim( theta = as.numeric(alpha0), f = applyDefaults(srllikb_fun,yd=yd,Wik=Wik,Zik=Zik), grad=applyDefaults(srllikb_grad,yd=yd,Wik=Wik,Zik=Zik),
                                       ui=diag(p), ci=rep(0,p), control = list(fnscale = 1) )
  }, error= function(e){ output_optim3$message <<- conditionMessage(e) } )
  convergence[3] = output_optim3$convergence
  if( output_optim3$convergence == 0 ){
    ## double-check feasibility
    if( sum(output_optim3$par>=0) == p ){ ## alpha >= 0 ?
      if( output_optim3$value < best_objective ){
        alpha1  = as.numeric( output_optim3$par )
        best_objective = output_optim3$value
        convergence[5] = 3 #third tested method did give the solution
      }
    }
  }

  ## BOBYQA, Bound Optimization by Quadatric Approximation
  output_optim4 <- list( par = alpha0, convergence=-1, value = srllikb_fun(alpha0, yd, Wik, Zik), message="" )

  tryCatch({
    output_optim4 =bobyqa( x0 = as.numeric(alpha0), fn = applyDefaults(srllikb_fun,yd=yd,Wik=Wik,Zik=Zik), lower=rep(0,p) )
  }, error= function(e){ output_optim4$message <<- conditionMessage(e) } )
  convergence[4] = output_optim4$convergence
  if( output_optim4$convergence > 0 ){ #the previous algorithms did possibly converge, so is this a better estimate?
    ## double-check feasibility
    if( sum(output_optim4$par>=0) == p ){ ## alpha >= 0 ?
      if( output_optim4$value < best_objective ){
        alpha1  <- output_optim4$par
        best_objective <- output_optim4$value
        convergence[5] = 4 #fourth tested method did give the solution
      }
    }
  }

  # ##report convergence properties
  # if( verbose ){
  #   if( output_optim$convergence == 0 & output_optim2$convergence== 0 ){
  #     print( paste0( "Objective L-BFGS-B: ", output_optim$value, " Objective PORT: ", output_optim2$objective) )
  #   } else {
  #     if( output_optim$convergence == 0 ) print( paste0( "Objective L-BFGS-B: ", output_optim$value, " PORT error ", output_optim2$message ) )
  #     if( output_optim2$convergence == 0 ) print( paste0( "Objective PORT: ", output_optim$value ) )
  #     if( output_optim$convergence != 0 & output_optim2$convergence != 0 ) {
  #       print( paste0( "L-BFGS-B error ", output_optim$message, ". PORT error ", output_optim2$message ) )
  #     }
  #   }
  # }

  if( convergence[1] * convergence[2] * convergence[3] != 0 & convergence[4] < 0 ){
    print( "*ERROR* none of the tested optimization methods did converge ")
    convergence[5] = 5
  }

  if( verbose ){
    print( paste0( "Final B-spline ", paste( round(alpha1,3), collapse = ","), " with objective ", srllikb_fun(alpha1, yd, Wik, Zik) ) )
  }

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
