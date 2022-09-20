#' Algorithm for function optimization
#'
#' Solves numerically for local optimum of multivariate functions using the Newton Raphson method
#'
#' @param f Multivariate function to optimize
#' @param grad Gradient vector of f. May be supplied numerically through the grad(vX,f) function using the numDeriv package.
#' @param hess Hessian matrix of f. May be supplied numerically through the hessian(vX,f) function using the numDeriv package.
#' @param dTol Convergence criterion. Default 1e-9.
#' @param max.iter Maximum number of iterations. Default 1000.
#' @param f Function to optimize
#' @param vX0 Vector of initial values
#' @param ... An additional argument passed to f, grad or hess
#' @return Local function optimum (a numeric vector)
#' 
#' @examples
#' f           = function(vX)  (1-vX[1])^2+100*(vX[2]-vX[1]^2)^2             
#' grad        = function(vX)  c(  2 * (vX[1]-1) - 400 * ( vX[2]-vX[1]^2) * vX[1],   200 * ( vX[2]-vX[1]^2)  )
#' hess        = function(vX)  matrix(c(2 - 400 * (vX[2] - vX[1]^2) + 800 * vX[1]^2, -400 * vX[1], -400 * vX[1], 200),nrow = 2, ncol = 2, byrow = TRUE)
#' vX0         = c(-3,2)                                                           
#' dTol        = 1e-9                                                                
#' max.iter    = 1000                                                              
#' vMaxNR      = NRmultivariateOpt(f, grad, hess, vX0, dTol, max.iter) 
#' @export
NRmultivariateOpt = function(f, grad, hess, vX0 = c(-3,2), dTol = 1e-9, max.iter = 1000, ...)  {
  vX = vX0#Initial guess
  vGrad = grad(vX)#Evaluate gradient in initial guess
  mHess = hess(vX)#Evaluate hessian in initial guess
  iter = 0
  while((max(abs(vGrad))) > dTol && (iter < max.iter)) {
    vX = vX - solve(mHess, vGrad)#Update maximum point using X_new =  grad(X_old) * Hess(X_old)^-1
    vGrad = grad(vX)#Evaluate gradient in new point
    mHess = hess(vX)#Evaluate hessian in new point
    iter = iter + 1
    cat("At iteration", iter, "x =", vX[1],"y =",vX[2],"\n")
  }
  if (iter == max.iter) {
    cat("Algorithm failed to converge\n")
    return(NULL)
  } else {
    cat("Algorithm converged at iteration",iter,"\n")
    cat("The function value is",f(vX),"at x =", vX[1],"and y =",vX[2],"\n")
    return(vX)
  }
}