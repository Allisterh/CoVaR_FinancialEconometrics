#' Algorithm for univariate root funding
#'
#' Solves numerically for the root of univariate functions using the Newton Raphson Method
#'
#' @param dX0 Initial value for right bound of search area
#' @param dTol Convergence criterion. Default 1e-9.
#' @param max.iter Maximum number of iterations. Default 1000.
#' @param f Univariate function for which the root is wanted
#' @param f_prime First derivative of f
#' @param ... An additional argument passed to f
#'
#' @return Function root (a double)
#'
#' @examples
#' f           = function(dX) 0.1 - (dX*(dX + 1.0)^15)/((dX + 1.0)^15 + 1) 
#' f_prime     = function(dX) numDeriv::grad(f,dX)                                 
#' dX0         = 0                                                  
#' dTol        = 1e-9                                                   
#' max.iter    = 1000                                                   
#' dRootNR     = NRrootUnivariate(f, f_prime, dX0, dTol, max.iter)
#' @export
NRrootUnivariate = function(f, f_prime, dX0, dTol = 1e-9, max.iter = 1000, ...)  {
  dX <- dX0
  fx <- f(dX, ...)
  iter <- 0
  while((abs(fx)) > dTol && (iter < max.iter)) {
    dX <- dX - f(dX, ...)/f_prime(dX, ...)
    fx <- f(dX, ...)
    iter <- iter + 1
    cat("At iteration", iter, "value of x is:", dX, "\n")
  }
  if (abs(fx) > dTol) {
    cat("Algorithm failed to converge\n")
    return(NULL)
  } else {
    cat("Algorithm converged\n")
    RootPlot = function(...) {
      vX = seq(dX-1,dX+1,1e-2)
      par(mfrow = c(1,1))
      plot(vX, f(vX, ...), type = "l", main = "Root search: Newton Raphson method")
      abline(h = 0, col = "red")
      abline(v = dX, col = "blue", lty = 2)
    }
    RootPlot()
    return(dX)
  }
}