#' Algorithm for univariate root funding
#'
#' Solves numerically for the root of univariate functions using the Secant Method
#'
#' @param dX0 Current first guess for function root
#' @param dX1 Current second guess for function root
#' @param dTol Convergence criterion. Default 1e-9.
#' @param max.iter Maximum number of iterations. Default 1000.
#' @param f Univariate function for which the root is wanted
#' @param ... An additional argument passed to f
#'
#' @return Function root (a double)
#'
#' @examples
#' f = function(dX) 0.1 - (dX*(dX + 1.0)^15)/((dX + 1.0)^15 + 1)       
#' dX0         = 0.2                                                    
#' dX1         = 0.3                                         
#' dTol        = 1e-9                                                
#' max.iter    = 1000                                                  
#' dRootSec    = SecantRoot(f, dX0, dX1, dTol, max.iter)    
#' @export
SecantRoot = function(f,dX0,dX1,dTol = 1e-9, max.iter = 1000, ...)  {
  fx = f(dX1)
  iter = 0
  while((abs(fx)) > dTol && (iter < max.iter)) {
    dOldX = dX1 #X_n-1
    dX1 = dX1 - f(dX1) * (dX1 - dX0) / ( f(dX1) - f(dX0) )
    dX0 = dOldX # Update X_n-1 = X_n
    fx = f(dX1)
    iter = iter + 1
    cat("At iteration", iter, "value of x is:", dX1, "\n")
  }
  if (abs(fx) > dTol) {
    cat("Algorithm failed to converge\n")
    return(NULL)
  } else {
    cat("Algorithm converged\n")
    RootPlot = function(...) {
      vX = seq(dX1 - 4, dX1 + 4,1e-2)
      plot(vX, f(vX, ...), type = "l",main = "Root search: Secant method")
      abline(h = 0, col = "red")
      abline(v = dX1, col = "blue", lty = 2)
    }
    RootPlot()
    return(dX1)
  }
}