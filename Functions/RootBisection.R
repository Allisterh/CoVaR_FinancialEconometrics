#' Algorithm for univariate root funding
#'
#' Solves numerically for the root of univariate functions using the Bisection Method
#'
#' @param dX.l Initial value for left bound of search area
#' @param dX.r Initial value for right bound of search area
#' @param dTol Convergence criterion. Default 1e-9.
#' @param f Function
#' @param ... An additional argument passed to f
#'
#' @return Function root (a double)
#'
#' @examples
#' f = function(dX) 0.1 - (dX*(dX + 1.0)^15)/((dX + 1.0)^15 + 1)   
#' dX.l        = -10                                                   
#' dX.r        = 10                                                   
#' dTol        = 1e-9                                              
#' dRootBi     = bisection(f, dX.l, dX.r, dTol)   
#' @export
bisection <- function(f, dX.l, dX.r, dTol = 1e-9,...) {

  if (dX.l >= dX.r) {
    cat("error: x.l >= x.r \n")
    return(NULL)
  }
  f.l <- f(dX.l, ...)
  f.r <- f(dX.r, ...)
  if (f.l == 0) {
    return(dX.l)
  } else if (f.r == 0) {
    return(dX.r)
  } else if (f.l*f.r > 0) { 
    cat("error: f(x.l)*f(x.r) > 0 \n")
    return(NULL)
  }
  
  iter <- 0
  while ((dX.r - dX.l) > dTol) {
    dX.m <- (dX.l + dX.r)/2
    f.m <- f(dX.m, ...)
    if (f.m == 0) {
      return(dX.m)
    } else if (f.l*f.m < 0) {   
      dX.r <- dX.m
      f.r <- f.m
    } else {
      dX.l <- dX.m
      f.l <- f.m
    }
    iter <- iter + 1
    cat("at iteration", iter, "the root lies between", dX.l, "and", dX.r, "\n")
  }
  
  RootPlot = function(...) {
    dRoot = (dX.l + dX.r)/2
    vX = seq(dRoot - 5,dRoot + 5,1e-2)
    plot(vX, f(vX, ...), type = "l",main = "Root search: Bisection method")
    abline(h = 0, col = "red")
    abline(v = dRoot, col = "blue", lty = 2)
  }
  
  RootPlot()
  return((dX.l + dX.r)/2)
  
}