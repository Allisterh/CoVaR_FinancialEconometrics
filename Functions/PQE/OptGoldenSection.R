#' Algorithm for function optimization
#'
#' Solves numerically for local optimum of univariate functions using the Golden Section Algorithm
#'
#' @param dXl Initial value for left bound of search area
#' @param dXr Initial value for right bound of search area
#' @param dXm Initial value for middle point. Best guess for function optimum.
#' @param dTol Convergence criterion. Default 1e-9.
#' @param max.iter Maximum number of iterations. Default 1000.
#' @param dRho Golden ratio. Should not be changed.
#' @param f Function to optimize
#' @param ... An additional argument passed to f
#' @return Local function optimum (a double)
#' 
#' @examples
#' f           = function(dX) ifelse(dX == 0, 0, abs(dX)*log(abs(dX)/2)*exp(-abs(dX)))
#' dXl         = -8                                                    
#' dXm         = -5                                                      
#' dXr         = -2                                                      
#' dRho        = (1+sqrt(5))/2                                            
#' dTol        = 1e-9                                                     
#' max.iter    = 1000                                           
#' dMaxGS        = GoldSec(f, dXl,dXm,dXr, dTol, max.iter, dRho)  
#' @export
GoldSec = function(f, dXl,dXm,dXr, dTol = 1e-9, max.iter = 1000,dRho = (1+sqrt(5))/2, ...)  {
  iter = 0
  while( dXr-dXl > dTol && iter < max.iter ) {
    
    boolRlargest = dXr-dXm > dXm - dXl
    
    if (boolRlargest) {
      
      dXy = dXm + (dXr - dXm) / (1+dRho)
      if (f(dXy) >= f(dXm)) {
        dXl = dXm
        dXm = dXy
      } else {
        dXr = dXy
      }
      
    } else {
      
      dXy = dXm - (dXm - dXl) / (1+dRho)
      if (f(dXy) >= f(dXm)) {
        dXr = dXm
        dXm = dXy
      } else {
        dXl = dXy
      }
    }
    iter <- iter + 1
    cat("At iteration", iter, "middle point is:", dXm, "\n")
  }
  if (dXr-dXl > dTol) {
    cat("Algorithm failed to converge\n")
    return(NULL)
  } else {
    cat("Algorithm converged\n")
    vX = seq(dXm-5, dXm+5, 0.001)
    plot(vX, f(vX), type = "l", main = "Function optimum using Golden Section")
    abline(v = dXm, col = "blue", lty = 2)
    #legend("bottomright", legend = c("f'(x)", "Local optimum"), lty = c(1, 2),
    #       lwd = c(1,1), col = c("black", "blue", "red"), cex = 0.6)
    return(dXm)
  }
}