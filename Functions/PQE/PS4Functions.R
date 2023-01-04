#Needed to make C++ functions work in the package
#' @useDynLib PS4Package
#' @importFrom Rcpp sourceCpp
NULL

#' Algorithm for function optimization
#'
#' Solves numerically for local optimum of univariate functions using the Golden Section Algorithm
#' and plots intermediate results.
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
#' dMaxGS        = GoldPlotIntermediateRes(f, dXl,dXm,dXr, dTol, max.iter, dRho)  
#' @export
GoldPlotIntermediateRes = function(f, dXl,dXm,dXr, dTol = 1e-1, max.iter = 1000,dRho = (1+sqrt(5))/2, ...)  {
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
    iter = iter + 1
    plot(vX, sapply(vX,f), type = "l")#Using sapply for the case where f does not accept vector inputs
    abline(v = dXl, col = "black", lty = 2)
    abline(v = dXm, col = "black", lty = 2)
    abline(v = dXr, col = "black", lty = 2)
    abline(v = dXy, col = "red", lty = 1)
    title(paste("Iteration",iter))
    Sys.sleep(0.5)
  }
  if (dXr-dXl > dTol) {
    cat("Algorithm failed to converge\n")
    return(NULL)
  } else {
    cat("Algorithm converged\n")
    return(dXm)
  }
}


#############################################################################
#             11: Simulating compund Poisson variable (Monte Carlo)         #
#############################################################################
#' Monte Carlo simulation
#' 
#' Function simulates data from the dynamic poisson model for a given set of parameter values
#' 
#' @param dLambda Poisson parameter \eqn{\lambda}. A double.
#' @param dMu Log-normal dis. parameter \eqn{\mu}. A double.
#' @param dSigma2 dMu Log-normal dis. parameter \eqn{\sigma^2}. A double.
#' @param iB Number of Monte-Carlo replications. An integer.
#' @return The vector of simulated total payments (vector of compound Poisson variables)
#' @export
MonteCarloR = function(dLambda, dMu, dSigma2, iB) {
  vTotalAmount = numeric(iB)# Initialize storage vector
  for (b in 1:iB) {
    dN = rpois(1,dLambda)# Draw one Poisson distributed RV
    vX = rlnorm(dN, dMu, sqrt(dSigma2))# Draw dN log-normally distributed RV
    vTotalAmount[b] = sum(vX) #Compound Poisson RV
  }
  return(vTotalAmount)
}
