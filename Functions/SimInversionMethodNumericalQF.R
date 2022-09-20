#' Algorithm for random variable simulation
#'
#' Simulates continues random variables using the Inversion Method
#'
#' @param iN Number of simulated variables
#' @param iSeed Seed 
#' @param CDF Cumulative density function for the simulated distribution
#' @param PDF Probability density function for the simulated distribution
#' @param dLbound Lower bound for the support of the simulated distribution. -Inf if unbounded.
#' @param dUbound Upper bound for the support of the simulated distribution. Inf if unbounded.
#' @param sDistribution Name of distribution; A string
#'
#' @return Vector of simulated variables
#' 
#' @examples
#' PDF           = function(vX) ifelse(vX >= 0 & vX <= 1, vX, ifelse(vX>1 & vX<=2,2-vX,0))                
#' CDF           = function(vX) ifelse(vX < 0, 0, ifelse(vX>=0 & vX<1,0.5*vX^2,ifelse(vX>=2,1,2*vX-0.5*vX^2-1))) 
#' iSeed         = 1                                       
#' iN            = 10000                                          
#' dLbound       = 0                                               
#' dUbound       = 2                                              
#' sDistribution = "Triangular PDF"                            
#' vX            = InversionSimCRVNumQF(CDF, PDF, iSeed, iN, dLbound, dUbound, sDistribution) 
#' @export
InversionSimCRVNumQF = function(CDF, PDF, iSeed, iN, dLbound, dUbound, sDistribution)   {
  
  set.seed(iSeed)
  InverseSim.AnyDis = function() {
    U = runif(iN)
    QuantileFunction = function(dY,dLower=dLbound,dUpper=dUbound) {
      optim(0, (function(dX) abs(CDF(dX) - dY)), 
            method = "L-BFGS-B", upper = dUbound,lower=dLbound)$par
    }
    return(sapply(U,QuantileFunction,dLbound,dUbound))   #Solve numerically for inverse of F for each element of U
  }
  vX = InverseSim.AnyDis()
  Xplot=vX[-10 <= vX & vX <=10]                              #Scale the axis
  dRangeHist = seq(min(Xplot)-1,max(Xplot)+1,0.1)  
  dRangePDF = seq(min(Xplot)-1,max(Xplot)+1,0.01)  
  h = hist(Xplot,breaks = dRangeHist,prob = TRUE)             #Define histogram
  
  plot(h, xlab = "x", ylab = "Density", freq = FALSE,col = "gray", main = sDistribution)
  lines(dRangePDF,PDF(dRangePDF),type = 'l', col = 'red', lwd = 1.2)
  legend("topright", legend = c("Simulated", "Theoretical"), lty = c(1, 1),
         lwd = c(5, 1.2), col = c("gray", "red"), cex = 0.6)
  return(vX)
}

