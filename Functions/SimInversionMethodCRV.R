#' Algorithm for random variable simulation
#'
#' Simulates continues random variables using the Inversion Method
#'
#' @param iN Number of simulated variables. An integer.
#' @param iSeed Seed. An integer.
#' @param CDFinverse Quantile function for the simulated distribution
#' @param PDF Probability density function for the simulated distribution
#' @param sDistribution Name of distribution. A string.
#'
#' @return Vector of simulated variables
#' 
#' @examples
#' PDF           = function(vX) return(2*vX)                      
#' CDFinverse    = function(vU) return(sqrt(vU))             
#' iSeed         = 1                                          
#' iN            = 10000                           
#' sDistribution = "Triangular PDF"        
#' vX            = InversionSimCRV(CDFinverse, PDF, iSeed, iN, sDistribution) 
#' @export
InversionSimCRV = function(CDFinverse, PDF, iSeed, iN, sDistribution)   {
  
  set.seed(iSeed)
  simX = function(...) {
    U = runif(iN)
    return(CDFinverse(U))
  }

  vX = simX()
  Xplot=vX[-10 <= vX & vX <=10]                           #Scale the axis
  dRangeHist = seq(min(Xplot)-1,max(Xplot)+1,0.1)  
  dRangePDF = seq(min(Xplot)-1,max(Xplot)+1,0.01)  
  h = hist(Xplot,breaks = dRangeHist,prob = TRUE)             #Define histogram
  par(mfrow = c(1,1))
  plot(h, xlab = "x", ylab = "Density", freq = FALSE,col = "gray", main = sDistribution)
  lines(dRangePDF,PDF(dRangePDF),type = 'l', col = 'red', lwd = 1.2)
  legend("topright", legend = c("Simulated", "Theoretical"), lty = c(1, 1),
         lwd = c(5, 1.2), col = c("gray", "red"), cex = 0.6)
  return(vX)
}

