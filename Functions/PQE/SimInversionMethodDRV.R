#' Algorithm for random variable simulation
#'
#' Simulates discrete random variables using the Inversion Method
#'
#' @param iN Number of simulated variables
#' @param iSeed Seed 
#' @param CDFinverse Quantile function for the simulated distribution
#' @param sDistribution Name of distribution; A string
#'
#' @return Vector of simulated variables
#' 
#' @examples
#' CDFinverse    = function(U) return( ifelse(U <= 1-0.5, 0, 1) )   
#' iSeed         = 1                                      
#' iN            = 10000                 
#' sDistribution = "Triangular PDF"                    
#' vX            = InversionSimDRV(CDFinverse, iSeed, iN, sDistribution)
#' @export

InversionSimDRV = function(CDFinverse, iSeed, iN, sDistribution)   {
  
  set.seed(iSeed)
  simX = function(...) {
    U = runif(iN)
    return(CDFinverse(U))
  }
  
  vX = simX()
  Xplot=vX[-10 <= vX & vX <=10]                              #Scale the axis
  dRange = seq(min(Xplot)-1,max(Xplot)+1,0.1)             #Grid to plot
  H = hist(Xplot,breaks = dRange,freq = FALSE)             #Define histogram
  H$density = H$counts / iN
  plot(H, xlab = "x", ylab = "Density", freq = FALSE,col = "gray", main = sDistribution)
  return(vX)
}
