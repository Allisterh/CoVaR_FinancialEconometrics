#' Simulate normals
#' 
#' Function simulates random variables from the normal distribution using uniformly distributed variables.
#' 
#' @param dMu Mean \eqn{\mu}. A double.
#' @param dSigma Standard deviation \eqn{\sigma}. A double.
#' @param iSize Number of simulated variables. An integer.
#' @param iSeed Seed. An integer.
#' @param bPlot A logical indicating whether to plot the simulated variables along with superimposed theoretical normal density.
#' @return The vector of simulated N(mu,Sigma^2) variables.
#' @export
Norm.SimUsingUniform = function(iSize, dMu, dSigma, iSeed, bPlot = FALSE) {
  Sim = function() {
    set.seed(iSeed)
    Z = numeric(iSize)
    for (i in 1:iSize) {
      Z[i] = sum(runif(12)) - 6
    }
    return(dSigma * Z + dMu)
  }
  vX = Sim()
  if (bPlot) {
    NormalDensity = function(vX) (dSigma*sqrt(2*pi))^-1*exp(-0.5*((vX-dMu)/dSigma)^2)
    Xplot=vX[-10 <= vX & vX <=10]    
    dRangeHist = seq(min(Xplot)-0.5, max(Xplot)+0.5,0.1)
    dRangePDF = seq(min(Xplot)-0.5, max(Xplot)+0.5,0.01)  
    par(mfrow = c(1,1))
    hist(Xplot, freq = FALSE,
         main = "Theoretical and simulated normal density",
         cex.main = 0.9, col = "gray", xlab = "x", breaks = dRangeHist)
    lines(dRangePDF,NormalDensity(dRangePDF),type = 'l', col = 'red', lwd = 1.2)
    legend("topright", legend = c("Simulated N(mu,sigma2)", "Theoretical N(mu,sigma2)"), lty = c(1, 1),
           lwd = c(5, 1.2), col = c("gray", "red"), cex = 0.7)
  }
  return(vX)
  
}