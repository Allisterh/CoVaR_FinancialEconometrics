#' Simulate normals
#' 
#' Function simulates random variables from the normal distribution using the Acceptance-Rejection algorithm
#' 
#' @param iN Number of simulated variables
#' @param iSeed Seed. An integer.
#' @param dC Supremum. A double. Leave out or set to NA to have the function calculate numerically.
#' @param SimPDF Probability density function for simulated variables
#' @param Envelope.Simulate Quantile function for envelope distribution used to simulate variables from the envelope distribution
#' @param EnvelopePDF Probability density function for envelope distribution
#' @param sDistribution Name of distribution. A string.
#' @param plotSupremum A logical indicating whether to plot the density ratio h = f/g and the supremum (either provided or numerically calculated). Default is FALSE.
#' @param bPlot A logical indicating whether to plot the simulated and theoretical normal density (and the supremum plot if chosen).
#' @return Vector of simulated variables
#'

#' @export
Norm.SimAccRej = function(dMu,dSigma,iN,iSeed,dC = NA,Envelope.Simulate, EnvelopePDF,
                        sDistribution, dLower = NA, dUpper = NA, plotSupremum = FALSE, bPlot = FALSE){

  N.01.density  = function(vX) (sqrt(2*pi))^-1*exp(-0.5*vX^2)
  N.muSigma.density  = function(vX) (dSigma*sqrt(2*pi))^-1*exp(-0.5*((vX-dMu)/dSigma)^2)
  
  DensityRatio = function(vX) 2*pmax( N.01.density(vX), 0 ) / EnvelopePDF(vX)
  
  if (is.na(dC)) {
    if ( any(is.na(c(dLower,dUpper))) ) {
      cat("Please either input the supremum or the upper and lower bounds for the support of the function h = f/g") 
      return(NULL)
    }
    dC = optimize(DensityRatio, maximum = TRUE, lower = dLower, upper = dUpper)$objective
  }
  
  AcceptReject.Sim = function() {
    vU <- rep(NA, iN)
    vY <- rep(NA, iN)
    vX <- rep(NA, iN)
    Unaccepted <- rep(TRUE, iN)
    
    while (any(Unaccepted)) {
      UnacceptedCount = sum(Unaccepted)
      vU = runif(UnacceptedCount)
      vY = Envelope.Simulate(iN = UnacceptedCount)
      
      vPhiY = N.01.density(vY) / EnvelopePDF(vY)/dC
      Accepted_ThisTimeLow = Unaccepted[Unaccepted] &
        ( vU <= 0.5 * vPhiY )
      Accepted_ThisTimeHigh = Unaccepted[Unaccepted] &
        ( vU > 0.5 * vPhiY & vU <= vPhiY )
      
      vX[Unaccepted][Accepted_ThisTimeLow] = -vY[Accepted_ThisTimeLow]
      vX[Unaccepted][Accepted_ThisTimeHigh] = vY[Accepted_ThisTimeHigh]
      Unaccepted[Unaccepted] = !Accepted_ThisTimeHigh & !Accepted_ThisTimeLow
    }
    return(dMu + dSigma * vX)
  }
  set.seed(iSeed)
  vX = AcceptReject.Sim() #Lyx for supremum derivation
  
  if (bPlot) {
  if (plotSupremum) par(mfrow = c(1,2)) else par(mfrow = c(1,1))
  Xplot=vX[-10 <= vX & vX <=10]    
  dRangeHist = seq(min(Xplot)-1,max(Xplot)+1,0.1)  
  dRangePDF = seq(min(Xplot)-1,max(Xplot)+1,0.01)        
  hist(Xplot, freq = FALSE,
       main = sDistribution,
       col = "gray",
       xlab = "x",
       cex.main = 1, breaks = dRangeHist)
  lines(dRangePDF,N.muSigma.density(dRangePDF),type = 'l', col = 'red', lwd = 1.2)
  legend("topright", legend = c("Simulated", "Theoretical"), lty = c(1, 1),
         lwd = c(5, 1.2), col = c("gray", "red"), cex = 0.6)
  
  if (plotSupremum) {
    supremumPlot = paste("Supremum =", round(dC, digits = 3), "(horizontal line)")
    vX = seq(min(Xplot),max(Xplot),0.0001)
    plot(vX,DensityRatio(vX), type = 'l', xlab = "x", ylab = "h = f / g",
         main = supremumPlot, cex.main = 0.85)
    abline(h = dC, col = 'red')
  }
  }
}