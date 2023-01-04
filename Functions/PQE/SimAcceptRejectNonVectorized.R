#' Acceptance-Rejection algorithm for random variable simulation
#'
#' Simulates continues random variables using the Acceptance-Rejection Method.
#' The function is NOT vectorized. The function can calculate the supremum of the density ratio numerically.
#' The function plots the histogram for the simulated variables along with superimposed theoretical density.
#' The function can further plot the ratio of densities and indicate the supremum.
#'
#' @param iN Number of simulated variables
#' @param iSeed Seed. An integer.
#' @param dC Supremum. A double. Leave out or set to NA to have the function calculate numerically.
#' @param SimPDF Probability density function for simulated variables
#' @param Envelope.Simulate Function that simulates variables from the envelope distribution
#' @param EnvelopePDF Probability density function for envelope distribution
#' @param sDistribution Name of distribution. A string.
#' @param plotSupremum A logical indicating whether to plot the density ratio h = f/g and the supremum (either provided or numerically calculated). Default is FALSE.
#' @param dLower Lower bound for PDF support. Required for numerical approximation of supremum only.
#' @param dUpper Upper bound for PDF support. Required for numerical approximation of supremum only.
#' @return Vector of simulated variables
#'
#' @examples
#' EnvelopePDF       = function(dX) ifelse(dX >= 0 & dX <= 2, 1/2, 0)
#' Envelope.Simulate = function() return(runif(1,0,2)) 
#' SimPDF            = function(dX) ifelse(dX >= 0 & dX <= 1, dX, ifelse(dX>1 & dX<=2,2-dX,0))
#' iSeed             = 1  
#' iN                = 10000
#' dC                = 2  
#' sDistribution     = "Bernoulli PDF"
#' vX                = AccRejSimCRV(iN,iSeed,dC,SimPDF,Envelope.Simulate,EnvelopePDF,sDistribution) 
#' @export
AccRejSimCRVNonVectorized = function(iN,iSeed,dC = NA,SimPDF,Envelope.Simulate, EnvelopePDF,
                                     sDistribution, dLower = NA, dUpper = NA, plotSupremum = FALSE){
  
  if (plotSupremum | is.na(dC)) DensityRatio = function(dX) SimPDF(dX) / EnvelopePDF(dX)
  if (is.na(dC)) {
    if ( any(is.na(c(dLower,dUpper))) ) {
      cat("Please either input the supremum or the upper and lower bounds for the support of the function h = f/g") 
      return(NULL)
    }
    dC = optimize(DensityRatio, maximum = TRUE, lower = dLower, upper = dUpper)$objective
  }
  
  AcceptReject.Sim = function() {
    
    vX <- rep(NA, iN) # Pre-allocate memory for output vector
    
    iRep = OutIdx = 0 #Total replications (iRep) and index of output value to fill (OutIdx)
    
    while ( any( is.na(vX) ) ) {
      iRep = iRep + 1 # Lead total iterations by one
      dU = runif(1,0,1) # Simulate one U(0,1) variable
      dY = Envelope.Simulate() # Simulate one variable from envelope distribution
      
      if (dU <= SimPDF(dY)/( dC * EnvelopePDF(dY) ) ) {
        OutIdx = OutIdx + 1 # Lead index number of output vector
        vX[OutIdx] = dY
      }
    }
    cat("Number of rejections =", iRep - iN, "\n")
    cat("Acceptance rate =", iN / iRep, "\n")
    cat("Supremum =", dC, "\n")
    return(vX)
  }
  set.seed(iSeed)
  vX = AcceptReject.Sim() #Lyx for supremum derivation
  
  if (plotSupremum) par(mfrow = c(1,2)) else par(mfrow = c(1,1))
  Xplot=vX[-10 <= vX & vX <=10]    
  dRangeHist = seq(min(Xplot)-1,max(Xplot)+1,0.1)  
  dRangePDF = seq(min(Xplot)-1,max(Xplot)+1,0.01)      
  hist(Xplot, freq = FALSE,
       main = sDistribution,
       col = "gray",
       xlab = "x",
       cex.main = 1, breaks = dRangeHist)
  lines(dRangePDF,SimPDF(dRangePDF),type = 'l', col = 'red', lwd = 1.2)
  legend("topright", legend = c("Simulated", "Theoretical"), lty = c(1, 1),
         lwd = c(5, 1.2), col = c("gray", "red"), cex = 0.6)
  
  if (plotSupremum) {
    supremumPlot = paste("Supremum =", round(dC, digits = 3), "(horizontal line)")
    vX = seq(min(Xplot),max(Xplot),0.01)
    plot(vX,DensityRatio(vX), type = 'l', xlab = "x", ylab = "h = f / g",
         main = supremumPlot, cex.main = 0.85)
    abline(h = dC, col = 'red')
  }
}
