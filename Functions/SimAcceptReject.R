#' Acceptance-Rejection algorithm for random variable simulation
#'
#' Simulates continues random variables using the Acceptance-Rejection Method.
#' The function is vectorized. The function can calculate the supremum of the density ratio numerically.
#' The function plots the histogram for the simulated variables along with superimposed theoretical density.
#' The function can further plot the ratio of densities and indicate the supremum.
#'
#' @param iN Number of simulated variables
#' @param iSeed Seed. An integer.
#' @param dC Supremum. A double. Leave out or set to NA to have the function calculate numerically.
#' @param SimPDF Vectorized probability density function for simulated variables
#' @param Envelope.Simulate Vectorized function that simulates variables from the envelope distribution
#' @param EnvelopePDF Vectorized probability density function for envelope distribution
#' @param sDistribution Name of distribution. A string.
#' @param plotSupremum A logical indicating whether to plot the density ratio h = f/g and the supremum (either provided or numerically calculated). Default is FALSE.
#' @param dLower Lower bound for PDF support. Required for numerical approximation of supremum only.
#' @param dUpper Upper bound for PDF support. Required for numerical approximation of supremum only.
#' @return Vector of simulated variables
#'
#' @examples
#' EnvelopePDF       = function(dX) ifelse(dX >= 0 & dX <= 2, 1/2, 0)
#' Envelope.Simulate = function(iN) return(runif(iN,0,2)) 
#' SimPDF            = function(dX) ifelse(dX >= 0 & dX <= 1, dX, ifelse(dX>1 & dX<=2,2-dX,0))
#' iSeed             = 1  
#' iN                = 10000
#' dC                = 2  
#' sDistribution     = "Bernoulli PDF"
#' vX                = AccRejSimCRV(iN,iSeed,dC,SimPDF,Envelope.Simulate,EnvelopePDF,sDistribution) 
#' @export
AccRejSimCRV = function(iN,iSeed,dC = NA,SimPDF,Envelope.Simulate, EnvelopePDF,
                        sDistribution, dLower = NA, dUpper = NA, plotSupremum = FALSE){
  
  if (plotSupremum | is.na(dC)) DensityRatio = function(vX) SimPDF(vX) / EnvelopePDF(vX)
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
    TotalUnaccepted = 0
    while (any(Unaccepted)) {
      UnacceptedCount = sum(Unaccepted)
      vU = runif(UnacceptedCount)
      vY = Envelope.Simulate(iN = UnacceptedCount)
      
      Accepted_ThisTime = Unaccepted[Unaccepted] &( vU <= ( SimPDF(vY) / EnvelopePDF(vY)/dC ) )
      
      vX[Unaccepted][Accepted_ThisTime] = vY[Accepted_ThisTime]
      Unaccepted[Unaccepted] = !Accepted_ThisTime
      TotalUnaccepted = TotalUnaccepted + sum(!Accepted_ThisTime)
    }
    cat("Number of rejections =", TotalUnaccepted, "\n")
    cat("Acceptance rate =", iN / (iN+TotalUnaccepted), "\n")
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