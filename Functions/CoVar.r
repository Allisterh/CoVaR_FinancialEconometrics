source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/mvNorm.r")
# computes the CoVaR
ComputeCoVaR <- function(dAlpha, vMu, mSigma) {

  dVaR = qnorm(dAlpha, vMu[2], sqrt(mSigma[2,2]))

  dCoVaR = uniroot(function(dCoVaR, dVaR, dAlpha) {

    PDFmvnorm(c(dCoVaR, dVaR), vMu, mSigma) - dAlpha^2

  }, lower = -10, upper = 10, dVaR = dVaR, dAlpha = dAlpha, extendInt = "yes", maxiter = 10000)$root

  return(dCoVaR)

}