
## density of a multivariate normal distribution
PDFmvnorm <- function(vY, vMu, mSigma, bLog = TRUE) {

  iN = length(vY)

  dM = as.numeric(t(vY - vMu) %*% solve(mSigma) %*% (vY - vMu))
  dPDF = -0.5*iN * log(2*pi) - 0.5*log(det(mSigma)) - 0.5 * dM

  if(!bLog) {
    dPDF = exp(dPDF)
  }
  return(dPDF)
}

## cdf of a multivariate normal distribution
CDFmvnorm <- function(vY, vMu, mSigma) {

  require(cubature)

  adaptIntegrate(dmvnorm, lowerLimit = rep(-50, length(vY)), upperLimit = vY,
                 vMu = vMu, mSigma = mSigma, bLog = FALSE, tol = 1e-03)$integral

}