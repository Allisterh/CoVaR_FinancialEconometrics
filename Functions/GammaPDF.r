# The density of a gamma distribution
GammaPDF <- function(dY, dMu, dA, bLog = TRUE) {
  dPDF = -lgamma(dA) + dA*log(dA) + (dA - 1)*log(dY) - dA*log(dMu) - dA*dY/dMu
  if (!bLog) {
    dPDF = exp(dPDF)
  }
  return(dPDF)
}