
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/ForcingVariable.r")
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/GammaPDF.r")

## Estimate the static model
Estimate_StaticGamma <- function(vY) {

  vPar = c("mu" = mean(vY),
           "a" = 160)

  optimizer = optim(vPar, function(vPar, vY) {

    dNLLK = -sum(GammaPDF(vY, dMu = vPar["mu"], dA = vPar["a"], bLog = TRUE))

    if(!is.finite(dNLLK)) {
      dNLLK = 1e10
    }

    return(dNLLK)

  }, vY = vY, method = "L-BFGS-B", lower = c(0.001, 0.001), upper = c(50, 300))

  vPar = optimizer$par
  dLLK = -optimizer$value

  iK = length(vPar)
  iT = length(vY)

  BIC = iK * log(iT) - 2*dLLK

  lOut = list()
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["optimizer"]] = optimizer
  lOut[["BIC"]] = BIC

  return(lOut)

}