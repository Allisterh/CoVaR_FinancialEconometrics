source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/ForcingVariable.r")
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/GammaPDF.r")

# The filter for the MEM model
Filter_MEM <- function(vY, dKappa, dEta, dPhi, dA) {

  iT = length(vY)

  vMu = numeric(iT)
  vLLK = numeric(iT)

  vMu[1] = dKappa/(1-dEta-dPhi)
  vLLK[1] = GammaPDF(vY[1], vMu[1], dA)

  for(t in 2:iT) {
    vMu[t] = dKappa + dEta*vY[t-1] + dPhi*vMu[t-1]
    vLLK[t] = GammaPDF(vY[t], vMu[t], dA)
  }

  lOut = list()
  lOut[["vY"]] = vY
  lOut[["vMu"]] = vMu
  lOut[["vLLK"]] = vLLK
  lOut[["dLLK"]] = sum(vLLK)

  return(lOut)

}

# This function estimates the MEM model
Estimate_MEM <- function(vY) {

  vPar = c("kappa" = mean(vY)*0.05,
           "eta" = 0.05,
           "phi" = 0.90,
           "a" = 1.0)

  optimizer = optim(vPar, function(vPar, vY) {

    dNLLK = -Filter_MEM(vY,
                        dKappa = vPar["kappa"],
                        dEta = vPar["eta"],
                        dPhi = vPar["phi"],
                             dA = vPar["a"])$dLLK

    if(!is.finite(dNLLK)) {
      dNLLK = 1e10
    }

    return(dNLLK)

  }, vY = vY, method = "L-BFGS-B",
  lower = c(0.1, 0.01, 0.01, 0.1),
  upper = c(10, 0.99, 0.99, 300))

  vPar = optimizer$par
  dLLK = -optimizer$value

  iK = length(vPar)
  iT = length(vY)

  BIC = iK * log(iT) - 2*dLLK

  Filter = Filter_MEM(vY,
                      dKappa = vPar["kappa"],
                      dEta = vPar["eta"],
                      dPhi = vPar["phi"],
                      dA = vPar["a"])

  lOut = list()
  lOut[["Filter"]] = Filter
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["BIC"]] = BIC

  return(lOut)
}