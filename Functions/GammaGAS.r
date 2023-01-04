
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/ForcingVariable.r")
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/GammaPDF.r")

# This function estimates the GammaGAS model
Estimate_GammaGAS <- function(vY) {

  vPar = c("omega" = log(mean(vY))*0.05,
           "alpha" = 0.05,
           "beta" = 0.95,
           "a" = 1.0)

  optimizer = optim(vPar, function(vPar, vY) {

    dNLLK = -Filter_GammaGAS(vY,
                    dOmega = vPar["omega"],
                    dAlpha = vPar["alpha"],
                    dBeta = vPar["beta"],
                    dA = vPar["a"])$dLLK

    if(!is.finite(dNLLK)) {
      dNLLK = 1e10
    }

    return(dNLLK)

  }, vY = vY, method = "L-BFGS-B",
  lower = c(-0.5, 0.001, 0.01, 0.1),
  upper = c(0.5, 1.5, 0.999, 300))

  vPar = optimizer$par
  dLLK = -optimizer$value

  iK = length(vPar)
  iT = length(vY)

  BIC = iK * log(iT) - 2*dLLK

  Filter = Filter_GammaGAS(vY,
                  dOmega = vPar["omega"],
                  dAlpha = vPar["alpha"],
                  dBeta = vPar["beta"],
                  dA = vPar["a"])

  lOut = list()
  lOut[["Filter"]] = Filter
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["BIC"]] = BIC

  return(lOut)
}


# The filter for the model
Filter_GammaGAS <- function(vY, dOmega, dAlpha, dBeta, dA) {

  iT = length(vY)

  vMu_tilde = numeric(iT)
  vMu = numeric(iT)
  vU = numeric(iT)

  vLLK = numeric(iT)

  vMu_tilde[1] = dOmega/(1-dBeta)
  vMu[1] = exp(vMu_tilde[1])
  vU[1] = ForcingVariable(vY[1], vMu[1], dA)

  vLLK[1] = GammaPDF(vY[1], vMu[1], dA)

  for(t in 2:iT) {

    vMu_tilde[t] = dOmega + dAlpha * vU[t-1] + dBeta * vMu_tilde[t - 1]
    vMu[t] = exp(vMu_tilde[t])
    vU[t] = ForcingVariable(vY[t], vMu[t], dA)

    vLLK[t] = GammaPDF(vY[t], vMu[t], dA)

  }

  lOut = list()
  lOut[["vY"]] = vY
  lOut[["vMu_tilde"]] = vMu_tilde
  lOut[["vMu"]] = vMu
  lOut[["vU"]] = vU
  lOut[["vLLK"]] = vLLK
  lOut[["dLLK"]] = sum(vLLK)

  return(lOut)

}
