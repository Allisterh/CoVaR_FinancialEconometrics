
## GARCH filter
GARCHFilter <- function(vY, dOmega, dAlpha, dBeta) {

  iT = length(vY)
  vSigma2 = numeric(iT)
  vLLK = numeric(iT)

  vSigma2[1] = dOmega/(1.0 - dAlpha - dBeta)
  vLLK[1] = dnorm(vY[1], 0, sqrt(vSigma2[1]), log = TRUE)

  for(t in 2:iT) {

    vSigma2[t] = dOmega + dAlpha * vY[t-1]^2 + dBeta * vSigma2[t-1]
    vLLK[t] = dnorm(vY[t], 0, sqrt(vSigma2[t]), log = TRUE)

  }

  lOut = list()
  lOut[["vY"]] = vY
  lOut[["vSigma2"]] = vSigma2
  lOut[["vLLK"]] = vLLK
  lOut[["dLLK"]] = sum(vLLK)

  return(lOut)

}

## This function estimates a GARCH model
Estimate_GARCH <- function(vY) {

  require(Rsolnp)

  vPar = c("omega" = var(vY)*0.05,
           "alpha" = 0.05,
           "beta" = 0.9)

  optimizer = solnp(vPar, function(vPar, vY) {

    dNLLK = -GARCHFilter(vY,
                dOmega = vPar["omega"],
                dAlpha = vPar["alpha"],
                dBeta = vPar["beta"])$dLLK

    return(dNLLK)

  }, vY = vY, LB = c(0.0001, 0.0001, 0.0001), UB = c(5, 0.99, 0.99),
  ineqfun = function(vPar, ...) {
    vPar["alpha"] + vPar["beta"]
  }, ineqLB = 0.0001, ineqUB = 0.999)

  vPar = optimizer$pars
  dLLK = -tail(optimizer$values, 1)

  Filter = GARCHFilter(vY,
                       dOmega = vPar["omega"],
                       dAlpha = vPar["alpha"],
                       dBeta = vPar["beta"])

  vZ = vY/sqrt(Filter$vSigma2)

  lOut = list()
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["Filter"]] = Filter
  lOut[["vZ"]] = vZ

  return(lOut)

}