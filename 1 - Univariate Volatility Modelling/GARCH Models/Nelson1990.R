# Proof of Nelson 1990 Theorem 1 and Theorem 2

library(rugarch)
library(Rsolnp)

# replicability
set.seed(69)


################################################################################
### Simulate from a GARCH(1, 1) process                                      ###
################################################################################

f_SimGarch <- function(iT, dOmega, dAlpha, dBeta) {
  ## Simulates `iT` number of GARCH(1,1) observation
  
  # placeholder for running variables
  vY <- numeric(iT)
  vSigma2 <- numeric(iT)
  vZ <- rnorm(iT)
  
  # variance initialized at unconditional value
  #vSigma2[1] <- dOmega / (1.0 - dAlpha - dBeta)
  vSigma2[1] <- 1
  # sample first observation
  vY[1] <- sqrt(vSigma2[1]) * vZ[1]
  
  for (t in 2:iT) {
    # sample volatility
    vSigma2[t] <- dOmega + dAlpha * vY[t - 1]^2 + dBeta * vSigma2[t - 1]
    
    # sample new observation
    vY[t] <- sqrt(vSigma2[t]) * vZ[t]
  }
  
  # return simulated series
  return(
    list(
      "returns" = vY,
      "vSigma2" = vSigma2
    )
  )
}


# y_t = sigma_t z_t , z_t ~ iid N(0,1)

ComputeElog <- function(dBeta, dAlpha) {
  int = integrate(function(dZ, dBeta, dAlpha) {
      log(dBeta + dAlpha * dZ^2)*dnorm(dZ)
  },lower = -7, upper = 7,
  dAlpha = dAlpha, dBeta = dBeta)
  dOut = int$value
  return(dOut)
}

dAlpha = 0.05
dBeta = 0.955

ComputeElog(dAlpha = dAlpha, dBeta = dBeta)

iT <- 1e4

lSim = f_SimGarch(iT = iT, dOmega = 0.0, dAlpha = dAlpha, dBeta = dBeta)
plot.ts(lSim$vSigma2)

