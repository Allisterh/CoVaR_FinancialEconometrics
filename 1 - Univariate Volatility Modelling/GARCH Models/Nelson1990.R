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

par(mfrow = c(4, 1))

iT <- 1e4

lSim = f_SimGarch(iT = iT, dOmega = 0.0, dAlpha = dAlpha, dBeta = dBeta)
plot(lSim$vSigma2, main = "dOmega = 0", type = 'l')


## dOmega > 0

dAlpha = 0.05
dBeta = 0.955

ComputeElog(dAlpha = dAlpha, dBeta = dBeta)

iT <- 1e4

lSim = f_SimGarch(iT = iT, dOmega = 0.1, dAlpha = dAlpha, dBeta = dBeta)
plot(lSim$vSigma2, main = "dOmega > 0,a + b > 1",
     type = 'l')

## Moment condition is satisfied
# dAlpha + dBeta = 1
dAlpha = 0.05
dBeta = 0.95

ComputeElog(dAlpha = dAlpha, dBeta = dBeta)

iT <- 1e4

lSim = f_SimGarch(iT = iT, dOmega = 0.1, dAlpha = dAlpha, dBeta = dBeta)
plot(lSim$vSigma2, main = "Moment Condition Satisfied, a + b = 1", type = 'l')

# dAlpha + dBeta < 1
dAlpha = 0.0499
dBeta = 0.95

ComputeElog(dAlpha = dAlpha, dBeta = dBeta)

iT <- 1e4

lSim = f_SimGarch(iT = iT, dOmega = 0.1, dAlpha = dAlpha, dBeta = dBeta)
plot(lSim$vSigma2, main = "Moment Condition Satisfied, a + b < 1", type = 'l',)


################################# Theorem 2.d ##################################
par(mfrow = c(1, 1))
dOmega = 0.1
dBeta = 0.9
dAlpha = 0.04
iT <- 1e3
lSimTrue = f_SimGarch(iT = iT, dOmega = dOmega, dAlpha = dAlpha, dBeta = dBeta)
plot(lSimTrue$vSigma2, main = "True Process",
     type = 'l',)


f_FilterGarch <- function(vY, dOmega, dAlpha, dBeta, dSigma0) {
  ## Simulates `iT` number of GARCH(1,1) observation
  
  # placeholder for running variables
  vSigma2 <- numeric(iT)
  
  # variance initialized at unconditional value
  vSigma2[1] <- dSigma0
  # sample first observation
  vY[1] <- sqrt(vSigma2[1]) * vZ[1]
  
  for (t in 2:iT) {
    # sample volatility
    vSigma2[t] <- dOmega + dAlpha * vY[t - 1]^2 + dBeta * vSigma2[t - 1]
  }
  
  # return simulated series
  return(
    list(
      "vSigma2" = vSigma2
    )
  )
}

iN <- 10

mSigma2 <- matrix(NA,iT,iN)
vSigma0 <- seq(1.5,3,length.out = iN)

for (n in 1:iN) {
  mSigma2[, n] = f_FilterGarch(vY = lSimTrue$returns, dOmega = dOmega,
                               dBeta = dBeta, dAlpha = dAlpha, dSigma0[n])
  
}






### 
# Leopoldo's code
# This function simulates iT observations from a GARCH(1,1) model
# defined as
# y_t = sigma_t z_t, z_t ~ iid N(0,1)
# sigma_t^2 = dOmega + dAlpha y_{t-1}^2 + dBeta * sigma_{t-1}^2
SimGARCH <- function(iT, dOmega, dAlpha, dBeta) {
  
  # empty vectors of conditional variances and returns 
  vSigma2 = numeric(iT)
  vY = numeric(iT)
  
  # iid N(0,1) shocks
  vZ = rnorm(iT)
  
  #initialization at an arbitrary value
  # when the process is weakly stationary
  # a more adequate initialization is 
  # E[sigma_t^2] = dOmega/(1- dAlpha - dBeta)
  vSigma2[1] = 1
  
  # initialization of vY
  vY[1] = sqrt(vSigma2[1]) * vZ[1]
  
  # main loop
  for (t in 2:iT) {
    #update the variance process
    vSigma2[t] = dOmega + dAlpha*vY[t-1]^2 + dBeta*vSigma2[t-1]
    #compute returns
    vY[t] = sqrt(vSigma2[t]) * vZ[t]
  }
  
  lOut = list()
  lOut[["vY"]] = vY
  lOut[["vSigma2"]] = vSigma2
  
  return(lOut)
}

# This function computes the moment condition E[log(beta + alpha*eps_t^2)]
# for eps_t ~ N(0,1) using numerical integration.
ComputeElog <- function(dAlpha, dBeta) {
  
  int = integrate(function(dZ, dAlpha, dBeta) {
    
    log(dBeta + dAlpha*dZ^2)*dnorm(dZ)
    
  }, lower = -7, upper = 7,
  #lower and upper bounds are arbitrary,
  #always check robustness of your choices.
  dAlpha = dAlpha, dBeta = dBeta)
  
  dOut = int$value
  
  return(dOut)
  
}

#Examples
set.seed(420)
iT = 1000

dOmega = 0.01

# a weakly stationary GARCH
dAlpha = 0.04
dBeta = 0.95

ComputeElog(dAlpha, dBeta)

lSim = SimGARCH(iT, dOmega, dAlpha, dBeta)
par(mfrow = c(1,2))
plot.ts(lSim$vY)
plot.ts(lSim$vSigma2)

# a strongly stationary GARCH
# which is not weakly stationary

dAlpha = 0.051
dBeta = 0.95

ComputeElog(dAlpha, dBeta)

lSim = SimGARCH(iT, dOmega, dAlpha, dBeta)
par(mfrow = c(1,2))
plot.ts(lSim$vY)
plot.ts(lSim$vSigma2)

# a nonstarionary GARCH

dAlpha = 0.06
dBeta = 0.95

ComputeElog(dAlpha, dBeta)

lSim = SimGARCH(iT, dOmega, dAlpha, dBeta)
par(mfrow = c(1,2))
plot.ts(lSim$vY)
plot.ts(lSim$vSigma2)

## when T increases the process "explodes"
lSim = SimGARCH(1e5, dOmega, dAlpha, dBeta)
par(mfrow = c(1,2))
plot.ts(lSim$vY)
plot.ts(lSim$vSigma2)






