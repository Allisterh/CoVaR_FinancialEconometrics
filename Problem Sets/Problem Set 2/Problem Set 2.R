
################################################################################
### Problem 1                                                                ###
################################################################################

# Solved in .tex/.pdf file

################################################################################
### Problem 2                                                                ###
################################################################################

# Write a code to estimate the SV model of slide 18 of Lecture 6 using the 24
# moments you have derived in the previous point.
# Simulate T = 1000 observations from the SV model with omega = 0, phi = 0.9,
# and sigma2_eta = 0.25. Set the seed to 123. Estimate the model parameters
# with the GMM estimator you derived.


## Point a) 
rm(list=ls()) 
setwd("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Problem Sets/Problem Set 2")


# See the file GMMEstimation_SV.R
source("GMMEstimation_SV.R")

## function to simulate from a SV model
SV_sim <- function(iT, dOmega, dPhi, dSigma2) {
  
  vU = rnorm(iT)
  vW = numeric(iT)
  
  # initialize sampling from the unconditional distribution
  vW[1] = rnorm(1, dOmega/(1-dPhi), sqrt(dSigma2/(1 - dPhi^2)))
  
  ## loop over t
  for (t in 2:iT) {
    vW[t] = dOmega + dPhi * vW[t-1] + sqrt(dSigma2) * rnorm(1)
  }
  
  # compute returns
  vSigma = exp(0.5 * vW);
  vR = vSigma * vU
  
  # output the results
  lOut = list()
  lOut[["vR"]] = vR
  lOut[["vLogVol"]] = log(vSigma)
  lOut[["vSigma"]] = vSigma
  
  return(lOut)
}

# Point c)

set.seed(123)

iT = 1000
dOmega = 0
dPhi = 0.9
dSigma2 = 0.25

lSim = SV_sim(iT, dOmega, dPhi, dSigma2) 

vR = lSim$vR

## plot the returns
plot.ts(vR)
## plot the variances
plot.ts(lSim$vSigma^2)

Fit_GMM = GMM_Estimator(vR)

Fit_GMM$par

########################################################
#                         EX 3                         #
########################################################

library(quantmod)

## download the series
mPrice = getSymbols("^GSPC", from = "2005-01-01", to = "2018-01-01", auto.assign = FALSE)

## compute log returns
vR = as.numeric((diff(log(mPrice[, 6]))*100)[-1])

## replace zeros with the empirical mean
vR[vR == 0] = mean(vR)

## estimate the model by GMM (may take a bit of time)
Fit_GMM = GMM_Estimator(vR)

## see estimated parameters
Fit_GMM$par

