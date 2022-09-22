
# The density of a gamma distribution
pdf <- function(dY, dMu, dA, bLog = TRUE) {

  dPDF = -lgamma(dA) + dA*log(dA) + (dA - 1)*log(dY) - dA*log(dMu) - dA*dY/dMu

  if (!bLog) {
    dPDF = exp(dPDF)
  }

  return(dPDF)

}

# this function returns u_t
ForcingVariable <- function(dY, dMu, dA) {

  dU = sqrt(dA)*(dY/dMu - 1)
  return(dU)

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

  vLLK[1] = pdf(vY[1], vMu[1], dA)

  for(t in 2:iT) {

    vMu_tilde[t] = dOmega + dAlpha * vU[t-1] + dBeta * vMu_tilde[t - 1]
    vMu[t] = exp(vMu_tilde[t])
    vU[t] = ForcingVariable(vY[t], vMu[t], dA)

    vLLK[t] = pdf(vY[t], vMu[t], dA)

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

# The filter for the MEM model
Filter_MEM <- function(vY, dKappa, dEta, dPhi, dA) {

  iT = length(vY)

  vMu = numeric(iT)
  vLLK = numeric(iT)

  vMu[1] = dKappa/(1-dEta-dPhi)
  vLLK[1] = pdf(vY[1], vMu[1], dA)

  for(t in 2:iT) {
    vMu[t] = dKappa + dEta*vY[t-1] + dPhi*vMu[t-1]
    vLLK[t] = pdf(vY[t], vMu[t], dA)
  }

  lOut = list()
  lOut[["vY"]] = vY
  lOut[["vMu"]] = vMu
  lOut[["vLLK"]] = vLLK
  lOut[["dLLK"]] = sum(vLLK)

  return(lOut)

}

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

# load quantmod
library(quantmod)

#download the series
VIX = getSymbols("^VIX", from = "2010-01-01", to = "2019-01-01", auto.assign = FALSE)

head(VIX)

#extract the levels
vY = as.numeric(VIX$VIX.Adjusted)

#estimate the model
Fit = Estimate_GammaGAS(vY)

## estimated parameters
Fit$vPar

plot(vY, type = "l")
lines(Fit$Filter$vMu, col = "blue")

## Estimate the static model
Estimate_StaticGamma <- function(vY) {

  vPar = c("mu" = mean(vY),
           "a" = 160)

  optimizer = optim(vPar, function(vPar, vY) {

    dNLLK = -sum(pdf(vY, dMu = vPar["mu"], dA = vPar["a"], bLog = TRUE))

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

## Estimate the static model
Fit_Static = Estimate_StaticGamma(vY)

Fit$vPar
Fit_Static$vPar

Fit$BIC
Fit_Static$BIC

## Estimante MEM
Fit_MEM = Estimate_MEM(vY)

Fit_MEM$vPar

lines(Fit_MEM$Filter$vMu, col = "red")

Fit_MEM$dLLK
Fit$dLLK

vVar = Fit_MEM$Filter$vMu^2 /Fit_MEM$vPar["a"]

par(mfrow = c(3, 1))

plot(vY, type = "l", main = "The VIX Series")
plot(Fit_MEM$Filter$vMu, type = "l", main = "The Filtered Mean")
plot(vVar, type = "l", main = "The Filtered Variance")

### Question 2

# Estimate the DCC model of Engle (2002) and the Patton's (2006) model

# Filter for the dynamic correlation of the DCC model of Engle (2002)
# mEta is a iN x iT matrix of standardized residuals
# dA and dB are the two parameters of the DCC model
# mQ is the unconditional correlation matrix compute as cor(mEta)
DCCFilter <- function(mEta, dA, dB, mQ) {

  iN = ncol(mEta)
  iT = nrow(mEta)

  # initialize the array for the correlations
  aCor = array(0, dim = c(iN, iN, iT))
  # initialize the array for the Q matrices
  aQ = array(0, dim = c(iN, iN, iT))

  ## initialization at the unconditional cor
  aCor[,, 1] = mQ
  aQ[,,1] = mQ

  #Compute the first likelihood contribution
  dLLK = mEta[1, , drop = FALSE] %*% solve(aCor[,, 1]) %*% t(mEta[1, , drop = FALSE]) -
    mEta[1, , drop = FALSE]%*% t(mEta[1, , drop = FALSE]) + log(det(aCor[,, 1]))

  #main loop
  for (t in 2:iT) {
    #update the Q matrix
    aQ[,, t] = mQ * (1 - dA - dB) + dA * t(mEta[t - 1, , drop = FALSE]) %*% mEta[t - 1, , drop = FALSE] +
      dB * aQ[,,t - 1]

    ## Compute the correlation as Q_tilde^{-1/2} Q Q_tilde^{-1/2}
    aCor[,, t] = diag(sqrt(1/diag(aQ[,, t]))) %*% aQ[,, t] %*% diag(sqrt(1/diag(aQ[,, t])))

    #augment the likelihood value
    dLLK = dLLK + mEta[t, , drop = FALSE] %*% solve(aCor[,, t]) %*% t(mEta[t, , drop = FALSE]) -
      mEta[t, , drop = FALSE] %*% t(mEta[t, , drop = FALSE]) + log(det(aCor[,, t]))
  }

  lOut = list()
  #remember to include the -1/2 term in the likelihood evaluation
  #see the equations in the corresponding lecture
  lOut[["dLLK"]] = -0.5 * dLLK
  lOut[["aCor"]] = aCor

  return(lOut)
}

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

# Function to estimate the DCC model
# y_t = Sigma_t^{1/2}z_t
# Sigma_t = D_t^{1/2} R_t D_t^{1/2}
# where D_t is a diagonal matrix with
# typical element D_iit = sigma_it^2.
# we define with eta_t = D_t^{-1/2} y_t
Estimate_DCC <- function(mY) {

  ## estimate the marginal GARCH models
  require(Rsolnp)

  #list where marginal models are stored
  lFit_univariate = list()

  #estimate the univariate GARCH models
  for(n in 1:ncol(mY)) {
    lFit_univariate[[n]] = Estimate_GARCH(mY[, n])
  }

  #Compute the residuals
  mEta = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    Fit$vZ
  }))

  #####################################################

  ## maximization of the DCC likelihood

  #initial parameters
  vPar = c(0.04, 0.9)

  #unconditional correlation
  mQ = cor(mEta)

  #maximize the DCC likelihood
  optimizer = solnp(vPar, fun = function(vPar, mEta, mQ) {

    Filter = DCCFilter(mEta, vPar[1], vPar[2], mQ)
    dNLLK = -as.numeric(Filter$dLLK)
    return(dNLLK)

  }, ineqfun = function(vPar, ...) {
    sum(vPar)
  }, ineqLB = 1e-4, ineqUB = 0.999,
  LB = c(1e-4, 1e-4), UB = c(0.999, 0.999),
  mEta = mEta, mQ = mQ)

  #Extract the estimated parameters
  vPar = optimizer$pars

  #Extract the likelihood of the correlation part
  dLLK_C = -tail(optimizer$values, 1)

  #Filter the dynamic correlation using the estimated parameters
  Filter = DCCFilter(mEta, vPar[1], vPar[2], mQ)

  #extract univariate volatilities
  mSigma = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(sqrt(Fit$Filter$vSigma2))
  }))

  #extract univariate estimated parameters
  mCoef = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(Fit$vPar)
  }))

  #compute the likelihood of the volatility  part
  dLLK_V = do.call(sum, lapply(lFit_univariate, function(Fit) {
    Fit$dLLK
  }))

  #compute the total likelihood
  dLLK = dLLK_V + dLLK_C

  ## Compute z_t
  aCor = Filter[["aCor"]]
  iT = nrow(mY)

  mZ = matrix(0, iT, ncol(mY))

  for (t in 1:iT) {
    mZ[t, ] = diag(1/mSigma[t, ]) %*% solve(chol(aCor[,,t])) %*% as.numeric(mY[t, ])
  }

  BIC = log(iT) * 8 - 2 * dLLK

  lOut = list()

  #output the results
  lOut[["dLLK"]] = dLLK
  lOut[["mCoef"]] = mCoef
  lOut[["vPar"]] = vPar
  lOut[["mSigma"]] = mSigma
  lOut[["aCor"]] = aCor
  lOut[["mEta"]] = mEta
  lOut[["mZ"]] = mZ
  lOut[["BIC"]] = BIC

  return(lOut)

}

# Function to estimate the CCC model
# y_t = Sigma_t^{1/2}z_t
# Sigma_t = D_t^{1/2} R D_t^{1/2}
# where D_t is a diagonal matrix with
# typical element D_iit = sigma_it^2.
# we define with eta_t = D_t^{-1/2} y_t
Estimate_CCC <- function(mY) {

  ## estimate the marginal GARCH models
  require(Rsolnp)

  #list where marginal models are stored
  lFit_univariate = list()

  #estimate the univariate GARCH models
  for(n in 1:ncol(mY)) {
    lFit_univariate[[n]] = Estimate_GARCH(mY[, n])
  }

  #Compute the residuals
  mEta = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    Fit$vZ
  }))

  mR = cor(mEta)

  #Filter the dynamic correlation using the estimated parameters
  Filter = DCCFilter(mEta, 0, 0, mR)

  dLLK_C = Filter$dLLK

  #extract univariate volatilities
  mSigma = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(sqrt(Fit$Filter$vSigma2))
  }))

  #extract univariate estimated parameters
  mCoef = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(Fit$vPar)
  }))

  #compute the likelihood of the volatility  part
  dLLK_V = do.call(sum, lapply(lFit_univariate, function(Fit) {
    Fit$dLLK
  }))

  #compute the total likelihood
  dLLK = dLLK_V + dLLK_C

  ## Compute z_t
  iT = nrow(mY)

  mZ = matrix(0, iT, ncol(mY))

  for (t in 1:iT) {
    mZ[t, ] = diag(1/mSigma[t, ]) %*% solve(chol(mR)) %*% as.numeric(mY[t, ])
  }

  BIC = log(iT) * 8 - 2 * dLLK

  lOut = list()

  #output the results
  lOut[["dLLK"]] = dLLK
  lOut[["mCoef"]] = mCoef
  lOut[["vPar"]] = vPar
  lOut[["mSigma"]] = mSigma
  lOut[["mR"]] = mR
  lOut[["mEta"]] = mEta
  lOut[["mZ"]] = mZ
  lOut[["BIC"]] = BIC

  return(lOut)

}

## download the series

GSPC = getSymbols("^GSPC", from = "2010-01-01", to = "2019-01-01", auto.assign = FALSE)
DJI = getSymbols("^DJI", from = "2010-01-01", to = "2019-01-01", auto.assign = FALSE)

mP = merge(GSPC$GSPC.Adjusted,
           DJI$DJI.Adjusted)

vY_GSPC = diff(log(as.numeric(mP[,1]))) * 100
vY_DJI = diff(log(as.numeric(mP[,2]))) * 100

mY = cbind(vY_GSPC, vY_DJI)

Fit_DCC = Estimate_DCC(mY)
Fit_CCC = Estimate_CCC(mY)

## Compute the minimum variance portfolio
Compute_MVP <- function(mSigma) {

  iN = ncol(mSigma)
  vOnes = rep(1, iN)

  vOmega = vOnes %*% solve(mSigma)
  vOmega = vOmega/sum(vOmega)

  return(vOmega)
}

iT = nrow(mY)
iN = ncol(mY)
#array with portfolio weights

aW = array(NA, dim = c(iT, iN, 2), dimnames = list(NULL, c("GPSC", "DJI"), c("DCC", "CCC")))

for (t in 1:iT) {

  mR_DCC = Fit_DCC$aCor[,,t]
  mD_DCC = diag(Fit_DCC$mSigma[t, ])
  mSigma_DCC = mD_DCC %*% mR_DCC %*% mD_DCC

  aW[t,,"DCC"] = Compute_MVP(mSigma_DCC)

  mR_CCC = Fit_CCC$mR
  mD_CCC = diag(Fit_CCC$mSigma[t, ])
  mSigma_CCC = mD_CCC %*% mR_CCC %*% mD_CCC

  aW[t,,"CCC"] = Compute_MVP(mSigma_CCC)

}

plot.ts(aW[,1,"DCC"])
lines(aW[,1,"CCC"], col = "red")


####

## density of a multivariate normal distribution
dmvnorm <- function(vY, vMu, mSigma, bLog = TRUE) {

  iN = length(vY)

  dM = as.numeric(t(vY - vMu) %*% solve(mSigma) %*% (vY - vMu))
  dPDF = -0.5*iN * log(2*pi) - 0.5*log(det(mSigma)) - 0.5 * dM

  if(!bLog) {
    dPDF = exp(dPDF)
  }
  return(dPDF)
}

## cdf of a multivariate normal distribution
pmvnorm <- function(vY, vMu, mSigma) {

  require(cubature)

  adaptIntegrate(dmvnorm, lowerLimit = rep(-50, length(vY)), upperLimit = vY,
                 vMu = vMu, mSigma = mSigma, bLog = FALSE, tol = 1e-03)$integral

}

# computes the CoVaR
ComputeCoVaR <- function(dAlpha, vMu, mSigma) {

  dVaR = qnorm(dAlpha, vMu[2], sqrt(mSigma[2,2]))

  dCoVaR = uniroot(function(dCoVaR, dVaR, dAlpha) {

    pmvnorm(c(dCoVaR, dVaR), vMu, mSigma) - dAlpha^2

  }, lower = -10, upper = 10, dVaR = dVaR, dAlpha = dAlpha, extendInt = "yes")$root

  return(dCoVaR)

}

### Compute the CoVaR

iT = nrow(mY)

aSigma_DCC = array(NA, dim = c(iN, iN, iT))
vMu = rep(0, iN)

mCoVaR = matrix(NA, iT, 2, dimnames = list(NULL, c("alpha = 1%", "alpha = 5%")))

for (t in 1:iT) {

  mR_DCC = Fit_DCC$aCor[,,t]
  mD_DCC = diag(Fit_DCC$mSigma[t, ])
  aSigma_DCC[,,t] = mD_DCC %*% mR_DCC %*% mD_DCC

  mCoVaR[t, "alpha = 1%"] = ComputeCoVaR(dAlpha = 0.01, vMu, aSigma_DCC[,,t])
  mCoVaR[t, "alpha = 5%"] = ComputeCoVaR(dAlpha = 0.05, vMu, aSigma_DCC[,,t])

}



