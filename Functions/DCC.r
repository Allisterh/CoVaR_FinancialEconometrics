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

# Two step estimation of DCC model with Student's t
# distributed shocks

Estimate_DCC_t <- function(mY) {

  # QML estimation of Sigma_t

  Fit_QML = Estimate_DCC(mY)

  # extract standardized residuals
  mZ = Fit_QML$mZ

  iT = nrow(mZ)
  iP = ncol(mZ)

  # Estimate a Multivariate Student's t
  # on the standardized residuals

  optimizer = optim(5, function(dNu, mZ, iT, iP) {

    dKern = 0.0
    for (t in 1:iT) {
      dKern = dKern + log(1.0 + c(t(mZ[t, ]) %*% (mZ[t, ]))/(dNu - 2.0))
    }
    dLLK = iT * (lgamma(0.5*(dNu + iP)) - lgamma(0.5*dNu) - 0.5 * iP * log(dNu - 2.0)) - 0.5 * (dNu + iP) * dKern

    return(-dLLK)

  }, method = "L-BFGS-B", lower = 2.01, upper = 50, mZ = mZ, iT = iT, iP = iP)

  dNu = optimizer$par

  lOut = list()
  lOut[["Fit_QML"]] = Fit_QML
  lOut[["dNu"]] = dNu

  return(lOut)

}



#DCC estimation
# StartingDCCPar is a vector of starting parameters for the DCC estimation
# we will use previous estimates of the DCC parameters as starting value during
# the for loop in the empirical part. This will speed up the estimation.
EstimateDCC <- function(vY1, vY2) {
  require(rugarch)

  #Model Specification
  ModelSpec = ugarchspec(mean.model = list(armaOrder = c(1, 0)))

  #Model estimation -- univariate
  Fit_1 = ugarchfit(ModelSpec, vY1)
  Fit_2 = ugarchfit(ModelSpec, vY2)

  #standardized residuas
  vZ_1 = residuals(Fit_1, standardize = TRUE)
  vZ_2 = residuals(Fit_2, standardize = TRUE)

  #unconditional correlation
  mR = cor(cbind(vZ_1, vZ_2))

  #Model estimation -- multivariate

  ## maximization of the DCC likelihood
  vPar = c(0.04, 0.9)

  #maximize the DCC likelihood
  optimizer = solnp(vPar, fun = function(vPar, mEta, mQ) {

    Filter = DCCFilter(mEta, vPar[1], vPar[2], mQ)
    dNLLK = -as.numeric(Filter$dLLK)
    return(dNLLK)

  }, ineqfun = function(vPar, ...) {
    sum(vPar)
  }, ineqLB = 1e-4, ineqUB = 0.999,
  LB = c(1e-4, 1e-4), UB = c(0.999, 0.999),
  mEta = cbind(vZ_1, vZ_2), mQ = mR)

  vPar = optimizer$pars

  ##prediction
  Forc_1 = ugarchforecast(Fit_1, n.ahead = 1)
  Forc_2 = ugarchforecast(Fit_2, n.ahead = 1)

  #one step ahead standard deviations
  vSigma1_tp1 = as.numeric(sigma(Forc_1))
  vSigma2_tp1 = as.numeric(sigma(Forc_2))
  #one step ahead mean
  vMu1_tp1 = as.numeric(fitted(Forc_1))
  vMu2_tp1 = as.numeric(fitted(Forc_2))

  #Filter DCC
  Filter_DCC = DCCFilter(mEta = cbind(vZ_1, vZ_2), vPar[1], vPar[2], mR)

  #prediction DCC
  mR_tp1 = Filter_DCC$aCor[,, length(vY1)]

  #one step ahead covariance matrix
  mSigma_tp1 = diag(c(vSigma1_tp1, vSigma2_tp1)) %*% mR_tp1 %*% diag(c(vSigma1_tp1, vSigma2_tp1))
  #one step ahead mean vector
  vMu_tp1 = c(vMu1_tp1, vMu2_tp1)

  #output
  lOut = list()

  lOut[["Fit_1"]] = Fit_1
  lOut[["Fit_2"]] = Fit_2
  lOut[["Filter_DCC"]] = Filter_DCC
  lOut[["optimizer"]] = optimizer
  lOut[["vDCCPar"]] = vPar
  lOut[["mR"]] = mR
  lOut[["mSigma_tp1"]] = mSigma_tp1
  lOut[["vMu_tp1"]] = vMu_tp1

  return(lOut)

}


# CCC estimation
EstimateCCC <- function(vY1, vY2) {
  require(rugarch)

  # Model Specification
  ModelSpec = ugarchspec(mean.model = list(armaOrder = c(1, 0)))

  # Model estimation
  Fit_1 = ugarchfit(ModelSpec, vY1)
  Fit_2 = ugarchfit(ModelSpec, vY2)

  # Standardized residuas
  vZ_1 = residuals(Fit_1, standardize = TRUE)
  vZ_2 = residuals(Fit_2, standardize = TRUE)

  # Unconditional correlation
  mR = cor(cbind(vZ_1, vZ_2))

  ## Prediction
  Forc_1 = ugarchforecast(Fit_1, n.ahead = 1)
  Forc_2 = ugarchforecast(Fit_2, n.ahead = 1)

  # One step ahead standard deviations
  vSigma1_tp1 = as.numeric(sigma(Forc_1))
  vSigma2_tp1 = as.numeric(sigma(Forc_2))
  # One step ahead mean
  vMu1_tp1 = as.numeric(fitted(Forc_1))
  vMu2_tp1 = as.numeric(fitted(Forc_2))

  # One step ahead covariance matrix
  mSigma_tp1 = diag(c(vSigma1_tp1, vSigma2_tp1)) %*% mR %*% diag(c(vSigma1_tp1, vSigma2_tp1))
  # One step ahead mean vector
  vMu_tp1 = c(vMu1_tp1, vMu2_tp1)

  # Output
  lOut = list()

  lOut[["Fit_1"]] = Fit_1
  lOut[["Fit_2"]] = Fit_2
  lOut[["mR"]] = mR
  lOut[["mSigma_tp1"]] = mSigma_tp1
  lOut[["vMu_tp1"]] = vMu_tp1

  return(lOut)

}