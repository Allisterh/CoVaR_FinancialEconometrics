# Estimate the DCC model of Engle (2002) and the Patton's (2006) model

                    ### DCC filter ###
# Filter for the dynamic correlation of the DCC model of Engle (2002)
# mEta is a iN x iT matrix of standardized residuals
# dA and dB are the two parameters of the DCC model
# mQ is the unconditional correlation matrix compute as cor(mEta)

DCCFilter <- function(mEta, dA, dB, mQ) {

  iN = ncol(mEta)
  iT = nrow(mEta)

  # Initialize the array for the correlations, R_t from the slides.
  aR = array(0, dim = c(iN, iN, iT + 1))
  # Initialize the array for the Q_t matrices
  aQ = array(0, dim = c(iN, iN, iT + 1))

  # Initialization at the unconditional correlation
  aR[, , 1] = mQ
  aQ[, , 1] = mQ

  # Compute the first likelihood contribution, slide 34
  dLLK = mEta[1, , drop = FALSE] %*% solve(aR[,, 1]) %*% t(mEta[1, , drop = FALSE]) -
    mEta[1, , drop = FALSE]%*% t(mEta[1, , drop = FALSE]) + log(det(aR[,, 1]))

  # Initiate the main loop
  for (t in 2:(iT + 1)) {
    #update the Q matrix, slide 33
    aQ[,, t] = mQ * (1 - dA - dB) + dA * t(mEta[t - 1, , drop = FALSE]) %*% mEta[t - 1, , drop = FALSE] +
      dB * aQ[,,t - 1]

    ## Compute the correlation as R_ t = Q_t_tilde^{-1/2} * Q_t * Q_t_tilde^{-1/2}
    aR[,, t] = diag(sqrt(1/diag(aQ[,, t]))) %*% aQ[,, t] %*% diag(sqrt(1/diag(aQ[,, t])))

    if (t <= iT) {
      # Augment the likelihood value, slide 35 correlation component:
      dLLK = dLLK + mEta[t, , drop = FALSE] %*% solve(aR[,, t]) %*% t(mEta[t, , drop = FALSE]) -
        mEta[t, , drop = FALSE] %*% t(mEta[t, , drop = FALSE]) + log(det(aR[,, t]))
    }
  }

  lOut = list()
  # Include the -1/2 term in the likelihood evaluation
  #see the equations in the corresponding lecture
  lOut[["dLLK"]] = -0.5 * dLLK
  lOut[["aR"]] = aR

  return(lOut)
}


                    ### DCC estimation ###
# StartingDCCPar is a vector of starting parameters for the DCC estimation
# we will use previous estimates of the DCC parameters as starting value during
# the for loop in the empirical part. This will speed up the estimation.
EstimateDCC <- function(vY1, vY2) {
  require(rugarch)

  # Model Specification
  ModelSpec = ugarchspec(mean.model = list(armaOrder = c(1, 0)))

  # Model estimation -- univariate
  Fit_1 = ugarchfit(ModelSpec, vY1)
  Fit_2 = ugarchfit(ModelSpec, vY2)

  # Standardized residuals, lecture 10, slide 31:
  # eta_t = D_t^{-1/2} * y_t
  vEta_1 = residuals(Fit_1, standardize = TRUE)
  vEta_2 = residuals(Fit_2, standardize = TRUE)

  # Unconditional correlation, lecture 10, slide 31:
  mR = cor(cbind(vEta_1, vEta_2))

  ### Model estimation -- multivariate ###
  ## Maximization of the DCC likelihood
  vPar = c(0.04, 0.9)

  # Maximize the DCC likelihood
  optimizer = solnp(vPar, fun = function(vPar, mEta, mQ) {

    Filter = DCCFilter(mEta, vPar[1], vPar[2], mQ)
    dNLLK = -as.numeric(Filter$dLLK)
    return(dNLLK)

  }, ineqfun = function(vPar, ...) {
    sum(vPar)
  }, ineqLB = 1e-4, ineqUB = 0.999,
  LB = c(1e-4, 1e-4), UB = c(0.999, 0.999),
  mEta = cbind(vEta_1, vEta_2), mQ = mR)

  vPar = optimizer$pars

  ### Prediction steps
  Forc_1 = ugarchforecast(Fit_1, n.ahead = 1)
  Forc_2 = ugarchforecast(Fit_2, n.ahead = 1)

  # One step ahead standard deviations
  vSigma1_tp1 = as.numeric(sigma(Forc_1))
  vSigma2_tp1 = as.numeric(sigma(Forc_2))
  # One step ahead mean
  vMu1_tp1 = as.numeric(fitted(Forc_1))
  vMu2_tp1 = as.numeric(fitted(Forc_2))

  ## Filter DCC
  Filter_DCC = DCCFilter(mEta = cbind(vEta_1, vEta_2), vPar[1], vPar[2], mR)

  #prediction DCC
  mR_tp1 = Filter_DCC$aR[,, length(vY1)]

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

                    ### CCC estimation ###
# 1. Define a function EstimateCCC with two arguments, vY1 and vY2.
# 2. Model specification, ModelSpec, that says that the mean model is an ARMA(1,0).
# 3. Function estimates the model for vY1 and for vY2.
# 4. Computes the standardized residuals for vY1 and for vY2.
#    and the unconditional correlation matrix mR.
# 5. Function forecasts one step ahead for vY1 and for vY2.
# 6. Computes the one step ahead mean vector vMu_tp and the covariance matrix mSigma_tp1.
# 7. The function outputs a list lOut with the estimated models, the unconditional correlation matrix,
#    the mean vector, and the covariance matrix.
EstimateCCC <- function(vY1, vY2) {
  require(rugarch)

  # Model Specification
  ModelSpec = ugarchspec(mean.model = list(armaOrder = c(1, 0)))

  # Model estimation
  Fit_1 = ugarchfit(ModelSpec, vY1)
  Fit_2 = ugarchfit(ModelSpec, vY2)

  # Standardized residuals
  # Lecture 10, slide 4: y_t = mu_t + epsilon_t <=> epsilon_t = y_t - mu_t = H_t^(-1/2) * z_t
  vEta_1 = residuals(Fit_1, standardize = TRUE)
  vEta_2 = residuals(Fit_2, standardize = TRUE)

  # Unconditional correlation
  # Lecture 10, slide 5: R_t = Corr_t-1(epsilon_t) = D_t^(-1/2) * H_t-1 * D_t^(-1/2)
  #                      D_t =  diag(H_t) 
  # Slide 27: R is estimated using the sample estimator of standardized residuals
  # eta_t = D_t^(-1/2) * y_t
  mR = cor(cbind(vEta_1, vEta_2))

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

# Two step estimation of DCC model with Student's t
# distributed shocks
Estimate_DCC_t <- function(mY) {
  # QML estimation of Sigma_t
  Fit_QML = Estimate_DCC(mY)

  # Extract standardized residuals
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


# Minimum variance portfolio
MinimumVariancePortfolio <- function(mSigma) {
  # Compute the inverse of the covariance matrix
  mSigmaInv = solve(mSigma)
  # Vector of ones
  vOnes = matrix(1, nrow = ncol(mSigma))
  #compute weights
  vOmega = mSigmaInv %*% vOnes
  #normalize weights
  vOmega = vOmega/sum(vOmega)
  return(vOmega)
}

# Efficient frontier
EfficientSet <- function(mSigma, vMu, dK) {

  mSigma_i = solve(mSigma)

  vI = matrix(1, nrow = 2)

  dA = as.numeric(t(vI) %*% mSigma_i %*% vMu)
  dB = as.numeric(t(vMu) %*% mSigma_i %*% vMu)
  dC = as.numeric(t(vI) %*% mSigma_i %*% vI)

  dD = dB*dC - dA^2

  vG = (dB * (mSigma_i %*% vI) - dA*(mSigma_i  %*% vMu))/dD
  vH = (dC*(mSigma_i %*% vMu) - dA * (mSigma_i  %*% vI))/dD

  vOmega = vG + vH*dK
  dSigma2_p = as.numeric(t(vOmega) %*% mSigma %*% vOmega)

  return(list("weights" = vOmega,
              "variance" = dSigma2_p,
              "mean" = dK))
}