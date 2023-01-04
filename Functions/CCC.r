
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/GARCH.r")

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
  # lOut[["vPar"]] = vPar
  lOut[["mSigma"]] = mSigma
  lOut[["mR"]] = mR
  lOut[["mEta"]] = mEta
  lOut[["mZ"]] = mZ
  lOut[["BIC"]] = BIC

  return(lOut)

}