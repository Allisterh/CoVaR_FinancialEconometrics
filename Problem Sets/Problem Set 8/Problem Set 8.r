rm(list=ls())
setwd("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Problem Sets/Problem Set 8")

################################################################################
### Problem 2                                                                ###
################################################################################

### Point i) Write a function that selects omega_t+1 by minimizing the one step 
#           ahead portfolio variance.
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

### Point ii) Write a function that selects omega_t+1 by minimizing the one step 
#             ahead portfolio variance subject to a minimum expected return of k%.

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

### Test the portfolio selection functions

# Assume these mean and variance
mSigma = matrix(c(5, -4,
                  -4, 15),2, byrow = TRUE)
vMu = c(2, 5)

# different levels of expected returns
vK = seq(0.01, 10.0, 0.01)

mFrontier = t(sapply(vK, function(dK) {
  set = EfficientSet(mSigma, vMu, dK)
  c(set$mean, set$variance)
}))

# plot of the Frontier
plot(mFrontier[,2], mFrontier[,1], type = "l", ylab = "variance", xlab = "mean")


if(!require(ggplot2)){install.packages('ggplot2')}
if(!require(ggthemes)){install.packages('ggthemes')}

df <- data.frame(
    x = mFrontier[,2],
    y = mFrontier[,1]
    )
# Plot the efficient frontier using ggplot2
ggplot(df, aes(x, y)) +
  geom_point() +
  labs(x = "Mean", y = "Variance") +
  theme_economist() +
  theme(legend.position = "none")
ggsave("./img/EfficientFrontier.pdf")


### Point iii) Write a function to estimate the CCC model in the Gaussian
            #  bivariate case assuming that the univariate models are 
            #  AR(1)–GARCH(1,1). The function should also return the one 
            #  step ahead prediction of the conditional mean and covariance matrix. 
            #  You can use the rugarch package to estimate the univariate models.

####################################################################################################
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

### Point iv) Write a function to estimate the DCC model in the Gaussian
           #  bivariate case assuming that the univariate models are 
           #  AR(1)–GARCH(1,1). The function should also return the one 
           #  step ahead prediction of the conditional mean and covariance matrix. 
           #  You can use the rugarch package to estimate the univariate models.


source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/DCC.r")

################################################################################
### Problem 3                                                                ###
################################################################################

