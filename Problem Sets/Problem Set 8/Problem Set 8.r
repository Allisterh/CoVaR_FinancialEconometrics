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

setwd("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Problem Sets/Problem Set 8")

source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/DCC.r")
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/ReplaceOutlierMean.r")

################################################################################
### Problem 3                                                                ###
################################################################################

# Consider the last 2000 financial returns of Hewlett–Packard (HPQ) and
# Procter & Gamble (PG) in the dji30ret dataset available in the rugarch package.

library(rugarch)
library(Rsolnp)
library(ggplot2)
library(ggthemes)
data("dji30ret")

#full sample size
iT = 2000

vY1 = tail(dji30ret[, "HPQ"] * 100, iT)
vY2 = tail(dji30ret[, "PG"] * 100, iT)

# Plot the cumulative returns of HPQ and PG using ggplot2
df <- data.frame(
    x = seq(1, 2000),
    y = c(cumsum(vY1),
          cumsum(vY2)),
    col = rep(c("HPQ", "PG"), each = 2000)
    )
ggplot(df, aes(x, y, color = col)) + 
        geom_line() + 
        labs(x = "Time", y = "Cumulative returns") + 
        theme_economist() + 
        theme(legend.title=element_blank())
ggsave("./img/CumulativeReturns.pdf")

## we (arbitrarily) set K such that we have an annualized expected return of 7%
dK = 7/225

#length of the out of sample period
# if IF = 1000 you get the answer of the Exercise set.
# here we use IF = 100 because it is faster to run.
iF = 100

aWeights = array(NA, dim = c(iF, 2, 2, 2),
                 dimnames = list(NULL, c("CCC", "DCC"), c("MVP", "FixMean"), c("omega1", "omega2")))

## it is computationally expensive !!! try with a fixed t before.

for (t in (iT - iF):(iT - 1)) {
  # Estimate models and make prediction for time t + 1
  Fit_CCC = EstimateCCC(vY1[(t - iT + iF + 1):t], vY2[(t - iT + iF + 1):t])
  Fit_DCC = EstimateDCC(vY1[(t - iT + iF + 1):t], vY2[(t - iT + iF + 1):t])

  # Compute Tangency Portfolio
  aWeights[t - iT + iF + 1, "CCC", "FixMean", ] = EfficientSet(mSigma = Fit_CCC$mSigma_tp1, vMu = Fit_CCC$vMu_tp1, dK)$weights
  aWeights[t - iT + iF + 1, "DCC", "FixMean", ] = EfficientSet(mSigma = Fit_DCC$mSigma_tp1, vMu = Fit_DCC$vMu_tp1, dK)$weights

  # Compute MVP
  aWeights[t - iT + iF + 1, "CCC", "MVP", ] = MinimumVariancePortfolio(mSigma = Fit_CCC$mSigma_tp1)
  aWeights[t - iT + iF + 1, "DCC", "MVP", ] = MinimumVariancePortfolio(mSigma = Fit_DCC$mSigma_tp1)

  cat(paste(t, "\n"))

}

# Array of portfolio returns computed according to different strategies
# and models.
mPortfolioReturns = matrix(NA, iF, 4, dimnames = list(NULL, c("CCC-MVP", "CCC-FixMean", "DCC-MVP", "DCC-FixMean")))

mPortfolioReturns[, "CCC-MVP"] = aWeights[, "CCC", "MVP", "omega1"] * tail(vY1, iF) +
  aWeights[, "CCC", "MVP", "omega2"] * tail(vY2, iF)

mPortfolioReturns[, "DCC-MVP"] = aWeights[, "DCC", "MVP", "omega1"] * tail(vY1, iF) +
  aWeights[, "DCC", "MVP", "omega2"] * tail(vY2, iF)

mPortfolioReturns[, "CCC-FixMean"] = aWeights[, "CCC", "FixMean", "omega1"] * tail(vY1, iF) +
  aWeights[, "CCC", "FixMean", "omega2"] * tail(vY2, iF)

mPortfolioReturns[, "DCC-FixMean"] = aWeights[, "DCC", "FixMean", "omega1"] * tail(vY1, iF) +
  aWeights[, "DCC", "FixMean", "omega2"] * tail(vY2, iF)

# Porfolio statistics
mStat = matrix(NA, 5, 4, dimnames = list(c("Mean", "SD", "SR", "Kurtosis", "Skewness"),
                                         c("CCC-MVP", "CCC-TP", "DCC-MVP", "DCC-TP")))

mStat["Mean", ] = colMeans(mPortfolioReturns)
mStat["SD", ] = apply(mPortfolioReturns, 2, sd)
mStat["SR", ] = mStat["Mean", ]/mStat["SD", ]
mStat["Kurtosis", ] = apply(mPortfolioReturns, 2, function(x) mean((x - mean(x))^4)/(sd(x)^4))
mStat["Skewness", ] = apply(mPortfolioReturns, 2, function(x) mean((x - mean(x))^3)/(sd(x)^3)) 

# Plot the portfolio weights
dfW <- data.frame(
    x = rep(seq(1, iF), 4),
    y = matrix(data = c(aWeights[, "CCC", "MVP", "omega1"],
                      aWeights[, "CCC", "FixMean", "omega1"],
                      aWeights[, "DCC", "MVP", "omega1"],
                      aWeights[, "DCC", "FixMean", "omega1"]),
           nrow = iF*4, ncol = 1, byrow = TRUE),
    col = rep(c("CCC-MVP", "CCC-FixMean", "DCC-MVP", "DCC-FixMean"), each = iF)
)

# Replace outliers in dfW with mean
dfW$y[dfW$y > 1.5] = mean(dfW$y[dfW$y < 1.5])
dfW$y[dfW$y < -1.5] = mean(dfW$y[dfW$y > -1.5])
 
ggplot(dfW, aes(x, y, color = col)) + 
        geom_line() + 
        labs(x = "Time", y = "Portfolio weights") + 
        theme_economist() + 
        theme(legend.title=element_blank())
ggsave("./img/PortfolioWeights.pdf")

# Plot the portfolio returns
dfR <- data.frame(
    x = rep(seq(1, iF), 4),
    y = matrix(data = c(mPortfolioReturns[,1],
                        mPortfolioReturns[,3],
                        mPortfolioReturns[,2],
                        mPortfolioReturns[,4]),
              nrow = 4*iF, ncol = 1, byrow = TRUE),
    # append CCC-TP and DCC-TP to CCC-MVP and DCC-MVP
    col = rep(c("CCC-MVP", "CCC-FixMean", "DCC-MVP", "DCC-FixMean"), each = iF)
    )
# Replace outliers in dfW with mean
dfR[] <- lapply(dfR, replace_outlier_with_mean)

dfR$y[dfR$y > 1.5] = mean(dfR$y[dfR$y < 1.5])
dfR$y[dfR$y < -1.5] = mean(dfR$y[dfR$y > -1.5])

ggplot(dfR, aes(x, y, color = col)) + 
        geom_line() + 
        labs(x = "Time", y = "Portfolio Returns") + 
        theme_economist() + 
        theme(legend.title=element_blank())
ggsave("./img/PortfolioReturns.pdf")



vDCCMVP <- c(aWeights[, "CCC", "MVP", "omega1"])


mW <- matrix(data = c(aWeights[, "CCC", "MVP", "omega1"],
                      aWeights[, "CCC", "FixMean", "omega1"],
                      aWeights[, "DCC", "MVP", "omega1"],
                      aWeights[, "DCC", "FixMean", "omega1"]),
           nrow = iF, ncol = 4, byrow = FALSE)



vW <- append(mW[,1], mW[,2])

