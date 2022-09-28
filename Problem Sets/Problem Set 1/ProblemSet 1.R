rm(list=ls()) 
setwd("/Users/tobiasbrammer/Library/Mobile Documents/
      com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/
      FinancialEconometrics/Problem Sets/Problem Set 1")

################################################################################
### Problem 1                                                                ###
################################################################################

# load list of tickers
vTickers <- read.csv("DJITicker.csv", sep = ";")

# running code in problem 2 without modifying data results in errors.
# some tickers have been delisted since, and we need to remove these
vTickers <- vTickers[!(vTickers$Symbol == "UTX" | vTickers$Symbol == "DWDP"), ]

################################################################################
### Problem 2                                                                ###
################################################################################

library(quantmod)

date_min <- "2021-09-16"
date_max <- "2022-09-16"

# Create dummy dataframe for first ticker in ticker-list
mData <- getSymbols(
  Symbols = vTickers$Symbol[1],
  from = date_min,
  to = date_max,
  src = "yahoo",
  symbol.lookup = TRUE,
  auto.assign = FALSE
)[, paste(vTickers$Symbol[1], "Adjusted", sep = ".")]

# remove "Adjusted" from ticker column-name
names(mData)[1] <- vTickers$Symbol[1]


# get remaining tickers in loop
for (sTicker in vTickers$Symbol[-1]) { # exclude the first ticker
  vReturns <- getSymbols(
    Symbols = c(sTicker),
    from = date_min,
    to = date_max,
    src = "yahoo",
    symbol.lookup = TRUE,
    auto.assign = FALSE
  )[, paste(sTicker, "Adjusted", sep = ".")]
  # append the column to the data-frame
  mData$sTicker <- vReturns
  # set name of column to the ticker-name
  names(mData)[ncol(mData)] <- sTicker}

# print last return observations
tail(mData)

pct_log_returns <- function(level_returns) {
  #' calculates % log returns and returns object of same shape as input
  #' Mathematical reasoning: r_{t}=(\ln(p_{t})-\ln(p_{t-1}))\cdot 100
  return(
    diff(log(level_returns)) * 100
  )
}

mReturns <- pct_log_returns(mData)
# Plot the returns of all tickers 
plot(mReturns, type = 'l', lwd = 0.75)
# This looks very cool!

rm(vReturns)
rm(sTicker)

################################################################################
### Problem 3                                                                ###
################################################################################

## create the matrix where to store descriptive statistics
DescStat = matrix(NA, 
                  nrow = nrow(vTickers),
                  ncol = 7,
                  ## the names of the rows and columns are defined with the
                  #  dimnames argument
                  dimnames = list(vTickers[, "Symbol"],
                                  c("mean", "median", "variance", "kurtosis",
                                    "skewness", "rho", "rho2"))
)

# function to calculate acf
autocorrelation <- function(vReturns, exponent = 1) {
  #' returns the autocorrelation coefficient for lag = 1
  # `[2]` due to acf calculates both lag=0 and lag=1
  return(
    acf(
      vReturns^exponent,
      lag = 1,
      na.action = na.pass,
      plot = FALSE
    )$acf[2]
  )
}

library(moments)

DescStat[, "mean"] <- colMeans(mReturns, na.rm = TRUE)
DescStat[, "median"] <- apply(mReturns, 2, median, na.rm = TRUE)
DescStat[, "variance"] <- apply(mReturns, 2, var, na.rm = TRUE)
DescStat[, "kurtosis"] <- apply(mReturns, 2, kurtosis, na.rm = TRUE)
DescStat[, "skewness"] <- apply(mReturns, 2, skewness, na.rm = TRUE)
DescStat[, "rho"] <- apply(mReturns, 2, autocorrelation)
DescStat[, "rho2"] <- apply(mReturns, 2, autocorrelation, exponent = 2)

# print descriptive statistics
DescStat

################################################################################
### Problem 4                                                                ###
################################################################################


vSP <- getSymbols(
  Symbols = "^GSPC",
  from = date_min,
  to = date_max,
  src = "yahoo",
  auto.assign = FALSE
)[, "GSPC.Adjusted"]

# Change the column name
names(vSP)[1] <- "S&P500"

# Calculate log returns
vSP_ret <- pct_log_returns(vSP)

fit_capm <- function(asset, factor) {
  # Remove NA
  asset <- asset[!is.na(asset)]
  factor <- factor[!is.na(factor)]
  
  # estimate OLS
  model_fit <- lm(asset ~ factor)
  
  # return matrix of coefficents
  mCAPM <- matrix(nrow = 3, ncol = 1)
  mCAPM[1,] <- model_fit$coefficients[1]
  mCAPM[2,] <- model_fit$coefficients[2]
  mCAPM[3,] <- (summary(model_fit)$sigma)^2
  
  return(mCAPM)
}

# calculate CAPM estimates
mCAPM <- t(apply(mReturns, 2, fit_capm, factor = vSP_ret))
colnames(mCAPM) <- c("alpha", "beta", "mse")

mCAPM

################################################################################
### Problem 5                                                                ###
################################################################################

fGARCHSim <- function(iT, dOmega, dAlpha, dBeta) {
  ## Initialize the vector of simulated returns and variances
  vY = numeric(iT)
  vSigma2 = numeric(iT)
  ## Initialize the variance at time t = 1 with its unconditional value
  vSigma2[1] = dOmega/(1.0 - dAlpha - dBeta)
  ## Sample the first observations
  vY[1] = rnorm(1, mean = 0, sd = sqrt(vSigma2[1]))
  ## Loop over iT. We start from t = 2 since t = 1 has already been sampled
  for (t in 2:iT) {
    #update the volatility
    vSigma2[t] = dOmega + dAlpha * vY[t - 1]^2 + dBeta * vSigma2[t - 1]
    #sample a new observarions
    vY[t] = rnorm(1, mean = 0, sd = sqrt(vSigma2[t]))
  }
  ## We return a list with two components: the sampled returns
  #  and the volatility.
  lOut = list()
  lOut[["Returns"]] = vY
  lOut[["Variance"]] = vSigma2
  ## Output lOut
  return(lOut)
} 

# simulate ARCH(1) by setting beta = 0
mGarchSim <- fGARCHSim(
  iT = 10000, dOmega = 0.3, dAlpha = 0.7, dBeta = 0
)

# plot simulated series
par(mfrow=c(1,2))
plot(
  mGarchSim$Variance,
  type = "l", ylab = "Conditional Variance", xlab = "Time",
  main="Conditional Variance"
)
plot(
  mGarchSim$Returns,
  type = "l", ylab = "Log Returns", xlab = "Time", main="Log Returns"
)


################################################################################
### Problem 6                                                                ###
################################################################################

## Function to evaluate the likelihood of an ARCH(1) model
ARCH.LLK <- function(vY, dOmega, dAlpha) {
  
  ## number of observations
  iT = length(vY)
  
  ## initialize the variance
  vSigma2 = numeric(iT)
  
  ## set the variance at time t = 1 at its unconditional value
  vSigma2[1] = dOmega/(1.0 - dAlpha)
  
  vSigma2[2:iT] = dOmega + dAlpha * vY[1:(iT - 1)]^2
  dLLK = sum(dnorm(vY, 0, sqrt(vSigma2), log = TRUE))
  
  # return the likelihood
  return(dLLK)
  
}

## Function to estimate an ARCH(1) model
#  We first code the objective function i.e. the negative log likelihood.
#  We specify a vector of coefficients to be estimated, the first coefficient
#  is omega, the second is alpha.

ObjectiveFunction <- function(vPar, vY) {
  
  dOmega = vPar[1]
  dAlpha = vPar[2]
  dLLK = ARCH.LLK(vY, dOmega, dAlpha)
  
  return(-dLLK)
}

## function to estimate the ARCH(1) model
EstimateARCH <- function(vY, ...) {
  
  # We set starting value for alpha equal to 0.1, and chose omega to target
  # the empirical variance by targeting the unconditional variance of the 
  # ARCH model.
  
  dAlpha = 0.3
  dOmega = var(vY) * (1.0 - dAlpha)
  
  # Precision constants (lower and upper)
  pre_l <- 1e-4 # 0.0001
  pre_u <- 1 - pre_l # 0.9999
  
  ## vector of starting parameters
  vPar = c(dOmega, dAlpha)
  
  ##optimization step
  optimizer = optim(vPar, fn = ObjectiveFunction, method = "L-BFGS-B", vY = vY,
                    ## note that we set suitable constraints on the model parameters
                    lower = c(0.00001, 0.0001), upper = c(10.0, 0.999)) 
  
  ## extract estimated parameters
  vPar = optimizer$par
  
  ## extract the likelihood computed at its maximum
  dLLK = -optimizer$value
  
  ## Compute the Average BIC
  iT = length(vY)
  BIC = (-2 * dLLK + log(iT) * length(vPar))/iT
  
  ## return a list with estimated parameters, likelihood value and BIC
  lOut = list()
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["BIC"]] = BIC
  
  return(lOut)
}

EstimateARCH(mGarchSim$Returns)

################################################################################
### Problem 7                                                                ###
################################################################################

# Set the seed to ensure replicability
set.seed(123)

## Number of replications (B = 500)
iB = 500
## Different sample sizes (T = 200, 500, 1000)
vT = c(200, 500, 1000)

## We all estimated results in an array
aCoef = array(NA, dim = c(iB, length(vT), 2, 2),
              dimnames = list(NULL, vT, c("omega", "alpha"), c("Correct", "Misspecified")))

## correctly specified model
for (iT in vT) { # Loop over the sample sized
  for (b in 1:iB) { # Loop over the replications
    # Simulate an ARCH model
    lSim = fGARCHSim(iT, dOmega = 0.3, dAlpha = 0.7, dBeta = 0.0)
    #estimate an ARCH model
    Fit = EstimateARCH(lSim[["Returns"]])
    #collect the estimated coefficients in the array
    aCoef[b, paste(iT), ,"Correct"] = Fit$vPar
    
  }
}

## misspecified model
for (iT in vT) { # Loop over the sample sized
  for (b in 1:iB) { # Loop over the replications
    
    # Simulate a GARCH model
    lSim = fGARCHSim(iT, dOmega = 0.3, dAlpha = 0.1, dBeta = 0.8)
    # Estimate an ARCH model
    Fit = EstimateARCH(lSim[["Returns"]])
    #collect the estimated coefficients in the array
    aCoef[b, paste(iT), ,"Misspecified"] = Fit$vPar
    
  }
}

## Plot the density of the estimated parameters
## we organize the plots in a 2x2 grid

## create a 2 x 2 grid of plot
par(mfrow = c(2, 2))

## plot omega in the correctly specified case
plot(density(aCoef[, "200", "omega", "Correct"]), col = "blue", main = "Omega - Correct Spec",
     xlab = "", ylab = "", lwd = 1.5, ylim = c(0, 20))
lines(density(aCoef[, "500", "omega", "Correct"]), col = "red", lwd = 1.5)
lines(density(aCoef[, "1000", "omega", "Correct"]), col = "purple", lwd = 1.5)
## add the true value
abline(v = 0.3)

## plot alpha in the correctly specified case
plot(density(aCoef[, "200", "alpha", "Correct"]), col = "blue", main = "Alpha - Correct Spec",
     xlab = "", ylab = "", lwd = 1.5, ylim = c(0, 10))
lines(density(aCoef[, "500", "alpha", "Correct"]), col = "red", lwd = 1.5)
lines(density(aCoef[, "1000", "alpha", "Correct"]), col = "purple", lwd = 1.5)
# add the true value
abline(v = 0.7)

## plot omega in the misspecified case
plot(density(aCoef[, "200", "omega", "Misspecified"]), col = "blue", main = "Omega - Misspecified Spec",
     xlab = "", ylab = "", lwd = 1.5, ylim = c(0, 5), xlim = c())
lines(density(aCoef[, "500", "omega", "Misspecified"]), col = "red", lwd = 1.5)
lines(density(aCoef[, "1000", "omega", "Misspecified"]), col = "purple", lwd = 1.5)
## add the true value
abline(v = 0.3)

## plot alpha in the misspecified specified case
plot(density(aCoef[, "200", "alpha", "Misspecified"]), col = "blue", main = "Alpha - Misspecified Spec",
     xlab = "", ylab = "", lwd = 1.5, ylim = c(0, 10))
lines(density(aCoef[, "500", "alpha", "Misspecified"]), col = "red", lwd = 1.5)
lines(density(aCoef[, "1000", "alpha", "Misspecified"]), col = "purple", lwd = 1.5)
## add the true value
abline(v = 0.1)


################################################################################
### Problem 8                                                                ###
################################################################################


# Load list of tickers
vTickers <- read.csv("DJITicker.csv", sep = ";")

# Remove delisted tickers
vTickers <- vTickers[!(vTickers$Symbol == "UTX" | vTickers$Symbol == "DWDP"), ]

library(quantmod)

date_min <- "2021-09-16"
date_max <- "2022-09-16"

lPrices = list()

for (sTicker in vTickers[, "Symbol"]) {
  ## download the price series of the ticker 
  ## Symbols = sTicker
  ## auto.assign = FALSE indicates that the result should be put in lPrices
  lPrices[[sTicker]] = getSymbols(
                                  Symbols = c(sTicker),
                                  from = date_min,
                                  to = date_max,
                                  src = "yahoo",
                                  symbol.lookup = TRUE,
                                  auto.assign = FALSE)
}

lRet = list()

pct_log_returns <- function(level_returns) {
  #' calculates % log returns and returns object of same shape as input
  #' Mathematical reasoning: r_{t}=(\ln(p_{t})-\ln(p_{t-1}))\cdot 100
  return(
    diff(log(level_returns)) * 100
  )
}

for (sTicker in vTickers[, "Symbol"]) {
  lRet[[sTicker]] = as.numeric(pct_log_returns(na.omit(lPrices[[sTicker]][, 6])))[-1]
}



## We store all the Bayesian Information Criteria in a matrix
mBIC = matrix(NA, length(lRet), 4, 
              dimnames = list(names(lRet), c("ARCH", "GARCH", "EGARCH", "GJRGARCH")))


## Estimate ARCH and store the BIC
for (sTicker in vTickers[, "Symbol"]) {
  ## Fir ARCH
  Fit = EstimateARCH(lRet[[sTicker]])
  ## store the BIC
  mBIC[sTicker, "ARCH"] = Fit$BIC
}


## install the rugarch package
# install.packages("rugarch")

#load the rugarch package
library(rugarch)

## Specify the three GARCH models
#GARCH
GARCHSpec = ugarchspec(variance.model = list(model = "sGARCH"), 
                       mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
#EGARCH
EGARCHSpec = ugarchspec(variance.model = list(model = "eGARCH"), 
                        mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
#GJRGARCH
GJRGARCHSpec = ugarchspec(variance.model = list(model = "gjrGARCH"), 
                          mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))

## Estimate the models and store the BIC
for (sTicker in vTickers[, "Symbol"]) {
  ## GARCH
  Fit_GARCH = ugarchfit(GARCHSpec, lRet[[sTicker]])
  ## EGARCH
  Fit_EGARCH = ugarchfit(EGARCHSpec, lRet[[sTicker]])
  ## GJRGARCH
  Fit_GJRGARCH = ugarchfit(GJRGARCHSpec, lRet[[sTicker]])
  
  ## extract the BIC
  mBIC[sTicker, "GARCH"] = infocriteria(Fit_GARCH)[2]
  mBIC[sTicker, "EGARCH"] = infocriteria(Fit_EGARCH)[2]
  mBIC[sTicker, "GJRGARCH"] = infocriteria(Fit_GJRGARCH)[2]
  
}

# Select the best model for each asset

vBest = apply(mBIC, 1, which.min)

vBest[] = colnames(mBIC)[vBest]

Selection = as.data.frame(vBest)
Selection


