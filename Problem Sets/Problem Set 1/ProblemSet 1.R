rm(list=ls()) 
setwd("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Problem Sets/Problem Set 1")

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
                  ## the names of the rows and columns are defined with the dimnames argument
                  dimnames = list(vTickers[, "Symbol"],
                                  c("mean", "median", "variance", "kurtosis", "skewness", "rho", "rho2"))
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
