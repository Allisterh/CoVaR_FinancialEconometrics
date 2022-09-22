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

mData <- pct_log_returns(mData)
# Plot the returns of all tickers 
plot(mData, type = 'l', lwd = 0.75)
# This looks very cool!

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
autocorrelation <- function(returns, exponent = 1) {
  #' returns the autocorrelation coefficient for lag = 1
  # `[2]` due to acf calculates both lag=0 and lag=1
  return(
    acf(
      returns^exponent,
      lag = 1,
      na.action = na.pass,
      plot = FALSE
    )$acf[2]
  )
}

library(moments)

DescStat[, "mean"] <- colMeans(mData, na.rm = TRUE)
DescStat[, "median"] <- apply(mData, 2, median, na.rm = TRUE)
DescStat[, "variance"] <- apply(mData, 2, var, na.rm = TRUE)
DescStat[, "kurtosis"] <- apply(mData, 2, kurtosis, na.rm = TRUE)
DescStat[, "skewness"] <- apply(mData, 2, skewness, na.rm = TRUE)
DescStat[, "rho"] <- apply(mData, 2, autocorrelation)
DescStat[, "rho2"] <- apply(mData, 2, autocorrelation, exponent = 2)

# print descriptive statistics
DescStat