################################################################################
#   Program to test the fit of ARCH(1), GARCH(1,1), e-GARCH(1,1),              #
#   and gjr-GARCH(1,1)                                                         #
#                                                                              #
#       Relies upon quantmod to collect returns of a list of tickers           #
#       Uses rugarch to specify and fit the different models to the returns    #
#                                                                              #
################################################################################

## install the packages
# install.packages("rugarch")
# install.packages("quantmod")

rm(list=ls()) 
setwd("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Problem Sets/Problem Set 1")

# Load list of tickers
vTickers <- read.csv("DJITicker.csv", sep = ";")

# Remove delisted tickers
vTickers <- vTickers[!(vTickers$Symbol == "UTX" | vTickers$Symbol == "DWDP"), ]

library(quantmod)

date_min <- "2021-09-16"
date_max <- "2022-09-16"

lPrices = list()

for (sTicker in vTickers[, "Symbol"]) {
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
  lRet[[sTicker]] = as.numeric(
    pct_log_returns(
      na.omit(
        lPrices[[sTicker]][, 6]
        )
      )
    )[-1] # exclude the first observation to avoid NA
}

## We store all the Bayesian Information Criteria in a matrix
mBIC = matrix(NA, length(lRet), 4, 
              dimnames = list(names(lRet), c("ARCH",
                                             "GARCH", 
                                             "EGARCH", 
                                             "GJRGARCH")
                              )
              )

#load the rugarch package
library(rugarch)

## Specify the GARCH models:
# ARCH
ARCHSpec = ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 0)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)
# GARCH
## The standard GARCH model from Bollerslev (1986)
GARCHSpec = ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)
# EGARCH - The exponential GARCH model
## Allows for asymmetric effects between positive and negative asset returns.
EGARCHSpec = ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)
# GJRGARCH - The GJR-GARCH model
## Models positive and negative shocks on the conditional variance
## asymmetrically via the use of the indicator function, I.
GJRGARCHSpec = ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)

## Estimate the models and store the BIC
for (sTicker in vTickers[, "Symbol"]) {
  ## GARCH
  Fit_ARCH = ugarchfit(ARCHSpec, lRet[[sTicker]])
  ## GARCH
  Fit_GARCH = ugarchfit(GARCHSpec, lRet[[sTicker]])
  ## EGARCH
  Fit_EGARCH = ugarchfit(EGARCHSpec, lRet[[sTicker]])
  ## GJRGARCH
  Fit_GJRGARCH = ugarchfit(GJRGARCHSpec, lRet[[sTicker]])
  
  ## extract the BIC
  mBIC[sTicker, "ARCH"] = infocriteria(Fit_ARCH)[2]
  mBIC[sTicker, "GARCH"] = infocriteria(Fit_GARCH)[2]
  mBIC[sTicker, "EGARCH"] = infocriteria(Fit_EGARCH)[2]
  mBIC[sTicker, "GJRGARCH"] = infocriteria(Fit_GJRGARCH)[2]
  
}

# Select the best model for each asset

vBest = apply(mBIC, 1, which.min)

vBest[] = colnames(mBIC)[vBest]

Selection = as.data.frame(vBest)
Selection

