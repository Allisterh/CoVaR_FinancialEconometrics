rm(list = ls()) # clear environment
cat("\014")     # clear console output
setwd("C:/Users/Christian Restrup/Dropbox/7. Semester/Financial Econometrics/?velser/Data+Additional Code")
######################################################
#                 COMPUTATIONAL PART                 #
######################################################

# First point
#Leopoldo
#Leepoldo has this code for cleaning the data:
# GSPC = getSymbols("^GSPC", from = "2010-01-01", to = "2019-01-01", auto.assign = FALSE)
# DJI = getSymbols("^DJI", from = "2010-01-01", to = "2019-01-01", auto.assign = FALSE)
# 
# mP = merge(GSPC$GSPC.Adjusted,
#            DJI$DJI.Adjusted   )
# 
# vY_GSPC = diff(log(as.numeric(mP[,1]))) * 100
# vY_DJI  = diff(log(as.numeric(mP[,2]))) * 100
# 
# plot.ts(vY_GSPC)
# lines(vY_DJI, col="blue")
# 
# dY <- vY_GSPC
# mY = cbind(vY_GSPC, vY_DJI)
# # Leo out



#The density of the GED distribution
log_of_p <- function(dY, dPsi, dUpsilon, log = TRUE) {
  
  
  if (dPsi < 1e-10) {
    dPsi = 1e-10
  }
  
  dLogPdf = -((1 + 1 / dUpsilon) * log(2) + log(dPsi) + lgamma(1 + 1/dUpsilon)) - abs(dY)^dUpsilon / (2*dPsi^dUpsilon) #(Eq.8)
  
  if (log == FALSE) {
    dLogPdf = exp(dLogPdf)
  }
  
  return(dLogPdf)
}


# The filter for the GAS-GED model
# vPar is a vector of 4 elements: omega, alpha, beta, upsilon.
GASGED_Filter <- function(vPar, vY) {
  
  iT = length(vY)
  vPsi = numeric(iT)
  vPsi_tilde = numeric(iT)
  vSigma = numeric(iT)
  
  dOmega = vPar[1]
  dAlpha = vPar[2]
  dBeta = vPar[3]
  dUpsilon = vPar[4]
  
  #initialize at the unconditional value
  vPsi_tilde[1] = dOmega/(1 - dBeta)
  vPsi[1] = exp(vPsi_tilde[1])
  
  #compute the volatility
  vSigma[1] = 2^(1/dUpsilon) * sqrt((gamma(3/dUpsilon))/(gamma(1/dUpsilon))) * vPsi[1] #(Eq.23)
  
  #initialize the likelihood
  dLLK = log_of_p(vY[1], vPsi[1], dUpsilon)
  
  ##loop over the data and compute the log likelihood
  for (t in 2:iT) {
    
    vPsi_tilde[t] = dOmega + dBeta * vPsi_tilde[t - 1] + dAlpha * (abs(vY[t - 1])^dUpsilon / vPsi[t-1]^dUpsilon - 1) #(Eq.24)
    vPsi[t] = exp(vPsi_tilde[t])                                                                                     # (Eq.11)
    vSigma[t] =  2^(1/dUpsilon) * sqrt((gamma(3/dUpsilon))/(gamma(1/dUpsilon))) * vPsi[t]                            #(Eq.23)
    
    dLLK = dLLK + log_of_p(vY[t], vPsi[t], dUpsilon)
  }
  
  lOut = list()
  
  lOut[["vPsi_tilde"]] = vPsi_tilde
  lOut[["vPsi"]] = vPsi
  lOut[["vSigma"]] = vSigma
  lOut[["dLLK"]] = dLLK
  
  return(lOut)
  
}

## This function evaluate the negative log likelihood
NegLogLikelihood <- function(vPar, vY) {
  
  Filter = GASGED_Filter(vPar, vY)
  
  dNLL = -Filter[["dLLK"]]
  
  #if for some reason the output is not finite
  # use a large number to penalize the likelihood
  if (!is.finite(dNLL)) {
    dNLL = 1e5
  }
  
  return(dNLL)
}


# Function to estimate the GAS-GED model
Estimate_GASGED <- function(vY) {
  
  # Initialization. You can modify this and initialize targeting
  # the unconditional standard deviation of vY
  vPar_Starting = c(0, 0.07, 0.8, 2)
  
  # maximize the likelihood 
  optimiser = optim(vPar_Starting, NegLogLikelihood, vY = vY,
                    method = "L-BFGS-B", lower = c(-2, 1e-3, 0.001, 0.5),
                    upper = c(2, 5, 0.9999, 4)
  )
  
  vPar = optimiser$par
  dLLK = -optimiser$value
  
  # run the filter
  FilteredValues = GASGED_Filter(vPar, vY)
  
  iT = length(vY)
  
  #compute the average BIC
  ABIC = (log(iT) * 4 - dLLK)/iT
  
  lOut = list()
  
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["ABIC"]] = ABIC
  lOut[["FilteredValues"]] = FilteredValues
  
  return(lOut)
}

# Second point

### compute VaR at level alpha

## two options
# 1) inverting the cdf

## The cdf is obtained by integration of the pdf
cdf <- function(dY, dPsi, dUpsilon) {
  
  dF = integrate(log_of_p, lower = -99, upper = dY, log = FALSE,
                 dPsi = dPsi, dUpsilon = dUpsilon)$value
  
  return(dF)
}

## the quantile function is obtained by solving F(y) = alpha
quantileFun <- function(dAlpha, dPsi, dUpsilon) {
  
  dVaR = uniroot(function(dVaR, dAlpha, dPsi, dUpsilon) {
    cdf(dVaR, dPsi, dUpsilon) - dAlpha # P(Y_t < dVaR|F_{t-1}) - alpha = 0
  }, dAlpha = dAlpha, dPsi = dPsi, dUpsilon = dUpsilon,
  interval = c(-20, 20), extendInt = "yes")$root
  
  return(dVaR)
}

## vectorized version as required in the exam
quantileFun_Vec  <- function(dAlpha, vPsi, dUpsilon) {
  
  iT = length(vPsi)
  vVaR = numeric(iT)
  
  for(t in 1:iT) {
    vVaR[t] = quantileFun(dAlpha, vPsi[t], dUpsilon)
  }
  
  return(vVaR)
  
}

## second way to obtain the quantile function using qdist() from rugarch

library(rugarch)
quantileFun2_Vec <- function(dAlpha, vSigma, dUpsilon) {
  
  vVaR = qdist(distribution = "ged", p = dAlpha,
               sigma = vSigma, shape = dUpsilon)
  
  return(vVaR)
  
}


######################################################
#                 EMPIRICAL  PART                    #
######################################################

############## POINT a) ##############################

library(rugarch)

#setwd("G:/Dropbox/Teaching/Aarhus_Financial Econometrics/2018/Exam/Questions/LC/data")

GSPC = getSymbols("^GSPC", from = "2010-01-01", to = "2019-01-01", auto.assign = FALSE)
DJI = getSymbols("^DJI", from = "2010-01-01", to = "2019-01-01", auto.assign = FALSE)

mP = merge(GSPC$GSPC.Adjusted,
           DJI$DJI.Adjusted   )

GSPC = diff(log(as.numeric(mP[,1]))) * 100
DJI  = diff(log(as.numeric(mP[,2]))) * 100

plot.ts(GSPC)
lines(DJI, col="blue")

mY = cbind(GSPC, DJI)

mRet = read.csv("data.csv", sep = ";", dec = ",", header = FALSE, row.names = 1)

colnames(mRet) = c("GM", "DJI30")

head(mRet)

# i)

# Specify the three models (have a look help(ugarchspec))

GARCH_Spec = ugarchspec(variance.model = list(model = "sGARCH"),
                        mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                        distribution.model = "ged")


EGARCH_Spec = ugarchspec(variance.model = list(model = "eGARCH"),
                         mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                         distribution.model = "ged")


GJRGARCH_Spec = ugarchspec(variance.model = list(model = "gjrGARCH"),
                           mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                           distribution.model = "ged")



## Estimates for GM
Fit_GARCH_GM = ugarchfit(GARCH_Spec, mRet[, "GM"])
Fit_EGARCH_GM = ugarchfit(EGARCH_Spec, mRet[, "GM"])
Fit_GJRGARCH_GM = ugarchfit(GJRGARCH_Spec, mRet[, "GM"])


## Estimates for DJI30
Fit_GARCH_DJI30 = ugarchfit(GARCH_Spec, mRet[, "DJI30"])
Fit_EGARCH_DJI30 = ugarchfit(EGARCH_Spec, mRet[, "DJI30"])
Fit_GJRGARCH_DJI30 = ugarchfit(GJRGARCH_Spec, mRet[, "DJI30"])

# ii) Compare the filtered volatilities

iT = nrow(mRet)

mSigma_GM = matrix(0, iT, 3, dimnames = list(NULL, c("GARCH", "EGARCH", "GJRGARCH")))
mSigma_DJI30 = matrix(0, iT, 3, dimnames = list(NULL, c("GARCH", "EGARCH", "GJRGARCH")))

## if you don't remenber the sigma function have a look at help(ugarchfit)
## where you find, for instance, that to extract estimated coefficients you can use
# coef(Fit_GARCH_GM)

mSigma_GM[, "GARCH"] = sigma(Fit_GARCH_GM)
mSigma_GM[, "EGARCH"] = sigma(Fit_EGARCH_GM)
mSigma_GM[, "GJRGARCH"] = sigma(Fit_GJRGARCH_GM)

head(mSigma_GM)

mSigma_DJI30[, "GARCH"] = sigma(Fit_GARCH_DJI30)
mSigma_DJI30[, "EGARCH"] = sigma(Fit_EGARCH_DJI30)
mSigma_DJI30[, "GJRGARCH"] = sigma(Fit_GJRGARCH_DJI30)

head(mSigma_DJI30)

plot.ts(mSigma_GM, plot.type = "single", col = c("black", "red", "blue"),
        main = "Volatility comparison for GM")
legend("topleft", legend = c("GARCH", "EGARCH", "GJRGARCH"),
       col = c("black", "red", "blue"), lty = 1)

plot.ts(mSigma_DJI30, plot.type = "single", col = c("black", "red", "blue"),
        main = "Volatility comparison for DJI30")
legend("topleft", legend = c("GARCH", "EGARCH", "GJRGARCH"),
       col = c("black", "red", "blue"), lty = 1)

## iii) Select the best model using BIC
## use the inforcriteria() function
## if you don't remenber the inforcriteria function have a look at help(ugarchfit)

vBIC_GM = numeric(3)
vBIC_GM[1] = infocriteria(Fit_GARCH_GM)[2]
vBIC_GM[2] = infocriteria(Fit_EGARCH_GM)[2]
vBIC_GM[3] = infocriteria(Fit_GJRGARCH_GM)[2]
sort(vBIC_GM)
which.min(vBIC_GM)
## We select EGARCH for GM because is the one with the minimum BIC

vBIC_DJI30 = numeric(3)
vBIC_DJI30[1] = infocriteria(Fit_GARCH_DJI30)[2]
vBIC_DJI30[2] = infocriteria(Fit_EGARCH_DJI30)[2]
vBIC_DJI30[3] = infocriteria(Fit_GJRGARCH_DJI30)[2]

which.min(vBIC_DJI30)
## We select EGARCH for DJI30 because is the one with the minimum BIC

# For Point b) we use the Estimate_GASGED and quantileFun2_Vec from Ex 2.

Fit_GAS_GM = Estimate_GASGED(mRet[, "GM"])
Fit_GAS_DJI30 = Estimate_GASGED(mRet[, "DJI30"])

#i)

plot(mSigma_GM[, "GARCH"], col = "black", type = "l",
     main = "GARCH - GAS Volatility comparison for GM")
lines(Fit_GAS_GM$FilteredValues$vSigma, col = "red")
legend("topleft", legend = c("GARCH", "GAS"),
       col = c("black", "red"), lty = 1)

plot(mSigma_DJI30[, "GARCH"], col = "black", type = "l",
     main = "GARCH - GAS Volatility comparison for DJI30",
     ylim = c(0.5, 6)
)
lines(Fit_GAS_DJI30$FilteredValues$vSigma, col = "red")
legend("topleft", legend = c("GARCH", "GAS"),
       col = c("black", "red"), lty = 1)

#ii)

vAlpha = c(0.01, 0.05)

aVaR = array(0, dim = c(iT, 2, 2, 2),
             dimnames = list(NULL, c("GARCH", "GAS"), c("1%", "5%"), c("GM", "DJI30")))

aVaR[, "GARCH", "1%", "GM"] = quantileFun2_Vec(0.01, mSigma_GM[, "GARCH"], coef(Fit_GARCH_GM)[4])
aVaR[, "GARCH", "5%", "GM"] = quantileFun2_Vec(0.05, mSigma_GM[, "GARCH"], coef(Fit_GARCH_GM)[4])

aVaR[, "GAS", "1%", "GM"] = quantileFun2_Vec(0.01, Fit_GAS_GM$FilteredValues$vSigma, Fit_GAS_GM$vPar[4])
aVaR[, "GAS", "5%", "GM"] = quantileFun2_Vec(0.05, Fit_GAS_GM$FilteredValues$vSigma, Fit_GAS_GM$vPar[4])

aVaR[, "GARCH", "1%", "DJI30"] = quantileFun2_Vec(0.01, mSigma_DJI30[, "GARCH"], coef(Fit_GARCH_DJI30)[4])
aVaR[, "GARCH", "5%", "DJI30"] = quantileFun2_Vec(0.05, mSigma_DJI30[, "GARCH"], coef(Fit_GARCH_DJI30)[4])

aVaR[, "GAS", "1%", "DJI30"] = quantileFun2_Vec(0.01, Fit_GAS_DJI30$FilteredValues$vSigma, Fit_GAS_DJI30$vPar[4])
aVaR[, "GAS", "5%", "DJI30"] = quantileFun2_Vec(0.05, Fit_GAS_DJI30$FilteredValues$vSigma, Fit_GAS_DJI30$vPar[4])

plot(aVaR[, "GARCH", "1%", "GM"], type = "l", col = "black", ylim = c(-25, 0))
lines(aVaR[, "GAS", "1%", "GM"], type = "l", col = "red")
lines(aVaR[, "GARCH", "5%", "GM"], type = "l", col = "black", lty = 2)
lines(aVaR[, "GAS", "5%", "GM"], type = "l", col = "red", lty = 2)
legend("bottomleft", legend = c("GARCH 1%", "GAS 1%", "GARCH 5%", "GAS 5%"),
       col = c("black", "red", "black", "red"), lty = c(1, 1, 2, 2))

plot(aVaR[, "GARCH", "1%", "DJI30"], type = "l", col = "black", ylim = c(-25, 0))
lines(aVaR[, "GAS", "1%", "DJI30"], type = "l", col = "red")
lines(aVaR[, "GARCH", "5%", "DJI30"], type = "l", col = "black", lty = 2)
lines(aVaR[, "GAS", "5%", "DJI30"], type = "l", col = "red", lty = 2)
legend("bottomleft", legend = c("GARCH 1%", "GAS 1%", "GARCH 5%", "GAS 5%"),
       col = c("black", "red", "black", "red"), lty = c(1, 1, 2, 2))

## iii) Select the best model between GAS-GED and GARCH usign BIC

vBIC_GM = numeric(2)
vBIC_GM[1] = infocriteria(Fit_GARCH_GM)[2]
vBIC_GM[2] = Fit_GAS_GM$ABIC

which.min(vBIC_GM)

vBIC_DJI30 = numeric(2)
vBIC_DJI30[1] = infocriteria(Fit_GARCH_DJI30)[2]
vBIC_DJI30[2] = Fit_GAS_DJI30$ABIC

which.min(vBIC_DJI30)

# We select the GAS model for both series

############## POINT b) ##############################

# Point c)

iT = nrow(mRet)

aSigma_GAS = array(0, dim = c(2, 2, iT))

## Compute the correlation matrix

mR = cor(mRet)

# the variance of the two series
aSigma_GAS[1,1, ] = Fit_GAS_GM$FilteredValues$vSigma^2
aSigma_GAS[2,2, ] = Fit_GAS_DJI30$FilteredValues$vSigma^2

aSigma_GAS[1,2, ] = Fit_GAS_GM$FilteredValues$vSigma * Fit_GAS_DJI30$FilteredValues$vSigma * mR[1, 2]
aSigma_GAS[2,1, ] = aSigma_GAS[1,2, ] 

# let's have a look at the covariance between GM and DJI
plot.ts(aSigma_GAS[1,2, ])

# point ii)
aSigma_GARCH = array(0, dim = c(2, 2, iT))

## Compute the correlation matrix

mR = cor(mRet)

# the variance of the two series
aSigma_GARCH[1,1, ] = sigma(Fit_GARCH_GM)^2
aSigma_GARCH[2,2, ] = sigma(Fit_GARCH_DJI30)^2

aSigma_GARCH[1,2, ] = sigma(Fit_GARCH_GM) * sigma(Fit_GARCH_DJI30) * mR[1, 2]
aSigma_GARCH[2,1, ] = aSigma_GARCH[1,2, ] 




## Compute the MVP weigths

mW_GAS   = matrix(0, iT, 2)
mW_GARCH = matrix(0, iT, 2)

#vector of ones
vIota = c(1, 1)

for(t in 1:iT) {
  
  mW_GAS[t, ] = solve(aSigma_GAS[,, t]) %*% vIota 
  mW_GAS[t, ] = mW_GAS[t, ]/sum(mW_GAS[t, ])
  
  mW_GARCH[t, ] = solve(aSigma_GARCH[,, t]) %*% vIota 
  mW_GARCH[t, ] = mW_GARCH[t, ]/sum(mW_GARCH[t, ])
  
}

# Point iv)

plot.ts(mW_GAS[, 1])
lines(mW_GARCH[, 1], col = "red")




