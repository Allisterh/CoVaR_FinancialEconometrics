rm(list=ls())
if(!is.null(dev.list())) dev.off()

setwd("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/")

library(quantmod)
library(ggplot2)
library(ggthemes)

## download the series
GSPC = getSymbols("^GSPC", from = "2018-01-01", to = "2019-01-01", auto.assign = FALSE)
DJI = getSymbols("^DJI", from = "2018-01-01", to = "2019-01-01", auto.assign = FALSE)

mP = merge(GSPC$GSPC.Adjusted,
           DJI$DJI.Adjusted)
  
vY_GSPC = diff(log(as.numeric(mP[,1]))) * 100
vY_DJI = diff(log(as.numeric(mP[,2]))) * 100

mY = cbind(vY_GSPC, vY_DJI)
# Import functions for DCC and CCC estimation
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/CCC.r")
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/DCC.r")

Fit_DCC = Estimate_DCC(mY)
Fit_CCC = Estimate_CCC(mY)



iT = nrow(mY)
iN = ncol(mY)
#array with portfolio weights

aW = array(NA, dim = c(iT, iN, 2), dimnames = list(NULL, c("GPSC", "DJI"), c("DCC", "CCC")))
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/MVP.r")
for (t in 1:iT) {

  mR_DCC = Fit_DCC$aCor[,,t]
  mD_DCC = diag(Fit_DCC$mSigma[t, ])
  mSigma_DCC = mD_DCC %*% mR_DCC %*% mD_DCC

  aW[t,,"DCC"] = Compute_MVP(mSigma_DCC)

  mR_CCC = Fit_CCC$mR
  mD_CCC = diag(Fit_CCC$mSigma[t, ])
  mSigma_CCC = mD_CCC %*% mR_CCC %*% mD_CCC

  aW[t,,"CCC"] = Compute_MVP(mSigma_CCC)
}

# Plot the weights over time ggplot

ggplot(data = data.frame(aW[,1,"DCC"]), aes(x = 1:iT, y = aW[,1,"DCC"], col = "DCC")) + 
  geom_line() +
  geom_line(data = data.frame(aW[,1,"CCC"]), aes(x = 1:iT, y = aW[,1,"CCC"], col = "CCC")) +
  labs(x = "Time", y = "Weight", title = "Portfolio weights for the DJI vs SP500") +
  theme_economist() +
  theme(legend.title = element_blank())
ggsave("./Old Exams/Weights.pdf")



### Compute the CoVaR ###
# This does not work, no unit root in the data.
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/CoVaR.r")

iT = nrow(mY)

aSigma_DCC = array(NA, dim = c(iN, iN, iT))
vMu = rep(0, iN)

mCoVaR = matrix(NA, iT, 2, dimnames = list(NULL, c("alpha = 1%", "alpha = 5%")))

for (t in 1:iT) {

  mR_DCC = Fit_DCC$aCor[,,t]
  mD_DCC = diag(Fit_DCC$mSigma[t, ])
  aSigma_DCC[,,t] = mD_DCC %*% mR_DCC %*% mD_DCC

  mCoVaR[t, "alpha = 1%"] = ComputeCoVaR(dAlpha = 0.01, vMu, aSigma_DCC[,,t])
  mCoVaR[t, "alpha = 5%"] = ComputeCoVaR(dAlpha = 0.05, vMu, aSigma_DCC[,,t])

}



