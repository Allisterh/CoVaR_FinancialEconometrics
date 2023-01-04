rm(list=ls())

# load quantmod
library(quantmod)
library(ggplot2)
library(ggthemes)

#download the series
VIX = getSymbols("^VIX", from = "2010-01-01", to = "2019-01-01", auto.assign = FALSE)

head(VIX)

#extract the levels
vY = as.numeric(VIX$VIX.Adjusted)

source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/GammaGAS.r")
Fit = Estimate_GammaGAS(vY)

## estimated parameters
Fit$vPar
Fit$Filter$vMu


## Plot vY and the filtered mean vMu in ggplot
ggplot(data = data.frame(vY, vMu = Fit$Filter$vMu), 
  aes(x = 1:length(vY), y = vY, color = "vY")) + 
  geom_line() + 
  geom_line(aes(x = 1:length(vY), y = Fit$Filter$vMu, color = "vMu")) +
  labs(title = "The VIX Series and the Filtered Mean", 
       x = "Time", y = "VIX") +
  theme_economist() +
  theme(legend.title=element_blank())       
ggsave("./Old Exams/VIX and Filtered Mean.pdf")




## Estimate the static model
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/StaticGamma.r")
Fit_Static = Estimate_StaticGamma(vY)

Fit$vPar
Fit_Static$vPar

Fit$BIC
Fit_Static$BIC

## Estimante MEM
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/GammaMEM.r")
Fit_MEM = Estimate_MEM(vY)

Fit_MEM$vPar



lines(Fit_MEM$Filter$vMu, col = "red")

Fit_MEM$dLLK
Fit$dLLK

vVar = Fit_MEM$Filter$vMu^2 /Fit_MEM$vPar["a"]

# Plot the VIX, filtered mean and variance in ggplot
ggplot(data = data.frame(vY, vMu = Fit_MEM$Filter$vMu, vVar = vVar), 
  aes(x = 1:length(vY), y = vY, color = "vY")) + 
  geom_line() + 
  geom_line(aes(x = 1:length(vY), y = Fit_MEM$Filter$vMu, color = "vMu")) +
  geom_line(aes(x = 1:length(vY), y = vVar, color = "vVar")) +
  labs(title = "The VIX Series and the Filtered Mean and Variance", 
       x = "Time", y = "VIX") +
  theme_economist() + 
  theme(legend.title=element_blank())
ggsave("./Old Exams/VIX, Filtered Mean and Variance.pdf")


