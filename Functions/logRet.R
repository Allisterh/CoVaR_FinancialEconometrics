# install.packages('rugarch')
library(rugarch)
library(quantmod)
#Load SP500 data and calculate log-returns
getSymbols("^GSPC", from="2005-01-01", to="2022-09-30")
vRet <- dailyReturn(Cl(GSPC), type='log')
plot(vRet)
