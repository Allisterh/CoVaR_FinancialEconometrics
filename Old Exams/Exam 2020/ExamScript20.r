##### Clearing space #####
rm(list = ls())

########################################
########## Computational Part ##########
########################################

####Exercise 1: Write a function to estimate the GAS–GED model derived earlier####

#A function for the score variable used in the GAS_GED_LLK function
ForcingVariableS <- function(dY, dPhi, dNu) {
  
  dS = (dNu * abs(dY) ^ dNu)/ (2 * dPhi ^ dNu) - 1.0
  return(dS)
  
}

#Likelihood of GAS GED model
#Used in the estimation function
GAS_GED_LLK <- function(vY, dOmega, dAlpha, dBeta, dNu) {
  
  #initializing the vempty vectors to store data
  iT = length(vY)
  vPhi = numeric(iT)
  vPhi_tilde = numeric(iT)
  vSigma = numeric(iT)
  vS = numeric(iT)
  vLLK = numeric(iT)
  
  
  #initialize at its unconditional value
  vPhi_tilde[1] = dOmega/(1.0 - dBeta)
  vPhi[1] = exp(vPhi_tilde[1])
  vS[1] = ForcingVariableS(vY[1], vPhi[1], dNu)
  
  #compute the log likelihood at time 1 /initial value 
  vLLK[1] = - log(2^(1+1/dNu)) - log(vPhi[1]) - lgamma(1+1/dNu) - 0.5*(abs(vY[1])^dNu)/(abs(vPhi[1])^dNu)
  
  #compute sigma at time 1
  vSigma[1] = 2^(1/dNu)*sqrt(gamma(3/dNu)/gamma(1/dNu)) * vPhi[1]
  
  for (t in 2:iT) {
    
    #update vPhi_tilde
    vPhi_tilde[t] = dOmega + dAlpha * vS[t - 1] + dBeta * vPhi_tilde[t - 1]
    
    #map
    vPhi[t] = exp(vPhi_tilde[t])
    
    #compute log likelihood at time t
    vLLK[t] = - log(2^(1+1/dNu)) - log(vPhi[t]) - lgamma(1+1/dNu) - 0.5*(abs(vY[t])^dNu)/(abs(vPhi[t])^dNu)
    
    #calculating score
    vS[t] = ForcingVariableS(vY[t], vPhi[t], dNu)
    
    #updating vSigma at time t
    vSigma[t] = 2^(1/dNu)*sqrt(gamma(3/dNu)/gamma(1/dNu)) %*% vPhi[t]
    
  }
  
  #Summing all the log likelihood contributions
  dLLK = sum(vLLK)
  
  #output: a list of the results
  lOut = list(vPhi = vPhi,
              dLLK = dLLK,
              vSigma = vSigma,
              vS = vS)
  
  return(lOut)
  
}

#Estimation of GAS GED
Estimate_GAS_GED <- function(vY) {
  
  #Setting starting values
  #within the parameter constraints
  vPar = c(omega = log(mean(vY)) * 0.1, #we isolate omega in the equation for unconditional phi
           alpha = 0.05,
           beta = 0.9, #i set beta quite high due to high persistence in the data
           nu = 2) #setting nu to start at the normal distribution
  
  #optimizer which maximizes the negatively LLK
  optimizer = optim(vPar, function(vPar, vY) {
    
    dnLLK = -GAS_GED_LLK(vY, vPar[1], vPar[2], vPar[3], vPar[4])$dLLK
    
    if (!is.finite(dnLLK)) {
      dnLLK = 1e10
    }
    
    return(dnLLK)
  
    #setting the boundaries from the derived parameter constraints
  }, lower = c(-3, 0.000001, -0.9999, 0.00001), upper = c(10, 2, 0.9999, 5),
  method = "L-BFGS-B", vY = vY)
  
  #extracting estimated parameters
  vPar = optimizer$par 
  #extracting the LLK, which is negative, as we wanted to maximize before 
  dLLK = -optimizer$value
  #Extracting more data from the output
  Filter = GAS_GED_LLK(vY = vY, dOmega = vPar[1],  dAlpha = vPar[2], dBeta = vPar[3], dNu = vPar[4])
  vPhi = Filter$vPhi
  vSigma = Filter$vSigma
  
  iK = length(vPar)
  iT = length(vY)
  #compute average BIC
  dABIC = (iK * log(iT) - 2*dLLK)/iT
  
  #output 
  lOut = list(vSigma = vSigma,
              vPhi = vPhi,
              vPar = vPar,
              dLLK = dLLK,
              dABIC = dABIC)
  
  return(lOut)
}

#####

####Exercise 2: Write a function to compute the Value–at–Risk (VaR) at level alpha = (0, 1)####

#Function to compute VaR
#Inputs:
#vY_t = vector of returns
#vSigma_t = vector of sigma
#dNu = double with nu
#dAlpha = double with the alpha level
#Output:
#VaR_t= vector of VaR
ComputeVAR <- function(vY, vSigma, dNu, dAlpha) {
  
  library(fGarch) #loading package
  
  #initializing vector
  iT = length(vSigma)
  VaR = numeric(iT)
  
  #looping through the different observations
  for (i in 1:iT) {
    
    #extracting the i'th element of the two vectors
    iY = vY[i]
    iSigma = vSigma[i]
    
    #using the ged distribution from fGarch package
    fx <- function(iY) dged(iY, nu = dNu, sd=iSigma, log = FALSE)
    
    #integrating over the cummulative density. 
    #from -inf up to the extracted iY
    cfx <- function(iY, dAlpha, FUN = fx) {
      
      integrate(function(iY, FUN = fx) FUN(iY), -Inf, iY, stop.on.error = FALSE)$value - (1 - dAlpha)
      
    }
    
    #using uniroot to find the root of the equation and then inserting into a vector
    VaR[i] <- -uniroot(cfx, c(0, 1), extendInt = "yes", tol = 0.001, dAlpha = dAlpha)$root
  }
  
  return(VaR)
  
  
}

#Alternative function to compute VAR using the rugarch package 
AltComputeVAR <- function(vSigma, dNu, dAlpha) {
  library(rugarch)
  
  dVAR = qdist(distribution = "ged", p = dAlpha, sigma = vSigma, 
               shape = dNu) 
  
  return(dVAR)
}

#####


########################################
########### Estimation Part ############
########################################

##########Part A##########
#Different model specifications

#####Exercise 1: Estimate these models (GARCH, EGARCH, gjr-GARCH) 
#####under GED assumptions for innovations (for both series) #####

#function to calculate LLK for the gas model
FunLLKGAS <- function(dY, dNu, dSigma) {
  
  #isolating phi in formula for sigma
  dPhi = dSigma / (2^(1/dNu)*sqrt(gamma(3/dNu)/gamma(1/dNu)))
  
  #using LLK function form GED
  dLLK = - log(2^(1+1/dNu)) - log(dPhi) - lgamma(1+1/dNu) - 0.5*(abs(dY)^dNu)/(abs(dPhi)^dNu)
  
  return(dLLK)
}

#GARCH(1,1) filter
GARCHFilter <- function(vY, dOmega, dAlpha, dBeta, dNu) {
  library(fGarch) #used to find the distribution of ged
  
  #initializing vector
  iT = length(vY)
  vSigma2 = numeric(iT)
  vLLK = numeric(iT)
  
  #calculating the unconditional variance 
  vSigma2[1] = dOmega/(1.0 - dAlpha - dBeta)
  #calculating likelihood contribution at time 1 with function created earlier
  vLLK[1] = FunLLKGAS(vY[1], dNu, sqrt(vSigma2[1]))
  
  #alternatively:
  #I used the dged function, alternatively one could make a function using the LLK from the theoretial part
  #vLLK[1] = dged(vY[1], 0, sqrt(vSigma2[1]), dNu, log = TRUE)
  
  for(t in 2:iT) {
    
    #calculating sigma2 at time t
    vSigma2[t] = dOmega + dAlpha * vY[t-1]^2 + dBeta * vSigma2[t-1]
    #calculating likelihood contribution at time t
    vLLK[t] = FunLLKGAS(vY[t], dNu, sqrt(vSigma2[t]))
    #vLLK[t] = dged(vY[t], 0, sqrt(vSigma2[t]), dNu,  log = TRUE)
    
  }
  
  #output: list of the results
  lOut = list(vSigma2 = vSigma2,
              vY = vY,
              vLLK= vLLK,
              dLLK = sum(vLLK))
  
  return(lOut)
  
}

#GARCH(1,1) estimation
Estimate_GARCH <- function(vY) {
  
  #load Rsolnp (used for solnp which is an optimizer that can use a ineq function)
  require(Rsolnp)
  
  #setting starting parameter values
  vPar = c("omega" = var(vY)*0.05, #isolating omega in expression for uncondtional variance
           "alpha" = 0.05, 
           "beta" = 0.9, #high persistence in the underlying data
           "nu" = 2) # #choosing nu starting at the normal distribution
  
  optimizer = solnp(vPar, function(vPar, vY) {
    
    dNLLK = -GARCHFilter(vY,
                         dOmega = vPar["omega"],
                         dAlpha = vPar["alpha"],
                         dBeta = vPar["beta"],
                         dNu = vPar["nu"])$dLLK
    
    return(dNLLK)
    #setting the boundaries from the derived parameter constraints
  }, vY = vY, LB = c(0.0001, 0.0001, 0.0001, 0.00001), UB = c(5, 0.99, 0.99, 5),
  ineqfun = function(vPar, ...) {
    vPar["alpha"] + vPar["beta"]
  }, ineqLB = 0.0001, ineqUB = 0.999)
  
  
  #extracting parameters
  vPar = optimizer$pars
  dLLK = -tail(optimizer$values, 1)
  
  Filter = GARCHFilter(vY,
                       dOmega = vPar["omega"],
                       dAlpha = vPar["alpha"],
                       dBeta = vPar["beta"],
                       dNu = vPar["nu"])
  
  #Getting volatility
  vSigma = sqrt(Filter$vSigma2)
  
  iK = length(vPar)
  iT = length(vY)
  
  #compute average BIC
  dBIC = (iK * log(iT) - 2*dLLK)
  
  #output: a list of the results
  lOut = list()
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["Filter"]] = Filter
  lOut[["dABIC"]] = dBIC/iT
  lOut[["vSigma"]] = vSigma
  
  return(lOut)

}

#EGARCH(1,1) filter
EGARCHFilter <- function(vY, dOmega, dBeta, dPhi, dPsi, dNu) {
  
  #initializing
  iT = length(vY)
  vSigma2 = numeric(iT)
  vLLK = numeric(iT)
  
  #calculating initial value of the variance
  vSigma2[1] = exp(dOmega/(1-dBeta))
  #calculating LLK contribution at time 1
  vLLK[1] = FunLLKGAS(vY[1], dNu, sqrt(vSigma2[1]))
  
  #Alternatively
  #vLLK[1] = dged(vY[1], 0, sqrt(vSigma2[1]), dNu, log = TRUE)
  
  #We know that the density of epsilon t is given as
  dLambda = sqrt( (2^(-2/dNu)) * ( gamma(1/dNu) / gamma(3/dNu) ) )
  dEpsD = 2^(1/dNu) * ( gamma(2/dNu) / gamma(1/dNu) ) * dLambda
  
  
  for(t in 2:iT) {
    #updating the variance
    vSigma2[t] = exp(dOmega + dBeta * log(vSigma2[t-1]) + dPhi * vY[t-1]/sqrt(vSigma2[t-1]) + dPsi * (abs(vY[t-1]/sqrt(vSigma2[t-1])) -  dEpsD) )
    #calculating LLK contribution at time t
    vLLK[t] = FunLLKGAS(vY[t], dNu, sqrt(vSigma2[t]))
    #vLLK[t] = dged(vY[t], 0, sqrt(vSigma2[t]), dNu,  log = TRUE)
    
  }
  
  #output: list of the results
  lOut = list()
  lOut[["vY"]] = vY
  lOut[["vSigma2"]] = vSigma2
  lOut[["vLLK"]] = vLLK
  lOut[["dLLK"]] = sum(vLLK)
  
  return(lOut)
  
}

#EGARCH(1,1) estimation
Estimate_EGARCH <- function(vY) {
  
  require(Rsolnp)
  
  #setting starting parameters
  vPar = c("omega" = log(var(vY))*0.1, #isolating omega in expression for unconditional variance
           "beta" = 0.9, #high persistence in the data
           "phi" = 0.05, #low value
           "psi" = 0.05, #low value
           "nu" = 2) #starting at normal dist
  
  optimizer = solnp(vPar, function(vPar, vY) {
    
    dNLLK = -EGARCHFilter(vY,
                         dOmega = vPar["omega"],
                         dBeta = vPar["beta"],
                         dPhi = vPar["phi"],
                         dPsi = vPar["psi"],
                         dNu = vPar["nu"])$dLLK
    
    return(dNLLK)
    #setting upper and lower bound for the a parameters
  }, vY = vY, LB = c(-1, 0.0001, -1, -1, 0.00001), UB = c(5, 0.99, 0.99, 0.99, 5))
  
  #extracting parameters
  vPar = optimizer$pars
  dLLK = -tail(optimizer$values, 1)
  
  Filter = EGARCHFilter(vY,
                       dOmega = vPar["omega"],
                       dBeta = vPar["beta"],
                       dPhi = vPar["phi"],
                       dPsi = vPar["psi"],
                       dNu = vPar["nu"])
  

  
  iK = length(vPar)
  iT = length(vY)
  
  #compute average BIC
  dBIC = (iK * log(iT) - 2*dLLK)
  
  #output a list of the results
  lOut = list()
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["Filter"]] = Filter
  lOut[["dABIC"]] = dBIC/iT
  
  
  return(lOut)
}

#gjr-GARCH(1,1) filter
GJRGARCHFilter <- function(vY, dOmega, dAlpha, dBeta, dGam, dNu) {
  
  #initializing vectors
  iT = length(vY)
  vSigma2 = numeric(iT)
  vLLK = numeric(iT)
  vS = numeric(iT)
  
  #setting starting value as the unconditional variance
  vSigma2[1] = dOmega/(1.0 - dAlpha - dBeta - 0.5 * dGam)
  #calculating LLK contribution at time 1
  vLLK[1] = FunLLKGAS(vY[1], dNu, sqrt(vSigma2[1]))
  #Alternatively
  #vLLK[1] = dged(vY[1], 0, sqrt(vSigma2[1]), dNu, log = TRUE)

  
  for(t in 2:iT) {
    #Making function for s (a dummy variable)
    dS = ifelse(vY[t-1]<0,1,0)
    #updating the variance
    vSigma2[t] = dOmega + dAlpha * vY[t-1]^2 + dGam * dS * vY[t-1]^2 + dBeta * vSigma2[t-1]
    #calculating LLK contribution at time t
    vLLK[t] = FunLLKGAS(vY[t], dNu, sqrt(vSigma2[t]))
    #vLLK[t] = dged(vY[t], 0, sqrt(vSigma2[t]), dNu,  log = TRUE)
    
  }
  
  #output: list of results
  lOut = list()
  lOut[["vY"]] = vY
  lOut[["vSigma2"]] = vSigma2
  lOut[["vLLK"]] = vLLK
  lOut[["dLLK"]] = sum(vLLK)
  
  return(lOut)
  
}

#gjr-GARCH(1,1) estimation
Estimate_GJRGARCH <- function(vY) {
  
  require(Rsolnp)
  
  #setting starting values for parameters
  vPar = c("omega" = var(vY)*0.025, #isolating omega in unconditional variance
           "alpha" = 0.05, #low value
           "beta" = 0.9, #high persistence in the data
           "gamma" = 0.05, #low value
           "nu" = 2) #starting at the normal distribution
  
  optimizer = solnp(vPar, function(vPar, vY) {
    
    dNLLK = -GJRGARCHFilter(vY,
                         dOmega = vPar["omega"],
                         dAlpha = vPar["alpha"],
                         dBeta = vPar["beta"],
                         dGam = vPar["gamma"],
                         dNu = vPar["nu"])$dLLK
    
    return(dNLLK)
    
    #setting boundaries under the derived constraints
  }, vY = vY, LB = c(0.0001, 0.0001, 0.0001, 0.0001, 0.00001), UB = c(5, 0.99, 0.99, 0.99, 5),
  ineqfun = function(vPar, ...) {
    vPar["alpha"] + vPar["beta"] + 0.5 * vPar["gamma"]
  }, ineqLB = 0.0001, ineqUB = 0.999)
  
  #extracting data
  vPar = optimizer$pars
  dLLK = -tail(optimizer$values, 1)
  
  Filter = GJRGARCHFilter(vY,
                       dOmega = vPar["omega"],
                       dAlpha = vPar["alpha"],
                       dBeta = vPar["beta"],
                       dGam = vPar["gamma"],
                       dNu = vPar["nu"])
  
  iK = length(vPar)
  iT = length(vY)
  
  #compute average BIC
  dBIC = (iK * log(iT) - 2*dLLK)
  
  
  #output: list of results
  lOut = list()
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["Filter"]] = Filter
  lOut[["dABIC"]] = dBIC/iT

  return(lOut)
  
}
#####


#####Exercise 2: Compare the filtered volatilities according to these three specifications
#####(GARCH, EGARCH, gjr-GARCH) in a figure. #####

#Downloading the data
library("quantmod") #Using the quantmod package
GSPC = getSymbols("^GSPC", from = "2007-01-03", to = "2019-01-01", auto.assign = FALSE)
DJI = getSymbols("^DJI", from = "2007-01-03", to = "2019-01-01", auto.assign = FALSE)

#mergning the two together
mP = merge(DJI$DJI.Adjusted,
           GSPC$GSPC.Adjusted)

#calculating the log differences 
vY_DJI = diff(log(as.numeric(mP[,1]))) * 100
vY_GSPC = diff(log(as.numeric(mP[,2]))) * 100

#column binding the two vectors together to a matrix
mY = cbind(vY_DJI, vY_GSPC)

#estimating the three models on DJI data
garchFIT_DJI <- Estimate_GARCH(mY[,1])
egarchFIT_DJI <- Estimate_EGARCH(mY[,1]) 
gjrgarchFIT_DJI <- Estimate_GJRGARCH(mY[,1])

#estimating the three models on S&P500 data
garchFIT_SPX <- Estimate_GARCH(mY[,2])
egarchFIT_SPX <- Estimate_EGARCH(mY[,2]) 
gjrgarchFIT_SPX <- Estimate_GJRGARCH(mY[,2])

#comparing the filtered volatilities in a figure
par(mfrow = c(2, 1))
par(mar=c(2,4,1,1))
plot.ts(garchFIT_SPX$vSigma, col = 1, ylab="filtered volatility", main = "S&P500")
lines(sqrt(egarchFIT_SPX$Filter$vSigma2), col = 2)
lines(sqrt(gjrgarchFIT_SPX$Filter$vSigma2), col = 4)
legend(2300, 6,legend = c("GARCH(1,1)", "EGARCH(1,1)", "gjr-GARCH(1,1)"), col = c(1,2,4), lty = 1, cex = 0.5, lwd = 3)
plot.ts(garchFIT_DJI$vSigma, col = 1, ylab="filtered volatility", main = "Dow Jones Industrial")
lines(sqrt(egarchFIT_DJI$Filter$vSigma2), col = 2)
lines(sqrt(gjrgarchFIT_DJI$Filter$vSigma2), col = 4)
legend(2300, 5.5, legend = c("GARCH(1,1)", "EGARCH(1,1)", "gjr-GARCH(1,1)"), col = c(1,2,4), lty = 1, cex = 0.5, lwd = 3)

#we notice that the filtered volatility under these three model specifications are quite similar

#####

#####Exercise 3: Select the best model using the BIC criteria.#####

#We make a list of the Average BIC on SPX and SJI data, respetively
#notice: I use the average, but it does not change anything, as it just divided by the same number
lBIC_SPX = list()
lBIC_SPX[["BIC_GARCH"]] = garchFIT_SPX$dABIC
lBIC_SPX[["BIC_EGARCH"]] = egarchFIT_SPX$dABIC
lBIC_SPX[["BIC_GJRGARCH"]] = gjrgarchFIT_SPX$dABIC
lBIC_SPX

lBIC_DJI = list()
lBIC_DJI[["BIC_GARCH"]] = garchFIT_DJI$dABIC
lBIC_DJI[["BIC_EGARCH"]] = egarchFIT_DJI$dABIC
lBIC_DJI[["BIC_GJRGARCH"]] = gjrgarchFIT_DJI$dABIC
lBIC_DJI

#From the BIC we would choose the EGARCH specification as it is the lowest (for both data sets)

#####

##########Part B##########
#Estimate the GAS-GED model on the DJI and SP500 series.

#####Exercise 1: Compare the filtered volatility of the GAS-GED model 
#####with that of the GARCH(1,1) model in a figure. #####

#Estimating using S&P 500 and DJI data
par(mfrow = c(2, 1))

garchFIT_DJI <- Estimate_GARCH(mY[,1])
gas_gedFIT_DJI = Estimate_GAS_GED(mY[,1])
garchFIT_SPX <- Estimate_GARCH(mY[,2])
gas_gedFIT_SPX = Estimate_GAS_GED(mY[,2])


#Plotting results
plot.ts(sqrt(garchFIT_SPX$Filter$vSigma2), ylab = "Filtered Volatility", main = "S&P500")
lines(gas_gedFIT_SPX$vSigma, col = "red")
legend(2200, 6, legend = c("GARCH(1,1)","GAS-GED"), col = 1:2, lty = 1, cex = 0.6, lwd = 3)

plot.ts(sqrt(garchFIT_DJI$Filter$vSigma2), ylab = "Filtered Volatility", main = "Dow Jones Industrial")
lines(gas_gedFIT_DJI$vSigma, col = "red")
legend(2200, 5.5, legend = c("GARCH(1,1)","GAS-GED"), col = 1:2, lty = 1, cex = 0.6, lwd = 3)

#####

#####Exercise 2: Compare VaR at levels alpha = 1% and alpha = 5% estimated with the 
#####GAS-GED and GARCH(1,1) models in a figure. #####

#Now we need to compute VAR with the function we made earlier
#par(mfrow = c(1, 1))
#testspx005 <- AltComputeVAR(mY[,1], vSigma = garchFIT_SPX$vSigma, dNu = garchFIT_SPX$vPar[4], dAlpha = 0.05)
#plot.ts(testspx005, main = "alpha=0.05", col = "1")
#lines(vSPX005_GARCH, col = "2")
#legend("bottomright", legend = c("GARCH(1,1)", "GAS-GED"), col = 1:2, lty = 1, cex = 0.5, lwd = 3 )


##FIRST for GSPC
#Vectors of VaR for the GAS-GED model for DJI with alpha =0.05 and 0.01, respectively
vDJI005_GASGED <- ComputeVAR(mY[,1], vSigma = gas_gedFIT_SPX$vSigma, dNu = gas_gedFIT_SPX$vPar[4], dAlpha = 0.05)
vDJI001_GASGED <- ComputeVAR(mY[,1], vSigma = gas_gedFIT_SPX$vSigma, dNu = gas_gedFIT_SPX$vPar[4], dAlpha = 0.01)
#Vectors of VaR for the GARCH model for DJI with alpha =0.05 and 0.01, respectively
vDJI005_GARCH <- ComputeVAR(mY[,1], vSigma = sqrt(garchFIT_SPX$Filter$vSigma2), dNu = garchFIT_SPX$vPar[4], dAlpha = 0.05)
vDJI001_GARCH <- ComputeVAR(mY[,1], vSigma = sqrt(garchFIT_SPX$Filter$vSigma2), dNu = garchFIT_SPX$vPar[4], dAlpha = 0.01)

#Vectors of VaR for the GAS-GED model for GSPC with alpha =0.05 and 0.01, respectively
vSPX005_GASGED <- ComputeVAR(mY[,2], vSigma = gas_gedFIT_SPX$vSigma, dNu = gas_gedFIT_SPX$vPar[4], dAlpha = 0.05)
vSPX001_GASGED <- ComputeVAR(mY[,2], vSigma = gas_gedFIT_SPX$vSigma, dNu = gas_gedFIT_SPX$vPar[4], dAlpha = 0.01)
#Vectors of VaR for the GARCH model for GSPC with alpha =0.05 and 0.01, respectively
vSPX005_GARCH <- ComputeVAR(mY[,2], vSigma = sqrt(garchFIT_SPX$Filter$vSigma2), dNu = garchFIT_SPX$vPar[4], dAlpha = 0.05)
vSPX001_GARCH <- ComputeVAR(mY[,2], vSigma = sqrt(garchFIT_SPX$Filter$vSigma2), dNu = garchFIT_SPX$vPar[4], dAlpha = 0.01)



#graph settings
par(mfrow = c(2, 1))
par(mar=c(4,4,1,1))

#Plotting VAR for SPX (alpha = 0.01, 0.05)
plot.ts(vSPX001_GARCH, ylab = "VaR SPX", main = "alpha=0.01", col = "1")
lines(vSPX001_GASGED, col = "2")
legend(2300, -12.25, legend = c("GARCH(1,1)", "GAS-GED"), col = 1:2, lty = 1, cex = 0.5, lwd = 3)
plot.ts(vSPX005_GARCH, ylab = "VaR SPX", main = "alpha=0.05", col = "1")
lines(vSPX005_GASGED, col = "2")
legend(2300, -7.5, legend = c("GARCH(1,1)", "GAS-GED"), col = 1:2, lty = 1, cex = 0.5, lwd = 3)

#Plotting VAR for DJI (alpha = 0.01, 0.05)
plot.ts(vDJI001_GARCH, ylab = "VaR DJI", main = "alpha=0.01", col = "1")
lines(vDJI001_GASGED, col = "2")
legend(2300, -12.25, legend = c("GARCH(1,1)", "GAS-GED"), col = 1:2, lty = 1, cex = 0.5, lwd = 3)
plot.ts(vDJI005_GARCH, ylab = "VaR DJI", main = "alpha=0.05", col = "1")
lines(vDJI005_GASGED, col = "2")
legend(2300, -7.5, legend = c("GARCH(1,1)", "GAS-GED"), col = 1:2, lty = 1, cex = 0.5, lwd = 3)

#####

#####Exercise 3: Select the best model between GAS-GED and GARCH(1,1) using BIC. #####

#Extracting 
lFit = list()
lFit[["SPX BIC GARCH(1,1)"]] = garchFIT_SPX$dABIC
lFit[["SPX BIC GAS GED"]] = gas_gedFIT_SPX$dABIC
lFit[["DJI BIC GARCH(1,1)"]] = garchFIT_DJI$dABIC
lFit[["DJI BIC GAS GED"]] = gas_gedFIT_DJI$dABIC
lFit

#The GARCH model specification gives us the lowest BIC for both dataset,
# and therefore, we would select that model.


#####


##########Part C##########
#Calculating covariance matrices and portfolio weights

#####Exercise 1: Compute the covariance matrix for each t when DJI and GSPC 
#####both follow the GAS-GED model##### 

#I made a function to calculate both the sigma and the weights, and will therefore use it multiple 
#times throughout Part C
CalculateWeights <- function(mY, model=Estimate_GARCH) {
  
  ## estimate the marginal GARCH models
  require(Rsolnp)
  
  #Initializing
  iT = nrow(mY)
  aS <- array(0, dim=c(2, 2, iT)) #array of Sigma
  mW <- matrix(NA, 2, iT, dimnames = list(c("DJI", "GSPC"), NULL)) #matrix of weights
  lFit_univariate = list()
  
  #estimate the univariate GARCH/GAS-GED models
  for(n in 1:ncol(mY)) {
    lFit_univariate[[n]] = model(mY[, n])
  }
  
  #Extract volatility from the given model specification
  mSigma = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    Fit$vSigma
  }))
  
  #Calculate correlation (R_t)
  mR = cor(mY)
  
  #Function for computing the minimum variance portfolio weights
  Compute_MVP <- function(mSigma) {
    #Computes the Minimum variance portfolio weights from the given Sigma matrix
    #Using formulas from lecture 15, slide 3
    
    iN = ncol(mSigma)
    vOnes = rep(1, iN)
    
    vOmega = vOnes %*% solve(mSigma)
    vOmega = vOmega/sum(vOmega)
    
    return(vOmega)
  }
  
  #Running a loop where we calculate the Sigma and weight for all t
  for (t in 1:iT) {
    mD = diag(mSigma[t, ])
    aS[,,t] = mD %*% mR %*% mD
    mW[,t] = Compute_MVP(aS[,,t])
  }
  
  #returning output which is 
  #an array (2x2xiT) of the Sigma 
  #a matrix (2xiT) of the portfolio weights
  lOut = list(aS = aS,
              mW = mW)
  
  return(lOut)
}

#extracting the Sigma from the function
#we choose the model to be estimated to be GAS-GED
#notice output is an array so WGAS_GED[,,t] is the covariance matrix a time t
aCov_GAS_GED = CalculateWeights(mY, model=Estimate_GAS_GED)$aS

#####

#####Exercise 2: Compute the covariance matrix for each t when DJI and GSPC 
#####both follow the GARCH(1,1) model#####

#Extracting the Sigma from the function
aCov_GARCH = CalculateWeights(mY, model=Estimate_GARCH)$aS

#####

#####Exercise 3: Compute portfolio weights for each point in time #####

#Extracting the weights from the function just created
mWGARCH = CalculateWeights(mY, model=Estimate_GARCH)$mW
mWGAS_GED = CalculateWeights(mY, model=Estimate_GAS_GED)$mW

#####

#####Exercise 4: Compare the portfolio weights of the two models in a figure #####

#Making plots of the weights for S&P500 and DJI, using the two model specifications GARCH and GAS-GED
par(mfrow = c(2, 1))
par(mar=c(2,4,1,1))
plot.ts(mWGARCH[2,], col = "black", ylab="w_GSPC", main = "S&P500")
lines(mWGAS_GED[2,], col = "red")
legend(2300, 3.4, legend = c("GARCH(1,1)", "GAS-GED"), col = 1:2, lty = 1, cex = 0.5, lwd = 3)
plot.ts(mWGARCH[1,], col = "black", ylab="w_DJI", main = "Dow Jones Industrial")
lines(mWGAS_GED[1,], col = "red")
legend(2300, 3.4, legend = c("GARCH(1,1)", "GAS-GED"), col = 1:2, lty = 1, cex = 0.5, lwd = 3)
#####
 