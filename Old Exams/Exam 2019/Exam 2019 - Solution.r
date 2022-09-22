rm(list = ls()) # clear environment
cat("\014")     # clear console output
setwd("C:/Users/Christian Restrup/Dropbox/7. Semester/Financial Econometrics/�velser/Data+Additional Code")
##########################################################################################
############################ Financial Econometrics Exam 2019 ############################
##########################################################################################

############################################################################
###################### Question 1: Computational part ######################
############################################################################

library("quantmod")
library("Rsolnp")
library("mvtnorm")

# The density function of the Scaled Gamma Distribution. 
# Input:        Vector of observations, vY, a vector of filtered means, vMu, scale parameter,
#               dA, and bLog, which indicate whether or not a log-density is considered.
# Output:       Vector of the (log) Gamma PDF, vLogGammaPDF.
GammaPDF <- function(vY, vMu, dA, bLog = TRUE) {
    
    vLogGammaPDF = - lgamma(dA) + dA * log(dA) + (dA - 1) * log(vY) - dA * log(vMu) - dA * vY / vMu   # (Eq.7)
                    #lgamma being the log gamma function, take the GAMMA pdf and take log, write it out
    if (!bLog) {
        
        vLogGammaPDF = exp(vLogGammaPDF)
    }
    
    return(vLogGammaPDF)
}

# Filter function of GAMMA-G AS model.
# Input:        Vector of observations, vY, the paramaters vMu depends upon,
#               dOmega, dAlpha and dBeta, as well as the scale parameter, dA.
# Output:       The loglikelihood of the filter as well as the filtered means vMu.
GAMMA_GAS_filter <- function(vY, dOmega, dAlpha, dBeta, dA) {
    
    iT      = length(vY)
    vMuT    = numeric(iT)
    vMu     = numeric(iT)
    vUT     = numeric(iT)
    
    # Initialising of Mu_tilde by the unconditional expectation
    vMuT[1]   = dOmega/(1 - dBeta)   # (Eq.9)
    vMu[1]    = exp(vMuT[1])         # (Eq.1) Mapping function
    
    # Updating step, u = u_tilde
    vUT[1]    = sqrt(dA)*(vY[1]/exp(vMuT[1]) - 1) # (Eq.3)
    
    dLLK = GammaPDF(vY[1], vMu[1], dA) # The first log likelihood contribution (Eq. 7)
  
    for (t in 2:iT) {
        
        vMuT[t]   = dOmega + dAlpha * vUT[t-1] + dBeta * vMuT[t-1] # (Eq.5)
        vMu[t]    = exp(vMuT[t])                                   # (Eq.1) mapping function
        vUT[t]    = sqrt(dA)*(vY[t]/exp(vMuT[t]) - 1)              # (Eq.6)
        
        dLLK      = dLLK + GammaPDF(vY[t], vMu[t], dA)             # Summing the loglikelihoods -> (Eq.8)
    }
    
    lOut = list()
    lOut[["vMu"]] = vMu
    lOut[["dLLK"]] = dLLK
    
    return(lOut)
}

# The estimation function of the GAMMA-GAS model.
# Input:        Vector of observations, vY.
# Output:       Filtered means, optimized parameter coefficients
#               for parameters: dOmega, dAlpha, dBeta and dA,
#               Loglikelihood as well as BIC.
GAMMA_GAS_estimation <- function(vY) {
    
    # Starting values, does not really fit into a system, but most of course lie within the lower and upper bound
    vPar = c(omega = 0.2, alpha = 0.05, beta = 0.9, a = 15  )
           #   vPar[1],     vPar[2],      vPar[3],  vPar[4]
    # Optimization of the filter
    optimization = solnp(vPar, function(vPar, vY) { # Nonlinear optimization using augmented Lagrange method
                      
        dNegLLK = - GAMMA_GAS_filter(vY, dOmega = vPar[1],
                                     dAlpha = vPar[2], 
                                     dBeta = vPar[3], 
                                     dA = vPar[4])$dLLK
        
        if (!is.finite(dNegLLK)) {
            dNegLLK = 1e05
        }
      
        return(dNegLLK)
    }, LB = c(-0.5,     0.001,      0.01,      0.1), # Contraints from question
       UB = c(0.5,      1.5,        0.999,     300),
            #{Omega},   {Alpha},    {Beta},    {a}
            #vPar[1],   vPar[2],    vPar[3], vPar[4]
    
    vY = vY)
    
    vPar    = optimization$pars  # Pars within optimization is a list of the optimal parameters, which is accessed this way 
    
    Filter  = GAMMA_GAS_filter(vY,  vPar[1],     vPar[2],    vPar[3],   vPar[4])
                                    #{Omega},   {Alpha},    {Beta},    {a}
                                    # The optimal values used in our filter function, which output vMu and the dllk
    lOut = list()
    lOut[["vMu"]]     = Filter$vMu                                  # From the filter list, take out the value vMu, which comes from the GAMMA-GAS filter
    lOut[["vPar"]]    = vPar                                        # The optimal set of parameters, those that maximize the negative likelihood
    lOut[["LLK"]]     = -tail(optimization$value, 1)                # The Last value of the optimazation
    lOut[["BIC"]]     = log(length(vY)) * 4 - 2 * lOut[["LLK"]]     # The BIC
    
    return(lOut)
}

# Gamma_GAS_Constraint - Filter
# Description:  Gamma Filter for the case with constant mean, dMu
# Input:        Vector of observations, vY, mean, dMu and scale parameter dA
#               of the gamma distribution.
# Output:       Loglikelihood
GAMMA_filter <- function(vY, dA, dMu) {  # Notice the HUGE difference from the GAMMA_GAS function
                                         # There is dMu as unput to the function here, but not in the GAMMA_GAS
                                         # Because we hold it constant here
    
    dLLK = GammaPDF(vY[1], dMu, dA)
    
    for (t in 2:length(vY)) {
        
        dLLK      = dLLK + GammaPDF(vY[t], dMu, dA)
    }
    
    lOut = list()
    lOut[["dLLK"]] = dLLK
    
    return(lOut)
}

# Gamma_GAS_Constraint - Estimation
# Description:  The estimation function of the GAMMA-GAS model with dMu constant. 
# Input:        Vector of observations, vY
# Output:       Optimized parameter coefficients for parameters: dMu and dA,
#               Loglikelihood as well as BIC.
GAMMA_GAS_Cons_estimation <- function(vY) {
    
    # Starting values
    vPar = c(a = 15, mu = 10)
    
    # Optimization of the filter
    optimization = solnp(vPar, function(vPar, vY) {
        
        dNegLLK = - GAMMA_filter(vY, dA = vPar[1],
                                 dMu = vPar[2])$dLLK
        
        if (!is.finite(dNegLLK)) {
            dNegLLK = 1e05
        }
        
        return(dNegLLK)
    }, LB = c(0.1, 0.0001), UB = c(300, 300),
    vY = vY)
    
    vPar    = optimization$pars
    
    Filter  = GAMMA_filter(vY, vPar[1], vPar[2])
    
    lOut = list()
    lOut[["vPar"]]    = vPar
    lOut[["LLK"]]     = -tail(optimization$value, 1)
    lOut[["BIC"]]     = log(length(vY)) * 2 - 2 * lOut[["LLK"]]
    
    return(lOut)
    
}

# Multiplicative Error Model - Filter
# Description:  MEM by Engle and Gallo (2006).
# Input:        Vector of observations, vY,parameters for the 
#               updating step of Mu, and the scale parameter, a.
# Output:       Filtered means and Loglikelihood.
MEM_filter <- function(vY, dKappa, dEta, dPhi, dA) {
    
    iT      = length(vY)
    vMu     = numeric(iT)
    
    # Initialising of Mu_tilde by the unconditional expectation
    vMu[1]  = dKappa / (1 - dEta - dPhi)  # (Eq.11)
    dLLK    = GammaPDF(vY[1], vMu[1], dA) # Initial likelihood
    
    for (t in 2:iT) {
        
        vMu[t]    = dKappa + dEta * vY[t - 1] + dPhi * vMu[t - 1] #(Eq.0)
        dLLK      = dLLK + GammaPDF(vY[t], vMu[t], dA) # Summing the likelihood
    }
    
    lOut = list()
    lOut[["vMu"]] = vMu
    lOut[["dLLK"]] = dLLK
    
    return(lOut)
}

# Multiplicative Error Model - Estimation
# Description:  MEM by Engle and Gallo (2006).
# Input:        Vector of observations, vY.
# Output:       Filtered means, optimized parameter coefficients
#               for parameters: dKappa, dEta, dPhi and dA,
#               the conditional variances of the model, vSigma2,
#               Loglikelihood as well as BIC.
MEM_estimation <- function(vY) {
    
    # Starting values
    vPar = c(kappa = 3, eta = 0.2, phi = 0.2, a = 15)
    
    # Optimization of the filter
    optimization = solnp(vPar, function(vPar, vY) {
        
        dNegLLK = - MEM_filter(vY, dKappa = vPar[1],
                               dEta = vPar[2], 
                               dPhi = vPar[3], 
                               dA = vPar[4])$dLLK
        
        if (!is.finite(dNegLLK)) {
            dNegLLK = 1e05
        }
        
        return(dNegLLK)
    }, ineqfun = function(vPar, ...) {
        sum(c(vPar[2],vPar[3])) 
    }, ineqLB = 0.01, ineqUB = 0.999,
    LB = c(0.1, 0.01, 0.01, 0.1), UB = c(10, 0.99, 0.99, 300),
    vY = vY)
    
    vPar    = optimization$pars
    
    Filter  = MEM_filter(vY, vPar[1], vPar[2], vPar[3], vPar[4])
    
    lOut = list()
    lOut[["vMu"]]     = Filter$vMu
    lOut[["vSigma2"]] = Filter$vMu^2 / vPar[4]
    lOut[["vSigma"]]  = sqrt(Filter$vMu^2 / vPar[4])
    lOut[["vPar"]]    = vPar
    lOut[["LLK"]]     = -tail(optimization$value, 1)
    lOut[["BIC"]]     = log(length(vY)) * 4 - 2 * lOut[["LLK"]]
    
    return(lOut)
}

############################################################################
######################## Question 1: Emperical part ########################
############################################################################

# Load data from Yahoo Finance
mData_Yahoo = getSymbols(Symbols = "^VIX",
                         from = "2010-01-01",
                         to = "2019-01-01",
                         auto.assign = FALSE)

# Adding dates as rownames to the adjusted VIX volatility
vDates  = as.Date(row.names(as.matrix(mData_Yahoo)))
vVIX    = as.numeric(mData_Yahoo[,6])

# Estimation of our 3 models by above functions.
Estimate_GAS        = GAMMA_GAS_estimation(vY = vVIX)
Estimate_GAS_CONS   = GAMMA_GAS_Cons_estimation(vY = vVIX)
Estimate_MEM        = MEM_estimation(vY = vVIX)
# Out of this is as
# Iter: 1 fn: 3756.5457	 Pars:    0.41005   0.94574   0.03229 162.87418

# Iter: 2 fn: 3756.5457	 Pars:    0.41000   0.94569   0.03234 162.85934
#             {llk}               {kappa}    {eta}    {phi}     {a}  
# comes out the same order as it goes in.
# solnp--> Completed in 2 iterations

vBIC = c(
    BIC_GAS           = Estimate_GAS$BIC,
    BIC_CONS          = Estimate_GAS_CONS$BIC
)

# Sort our models by BIC - Estimate_GAS preferable.
sort(vBIC)

vLLK = c(
    LLK_GAS           = Estimate_GAS$LLK,
    LLK_CONS          = Estimate_GAS_CONS$LLK,
    LLK_MEM           = Estimate_MEM$LLK
)

# Sort our models by Loglikelihood - Estimate_MEM preferable.
sort(vLLK)

#######################################
########### MEM 3 x 1 plot ############
#######################################

# Plot MEM models as it maximixes the loglikelihood - 3 x 1 plot
par(mfrow = c(3, 1))
par(mar = c(3,3,3,2))
plot(vDates, vVIX,                   t = "l", main = "VIX")                                                 #1
plot(vDates, Estimate_MEM$vMu,       t = "l", main = expression(paste("Filtered Conditional Means, ",mu)))  #2
plot(vDates, Estimate_MEM$vSigma2,   t = "h", main = expression(paste("Conditional Variance, ", sigma^2)))  #3

# 1 and 2 are equal?? why?



############################################################################
################################ Question 2 ################################
############################################################################

# Load data as well as dates (dates loaded from Yahoo Finance to map them against data)
mData   = read.table(file = "data_new.csv", header = TRUE, sep = ";", dec = ",")

mGSPC   = getSymbols(Symbols = "^GSPC",
                     from = "2007-01-04",     # Note the first date we start from would be the 4th, as returns is based on previous day.
                     to = "2019-01-01",
                     auto.assign = FALSE)

vDates_GSPC  = as.Date(row.names(as.matrix(mGSPC)))
row.names(mData) = vDates_GSPC

# Returns plotted against time
par(mfrow = c(2,1))
par(mar = c(5,5,2,2))
plot(vDates_GSPC, mData[,1], main ="GSPC returns",  xlab = "Time", ylab = "Return", type = "l")
plot(vDates_GSPC, mData[,2], main ="DJI30 returns", xlab = "Time", ylab = "Return", type = "l")

# Returns plotted against each other
par(mfrow = c(1,1))
par(mar = c(5,5,3,3))
plot(mData[,1],mData[,2], main ="GSPC plotted against DJI30", xlab = "Return GSPC", ylab = "Return DJI30", xlim = c(-10,10), ylim = c(-10,10))



############################################################################
###################### Question 2: Computational part ######################
############################################################################

# GARCH(1,1) filter
# Input:        Vector of observations, vY, as well as parameter values
#               in the updating step for sigma^2.
# Output:       Estimated volatility of the garch(1,1) model and likelihood.



GARCH_Filter <- function(vY, dOmega, dAlpha, dBeta) {    # Exact same function, as the GarchSim in Exercise 1.
    
    iT      = length(vY)
    vSigma2 = numeric(iT)
    
    # Initialising of Sigma^2 by the unconditional expectation (assuming stationarity)
    vSigma2[1]  = dOmega / (1 - dAlpha - dBeta)  # Unconditional Variance
    
    dLLK = dnorm(vY[1], mean = 0, sd = sqrt(vSigma2[1]), log = TRUE) # Initialise the likelihood
    
    for (t in 2:iT) {
        
        vSigma2[t]  = dOmega + dAlpha * vY[t - 1]^2 + dBeta * vSigma2[t - 1] # running over the simple GARCH(1,1) model.
        dLLK = dLLK + dnorm(vY[t], mean = 0, sd = sqrt(vSigma2[t]), log = TRUE)
    }
    
    lOut = list()
    lOut[["vSigma2"]] = vSigma2
    lOut[["vSigma"]]  = sqrt(vSigma2)
    lOut[["dLLK"]]    = dLLK
    
    return(lOut)
}

# GARCH(1,1) estimation
# Input:        Vector of observations.
# Output:       Filtered volatility, optimized parameter coefficients
#               for parameters: dOmega, dAlpha and dBeta,
#               Loglikelihood as well as BIC.
GARCH_Estimation <- function(vY) {
    
    # Starting values
    vPar = c(dOmega = 3, dAlpha = 0.2, dBeta = 0.2)
    
    # Optimization of the filter.
    optimization = solnp(vPar, function(vPar, vY) {
        
        dNegLLK = - GARCH_Filter(vY, dOmega = vPar[1],
                                     dAlpha = vPar[2], 
                                     dBeta = vPar[3])$dLLK
        
        if (!is.finite(dNegLLK)) {
            dNegLLK = 1e05
        }
        
        return(dNegLLK)
    }, ineqfun = function(vPar, ...) {
        sum(c(vPar[2],vPar[3])) 
    }, ineqLB = 0.0001, ineqUB = 0.99,
    LB = c(0.0001, 0.0001, 0.0001), UB = c(100, 0.99, 0.99),
    vY = vY)
    
    vPar    = optimization$pars
    
    Filter  = GARCH_Filter(vY, vPar[1], vPar[2], vPar[3])
    
    lOut = list()
    lOut[["vSigma2"]] = Filter$vSigma2
    lOut[["vSigma"]]  = Filter$vSigma
    lOut[["vPar"]]    = vPar
    lOut[["LLK"]]     = -tail(optimization$value, 1)
    lOut[["BIC"]]     = log(length(vY)) * 3 - 2 * lOut[["LLK"]]
    
    return(lOut)
}

# Estimation of the Constant Conditional Correlation (CCC) Model
# Input:        Vector of observations.
# Output:       Optimized parameter coefficients and volatility
#               for the marginal models, correlation (constant) and
#               covariance matrix for the CCC model, standardized returns,
#               Loglikelihood as well as BIC.
Estimate_CCC <- function(mY) { # Almost the same as the CCC in exercise 6, but without the use of "rugarch"
    
    # List for storage
    lFit_univariate = list()
    
    # Univariate GARCH models placed in the same list
    for(n in 1:ncol(mY)) {
        lFit_univariate[[n]] = GARCH_Estimation(mY[, n])
    }
    
    # Standardized returns - initialising with matrix
    mEta = matrix(0,nrow = nrow(mData), ncol = ncol(mY))
    
    for(n in 1:ncol(mY)) {
        mEta[,n] <- mData[,n] / lFit_univariate[[n]][[2]] 
    }
    
    # Extract univariate volatilities
    mSigma = matrix(0,nrow = nrow(mData), ncol = ncol(mY))
    
    for(n in 1:ncol(mY)) {
        mSigma[,n] <- lFit_univariate[[n]][[2]]
    }
    
    # Extract univariate estimated parameters
    mCoef = matrix(0,nrow = ncol(mData), ncol = 3)
    
    for(n in 1:ncol(mY)) {
        mCoef[n,] <- lFit_univariate[[n]]$vPar
    }
    
    # Likelihood of the volatility  part
    dLLK_V = 0 
    
    for(n in 1:ncol(mY)) {
        dLLK_V <- dLLK_V + lFit_univariate[[n]][[4]]   #(Eq.24)
    }
    
    # Correlation and Covariances
    iT = nrow(mY)
    iN = ncol(mEta)
    
    mR = cor(mEta)                  #(Eq.15)
    aS = array(NaN, dim = c(2,2,iT))
    
    for (t in 1:iT) {
        mD = diag(c(as.numeric(mSigma[t,1]), as.numeric(mSigma[t,2])), ncol = 2) # (Eq. 14)
        aS[,,t] = mD %*% mR %*% mD #(Eq.13)
    }
    
    
    # Likelihood of the correlation part
    dLLK_C = 0
    
    for (t in 1:iT) {
        
        dLLK_C      = dLLK_C + mEta[t,, drop = FALSE] %*% solve(mR) %*% t(mEta[t,, drop = FALSE]) - 
            mEta[t,, drop = FALSE] %*% t(mEta[t,, drop = FALSE]) + log(det(mR)) # (Eq.25)
    }
    
    dLLK_C = -0.5 * dLLK_C
    
    # Total likelihood
    dLLK = dLLK_V + dLLK_C # (Eq.23)
    
    lOut = list()
    lOut[["dLLK"]]    = dLLK
    lOut[["BIC"]]     = log(iT) * 8 - 2 * dLLK
    lOut[["mCoef"]]   = mCoef
    lOut[["mSigma"]]  = mSigma
    lOut[["mCor"]]    = mR
    lOut[["aS"]]      = aS
    lOut[["mEta"]]    = mEta
    
    return(lOut)
}

# Filter of the Dynamic Conditional Correlation (DCC) Model
# Input:        Standardized returns, coefficients; alpha and beta 
#               and unconditional correlation of the standardized returns.
# Output:       Loglikelihood and the correlation matrix.    
DCC_Filter <- function(mEta, dA, dB, mQ) {
                                  #{Q_bar}
    
    iN = ncol(mEta)
    iT = nrow(mEta)
    
    # Initialize the array for the correlations
    aCor = array(0, dim = c(iN, iN, iT))
    
    # Initialize the array for the Q matrices
    aQ = array(0, dim = c(iN, iN, iT))
    
    # Initialise at the unconditional correlation
    aCor[,,1]   = mQ
    aQ[,,1]     = mQ
    
    # Initialise the likelihood.
    dLLK = mEta[1, , drop = FALSE] %*% solve(aCor[,, 1]) %*% t(mEta[1, , drop = FALSE]) - 
        mEta[1, , drop = FALSE] %*% t(mEta[1, , drop = FALSE]) + log(det(aCor[,, 1])) # Initial part(without sum), of (Eq.25)
    
    for (t in 2:iT) {  # (Eq.18)
        aQ[,,t]   = mQ * (1 - dA - dB) + dA * t(mEta[t - 1, , drop = FALSE]) %*% mEta[t - 1, , drop = FALSE] +  dB * aQ[,, t - 1] 
              #  {Q_bar}
              # Compute the correlation as Q_tilde^{-1/2} Q Q_tilde^{-1/2}  (Eq.17)
        aCor[,, t] = diag(sqrt(1/diag(aQ[,, t]))) %*% aQ[,, t] %*% diag(sqrt(1/diag(aQ[,, t]))) 
              #  
        dLLK      = dLLK + mEta[t,, drop = FALSE] %*% solve(aCor[,, t]) %*% t(mEta[t,, drop = FALSE]) - 
            mEta[t,, drop = FALSE] %*% t(mEta[t,, drop = FALSE]) + log(det(aCor[,, t]))
    }
    
    lOut = list()
    lOut[["dLLK"]] = -0.5 * dLLK
    lOut[["aCor"]] = aCor
    
    return(lOut)
}

# Estimation of the Dynamic Conditional Correlation (DCC) Model
# Input:        Vector of observations.
# Output:       Optimized parameter coefficients and volatility
#               for the marginal models, parameters, correlation (constant) and
#               covariance matrix for the DCC model, standardized returns,
#               Loglikelihood as well as BIC.
Estimate_DCC <- function(mY) {
    
    # List where marginal models are stored
    lFit_univariate = list()
    
    # Estimate the univariate GARCH models
    for(n in 1:ncol(mY)) {
        lFit_univariate[[n]] = GARCH_Estimation(mY[, n])
    }
    
    # Standardized returns
    mEta = matrix(0,nrow = nrow(mData), ncol = ncol(mY))
    
    for(n in 1:ncol(mY)) {
        mEta[,n] <- mData[,n] / lFit_univariate[[n]][[2]]
    }
    
    ###### Maximization of the DCC likelihood ######
    
    # Initial parameters
    vPar = c(0.04, 0.9)
    
    # Unconditional correlation
    mQ = cor(mEta)
    
    # Maximize the DCC likelihood
    optimizer = solnp(vPar, fun = function(vPar, mEta, mQ) {
        
        Filter = DCC_Filter(mEta, vPar[1], vPar[2], mQ)
        dNLLK = -as.numeric(Filter$dLLK)
        
        return(dNLLK)
    }, ineqfun = function(vPar, ...) {
        sum(vPar)
    }, ineqLB = 1e-4, ineqUB = 0.999, LB = c(1e-4, 1e-4), UB = c(0.999, 0.999), mEta = mEta, mQ = mQ)
    
    # Extract the estimated parameters
    vPar = optimizer$pars
    
    # Extract the likelihood of the correlation part
    dLLK_C = -tail(optimizer$values, 1)
    
    # Filter the dynamic correlation using the estimated parameters
    Filter = DCC_Filter(mEta, vPar[1], vPar[2], mQ)
    
    # Extract univariate volatilities
    mSigma = matrix(0,nrow = nrow(mData), ncol = ncol(mY))
    
    for(n in 1:ncol(mY)) {
        mSigma[,n] <- lFit_univariate[[n]][[2]]
    }
    
    # Extract univariate estimated parameters
    mCoef = matrix(0,nrow = ncol(mData), ncol = 3)
    
    for(n in 1:ncol(mY)) {
        mCoef[n,] <- lFit_univariate[[n]]$vPar
    }
    
    # Compute the likelihood of the volatility  part
    dLLK_V = 0 
    
    for(n in 1:ncol(mY)) {
        dLLK_V <- dLLK_V + lFit_univariate[[n]][[4]]
    }
    
    # Compute the total likelihood
    dLLK = dLLK_V + dLLK_C
    
    # Compute the covariances
    aCor = Filter[["aCor"]]
    iT = nrow(mY)
    
    aS = array(NaN, dim = c(2,2,iT))
    
    for (t in 1:iT) {
        mD = diag(c(as.numeric(mSigma[t,1]), as.numeric(mSigma[t,2])), ncol = 2)
        aS[,,t] = mD %*% aCor[,,t] %*% mD
    }
    
    lOut = list()
    lOut[["dLLK"]]    = dLLK
    lOut[["mCoef"]]   = mCoef
    lOut[["vPar"]]    = vPar
    lOut[["mSigma"]]  = mSigma
    lOut[["aCor"]]    = aCor
    lOut[["aS"]]      = aS
    lOut[["mEta"]]    = mEta
    lOut[["BIC"]]     = log(iT) * 8 - 2 * dLLK
    
    return(lOut)
}

# Estimation of the Minimum Variance Portfolio
# Input:        Covariance matrix, Sigma from our DCC and CCC models.
# Output:       Array of weights, omega and 1-omega, in two assets to form MVP.
MVP <- function(aSigma) {
    
    iT      = length(as.numeric(aSigma[1,1,]))
    aOmega  = array(NaN, dim = c(2,1,iT))
    vOnes   = matrix(1, nrow = ncol(aSigma))
    
    for (t in 1:iT) {
        
        mSigma = aSigma[,,t]
        mSigmaInv = solve(mSigma)
        
        aOmega[,,t] = mSigmaInv %*% vOnes
        aOmega[,,t] = aOmega[,,t] / sum(aOmega[,,t])
    }
    
    return(aOmega)
}

############################################################################
######################## Question 2: Emperical part ########################
############################################################################

# Marginal Densities - GARCH(1,1) models
GSPC_Fit_GARCH = GARCH_Estimation(mData[,1])
DJI_Fit_GARCH  = GARCH_Estimation(mData[,2])

# Estimation of our DCC and CCC models based on data
Fit_DCC   = Estimate_DCC(mData)
Fit_CCC   = Estimate_CCC(mData)

# Model Comparison (DCC vs CCC) - Correlations
par(mfrow=c(1,1))
par(mar = c(5,5,3,3))
plot(vDates_GSPC, Fit_DCC$aCor[1,2,], main = "Correlation of DCC and CCC models", 
     ylab = "Correlation", xlab = "Time", type = "l", ylim = c(0.75,1))
lines(vDates_GSPC, rep(Fit_CCC$mCor[1,2], length(vDates_GSPC)), col = "red")
legend("bottomleft", legend = c("DCC", "CCC"), col = c("black", "red"), lty = 1:1, cex = 0.5)

#######################################
##### Minimum Variance Portfolios #####
#######################################

# Estimation of our Minimum Variance Portfolios
MVP_CCC <- MVP(Fit_CCC$aS)
MVP_DCC <- MVP(Fit_DCC$aS)

# Plot of MVP - Weight in GSPC
plot(vDates_GSPC, MVP_DCC[1,1,], main = "MVP - Weight on GSPC",
     ylab = "Weight", xlab = "Time", type = "l", ylim = c(-3,3))
lines(vDates_GSPC, MVP_CCC[1,1,], col = "red")
abline(a = 0, b = 0, col = "grey", lty = "dashed")
legend("bottomright", legend = c("DCC", "CCC"), col = c("black", "red"), lty = 1:1, cex = 0.5)

# Plot of MVP - Weight in DJI
plot(vDates_GSPC, MVP_DCC[2,1,], main = "MVP - Weight on DJI",
     ylab = "Weight", xlab = "Time", type = "l", ylim = c(-2.5,4))
lines(vDates_GSPC, MVP_CCC[2,1,], col = "red")
abline(a = 0, b = 0, col = "grey", lty = "dashed")
legend("bottomright", legend = c("DCC", "CCC"), col = c("black", "red"), lty = 1:1, cex = 0.5)

#######################################
################ CoVaR ################
#######################################
iT = nrow(mData)
vAlpha = c(0.01, 0.05)

# Make room for VaR for DJI for both values of alpha
vVaR_AlfaLow    = numeric(iT)
vVaR_AlfaHigh   = numeric(iT)

# Make room for thee computed CoVaR
CoVaR_CCC = matrix(0,iT,length(vAlpha))
CoVaR_DCC = matrix(0,iT,length(vAlpha))

for (t in 1:iT) {
    
    # vAlpha = 0.01
    vVaR_AlfaLow[t]   = uniroot(f = function(dVaR){
        pnorm(dVaR, 0, DJI_Fit_GARCH$vSigma[t]) - vAlpha[1]
    }, interval   = c(-100,0))$root
    
    CoVaR_CCC[t,1] = uniroot(f = function(dCoVar) {
        pmvnorm(lower = c(-Inf, -Inf), upper = c(dCoVar, vVaR_AlfaLow[t]), mean = rep(0, length = ncol(mData)), sigma = Fit_CCC$aS[,,t]) - vAlpha[1]^2
    }, interval = c(-100, 100))$root
    
    CoVaR_DCC[t,1] = uniroot(f = function(dCoVar) {
        pmvnorm(lower = c(-Inf, -Inf), upper = c(dCoVar, vVaR_AlfaLow[t]), mean = rep(0, length = ncol(mData)), sigma = Fit_DCC$aS[,,t]) - vAlpha[1]^2
    }, interval = c(-100, 100))$root
    
    # vAlpha = 0.05
    vVaR_AlfaHigh[t]   = uniroot(f = function(dVaR){
        pnorm(dVaR, 0, DJI_Fit_GARCH$vSigma[t]) - vAlpha[2]
    }, interval   = c(-100,0))$root
    
    CoVaR_CCC[t,2] = uniroot(f = function(dCoVar) {
        pmvnorm(lower = c(-Inf, -Inf), upper = c(dCoVar, vVaR_AlfaHigh[t]), mean = rep(0, length = ncol(mData)), sigma = Fit_CCC$aS[,,t]) - vAlpha[2]^2
    }, interval = c(-100, 100))$root
    
    CoVaR_DCC[t,2] = uniroot(f = function(dCoVar) {
        pmvnorm(lower = c(-Inf, -Inf), upper = c(dCoVar, vVaR_AlfaHigh[t]), mean = rep(0, length = ncol(mData)), sigma = Fit_DCC$aS[,,t]) - vAlpha[2]^2
    }, interval = c(-100, 100))$root
    
}

# Comparing values of alpha
par(mfrow = c(2,1))
par(mar = c(3,3,2,2))
plot(vDates_GSPC,CoVaR_DCC[,1], main = expression(paste("CoVaR of GSPC given DJI, ", alpha, "= 0.01 vs. ", alpha, "= 0.05")), type = "l")
lines(vDates_GSPC, CoVaR_DCC[,2], col = "red")
legend("bottomright", legend = c(expression(paste("DCC, ", alpha, " = 0.01")), expression(paste("DCC, ", alpha, " = 0.05"))), col = c("black", "red"), 
       lty = 1:1, cex = 0.5)

plot(vDates_GSPC,CoVaR_CCC[,1], main = expression(paste("CoVaR of GSPC given DJI, ", alpha, "= 0.01 vs. ", alpha, "= 0.05")), type = "l")
lines(vDates_GSPC, CoVaR_CCC[,2], col = "red")
legend("bottomright", legend = c(expression(paste("CCC, ", alpha, " = 0.01")), expression(paste("CCC, ", alpha, " = 0.05"))), col = c("black", "red"), 
       lty = 1:1, cex = 0.5)

# Comparing models - NOT FEATURED IN THE REPORT
plot(vDates_GSPC,CoVaR_DCC[,1], main = expression(paste("CoVaR of GSPC given DJI, ", alpha, "= 0.01")), type = "l")
lines(vDates_GSPC, CoVaR_CCC[,1], col = "red")
legend("bottomright", legend = c("DCC", "CCC"), col = c("black", "red"), lty = 1:1, cex = 0.5)

plot(vDates_GSPC, CoVaR_DCC[,2], main = expression(paste("CoVaR of GSPC given DJI, ", alpha, "= 0.05")), type = "l")
lines(vDates_GSPC, CoVaR_CCC[,2], col = "red")
legend("bottomright", legend = c("DCC", "CCC"), col = c("black", "red"), lty = 1:1, cex = 0.5)


#######################################
# CoVaR= - NOT FEATURED IN THE REPORT #
#######################################
iT = nrow(mData)

CoVaR_equal_DCC <- matrix(NA, iT, 2)
CoVaR_equal_CCC <- matrix(NA, iT, 2)

for (t in 1:iT) {
    CoVaR_equal_DCC[t,1] = GSPC_Fit_GARCH$vSigma[t] * (Fit_DCC$aCor[1,2,t] * qnorm(vAlpha[1], mean = 0, sd = 1) +
                                                           qnorm(vAlpha[1], mean = 0, sd = 1) * sqrt(1-Fit_DCC$aCor[1,2,t]^2))
    
    CoVaR_equal_CCC[t,1] = GSPC_Fit_GARCH$vSigma[t] * (Fit_CCC$mCor[1,2] * qnorm(vAlpha[1], mean = 0, sd = 1) + 
                                                           qnorm(vAlpha[1], mean = 0, sd = 1) * sqrt(1-Fit_CCC$mCor[1,2]^2))
    
    CoVaR_equal_DCC[t,2] = GSPC_Fit_GARCH$vSigma[t] * (Fit_DCC$aCor[1,2,t] * qnorm(vAlpha[2], mean = 0, sd = 1) +
                                                           qnorm(vAlpha[2], mean = 0, sd = 1) * sqrt(1-Fit_DCC$aCor[1,2,t]^2))
    
    CoVaR_equal_CCC[t,2] = GSPC_Fit_GARCH$vSigma[t] * (Fit_CCC$mCor[1,2] * qnorm(vAlpha[2], mean = 0, sd = 1) + 
                                                           qnorm(vAlpha[2], mean = 0, sd = 1) * sqrt(1-Fit_CCC$mCor[1,2]^2))
}

# Comparing values of alpha
plot(vDates_GSPC, CoVaR_equal_DCC[,1], main = expression(paste("CoVaR= of GSPC given DJI, ", alpha, "= 0.01")), type = "l")
lines(vDates_GSPC, CoVaR_equal_DCC[,2], col = "red")
legend("bottomright", legend = c(expression(paste("DCC, ", alpha, " = 0.01")), expression(paste("DCC, ", alpha, " = 0.05"))), 
       col = c("black", "red"), lty = 1:1, cex = 0.5)

plot(vDates_GSPC, CoVaR_equal_CCC[,1], main = expression(paste("CoVaR= of GSPC given DJI, ", alpha, "= 0.05")), type = "l")
lines(vDates_GSPC, CoVaR_equal_CCC[,2], col = "red")
legend("bottomright", legend = c(expression(paste("CCC, ", alpha, " = 0.01")), expression(paste("CCC, ", alpha, " = 0.05"))), 
       col = c("black", "red"), lty = 1:1, cex = 0.5)

