##### Filter for the Patton (2006) model #####
# mU is a iT x iN matrix with estimated PIT
# CopType is either "t" or "norm" for the t and normal copula, respectively.
# dOmeha, dAlpha, dBeta are the parameters driving the 
# copula correlation dynamic. dNu is the degree of freedom patameter
# in the t copula case. 

####################################################################################################
# 1. Defines a function called PattonFilter that takes as input a matrix of uniformly distributed
#    pseudo-random numbers (mU), the copula type (CopType), the parameters of the Patton's model 
#    (dOmega, dAlpha, dBeta, dNu) and returns the likelihood contribution (dLLK) and the
#    correlation dynamic (vCor).
# 2. Defines the Lambda mapping function (LambdaFun) that takes as input a vector of real 
#    numbers (x) and returns the vector of real numbers defined as 
#    \Lambda(x)=\frac{1-e^{-x}}{1+e^{-x}}.
# 3. Computes the Probability Integral Transformation (mU_Tr) of mU.
# 4. Computes the number of observations (iT).
# 5. Initializes the correlation dynamic (vCor) and the likelihood contribution (dLLK).
# 6. Computes the first value of the correlation dynamic (vCor[1]) and the first 
#    likelihood contribution (dLLK).
# 7. Defines the main loop that computes the correlation dynamic (vCor) and the 
#    likelihood contribution (dLLK) at each time t.
# 8. Computes the correlation (vCor[t]) using the Patton's (2006) recursion.
# 9. Computes the likelihood contribution (dLLK) at time t.
# 10.Returns the likelihood contribution (dLLK) and the correlation dynamic (vCor).
####################################################################################################



PattonFilter <- function(mU, CopType, dOmega, dAlpha, dBeta, dNu, M = 1) {
  
  # Lambda mapping function (modified logistic function)
  # Returns \Lambda(x)=\frac{1-e^{-x}}{1+e^{-x}}, see slide 32 of lecture 11
  LambdaFun <- function(x) {
    (1 - exp(-x))/(1 + exp(-x))
  }
  
  # Probability Integral Transformation
  if (CopType == "norm") {
    mU_Tr = qnorm(mU)
  } 
  if (CopType == "t") {
    mU_Tr = qt(mU, dNu)
  }
  
  # Number of observations
  iT = nrow(mU)

  # Initialize the correlation dynamic
  vCor = numeric(iT) # Vector of correlation coefficients
  vCor[1] = cor(mU_Tr[, 1], mU_Tr[, 2]) # Empirical correlation
  
  # Compute the first likelihood contribution
  if (CopType == "norm") {
    norm.cop <- normalCopula(vCor[1]) # Normal copula from package copula
    dLLK = dCopula(mU[1, ], norm.cop, log = TRUE) # Copula density
  }
  if (CopType == "t") {
    t.cop = tCopula(vCor[1], df = dNu) # t copula from package copula
    dLLK = dCopula(mU[1, ], t.cop, log = TRUE) # Copula density
  }
  
  ## Main loop ##
  for (t in 2:iT) {
    for (j in 1:M) {
      # Compute the correlation using the Patton's (2006) recursion, see slide 32 of lecture 11
      vCor[t] = LambdaFun(dOmega + dBeta * vCor[t - 1] + dAlpha * sum(mU_Tr[t - j, 1] * mU_Tr[t - j, 2]))
    }
    # Update the correlation using the Patton's (2006) recursion, see slide 32 of lecture 11
    #vCor[t] = LambdaFun(dOmega + dBeta * vCor[t - 1] + dAlpha * mU_Tr[t - 1, 1] * mU_Tr[t - 1, 2])
    
    # Compute the likelihood contribution at time t
    if (CopType == "norm") {
      norm.cop <- normalCopula(vCor[t])
      dLLK = dLLK + dCopula(mU[t, ], norm.cop, log = TRUE)
    }
    if (CopType == "t") {
      t.cop = tCopula(vCor[t], df = dNu)
      dLLK = dLLK + dCopula(mU[t, ], t.cop, log = TRUE)
    }
  }
  
  # Output the result
  lOut = list()
  lOut[["dLLK"]] = dLLK
  lOut[["vCor"]] = vCor
  
  return(lOut)
}

##### Main function to estimate the Patton model with t or normal copula #####
### mY are the observations. The univariate models are GARCH with t distributed errors ###

####################################################################################################
# 1. The function "Estimate_Patton" estimates the Patton model for the returns of the portfolio.
#    The function takes the portfolio returns "mY" and the copula type "CopType" as arguments.
#    The function returns a list of the estimated parameters, the log likelihood, and the BIC. 
# 2. The function "ugarchspec" specifies the univariate GARCH model. We set the mean model
#    to be ARMA(0,0) and the distribution model to be "std" (standard normal).
# 3. The function "ugarchfit" estimates the univariate GARCH models using the returns 
#    of the portfolio "mY".
# 4. The function "pit" computes the probabilities of the outcomes of the univariate GARCH models.
# 5. The function "cor" computes the correlation between the returns of the two assets. 
# 6. The function "solnp" is a solver for non-linear programming problems. 
#    We use the solver to maximize the log likelihood of the Patton model.
#    We impose the constraint alpha + beta < 1 to avoid explosive patterns in the 
#    correlation parameter. 
# 7. The function "PattonFilter" computes the filtered correlation parameter.
####################################################################################################

Estimate_Patton <- function(mY, CopType, M = 1) {
  
  ## estimate the marginal GARCH models
  require(rugarch) #load the rugarch package
  require(copula)
  require(Rsolnp)
  
  # Specify the univariate GARCH models
  SpecGARCH = ugarchspec(mean.model = list(armaOrder = c(0, 0)), distribution.model = "std")
  
  # Estimate the univariate GARCH models
  lFit_univariate = list()
  
  # Run ugarchfit for each column of mY.
  for(n in 1:ncol(mY)) {
    lFit_univariate[[n]] = ugarchfit(SpecGARCH, mY[, n])
  }
  
  ## Extract the Probability Integral Transformation
  # 1. Get a list of univariate linear models (lFit_univariate)
  # 2. For each model in the list, extract the probabilities of each of the 
  #    outcomes (pit(Fit))
  # 3. Convert the list of probabilities into a matrix (do.call(cbind, ...))
  mU = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(pit(Fit))
  }))
  
  ## Robustify
  mU[mU > 0.999] = 0.999
  mU[mU < 0.001] = 0.001
  
  ## Maximization of the Copula likelihood in the t and norm copula cases
  # here we impose the constraint alpha + beta < 1 to avoid explosive patterns
  # in the correlation parameter
  if (CopType == "norm") { # Normal copula case
    # Approximated unconditional correlation
    dCor_app = cor(mU[, 1], mU[, 2]) * 0.16 
    # Approximated unmapped unconditional correlation
    dOmega_starting = log(dCor_app + 1) - log(1 - dCor_app) 
    vPar = c(dOmega_starting, 0.04, 0.8)
    
    optimizer = solnp(vPar, fun = function(vPar, mU) {
      
      Filter = PattonFilter(mU, CopType = "norm", M = M, vPar[1], vPar[2], vPar[3], dNu = NA)
      dNLLK = -as.numeric(Filter$dLLK)
      
      if (!is.finite(dNLLK)) {
        dNLLK = 1e4
      }
      
      if (!is.numeric(dNLLK)) {
        dNLLK = 1e4
      }

      return(dNLLK)
    }, 
    LB = c(-3, -0.999, 1e-4), UB = c(3, 0.999, 0.9999), 
    mU = mU)
  }
  
  if (CopType == "t") { # t copula case
    ## Unmap initial value
    # Approximated unconditional correlation
    dCor_app = cor(mU[, 1], mU[, 2]) * 0.16
    # Approximated unmapped unconditional correlation
    dOmega_starting = log(dCor_app + 1) - log(1 - dCor_app) 
    vPar = c(dOmega_starting, 0.04, 0.8, 5)
    
    optimizer = solnp(vPar, fun = function(vPar, mU) {
      
      Filter = PattonFilter(mU, CopType = "t", M = M, vPar[1], vPar[2], vPar[3], dNu = vPar[4])
      dNLLK = -as.numeric(Filter$dLLK)
      
      if (!is.finite(dNLLK)) {
        dNLLK = 1e4
      }
      
      if (!is.numeric(dNLLK)) {
        dNLLK = 1e4
      }
      
      return(dNLLK)
      
    },  
    LB = c(-3, -0.999, 1e-4, 2.01), UB = c(3, 0.999, 0.9999, 30), 
    mU = mU)
  }
  
  vPar = optimizer$pars
  dLLK_C = -tail(optimizer$values, 1)
  
  # compute the filtered correlation parameter
  if (CopType == "norm") {
    Filter = PattonFilter(mU, CopType = "norm", M = M, vPar[1], vPar[2], vPar[3], dNu = NA)
  }
  if (CopType == "t") {
    Filter = PattonFilter(mU, CopType = "t", M = M, vPar[1], vPar[2], vPar[3], dNu = vPar[4])
  }
  
  #extract univariate volatilities
  mSigma = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(sigma(Fit))
  }))
  
  #extract univariate estimated parameters
  mCoef = do.call(cbind, lapply(lFit_univariate, function(Fit) {
    as.numeric(coef(Fit))
  }))
  
  #compute the likelihood od the univariate models
  dLLK_V = do.call(sum, lapply(lFit_univariate, function(Fit) {
    as.numeric(likelihood(Fit))
  }))
  
  #compute the total likelihood of the model
  dLLK = dLLK_V + dLLK_C
  
  if (CopType == "norm") {
    iK = 9 # 3 pars for each marginal + 3 for the copula
  }
  if (CopType == "t") {
    iK = 10 # 3 pars for each marginal + 4 for the copula
  }
  
  iT = nrow(mY)
  
  BIC = log(iT) * iK - 2 * dLLK
  
  lOut = list()
  
  #output the results
  lOut[["dLLK"]] = dLLK
  lOut[["mCoef"]] = mCoef
  lOut[["vPar"]] = vPar
  lOut[["mSigma"]] = mSigma
  lOut[["vCor"]] = Filter[["vCor"]]
  lOut[["BIC"]] = BIC
  
  return(lOut)
}