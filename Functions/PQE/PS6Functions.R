#Needed to make C++ functions work in the package
#' @useDynLib PS6Package
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------------------------------#
#  Problem 1)   Simulating log returns from GJR-GARCH(1,1) distribution  #
#------------------------------------------------------------------------#
#' Simulate data from the GJR-GARCH(1,1) model
#' 
#' This function simulates data from the GJR-GARCH(1,1) model
#' @param iT Number of time periods to simulate. An integer.
#' @param dOmega True value of \eqn{\omega}. A double
#' @param dAlpha True value of \eqn{\alpha}. A double.
#' @param dGamma True value of \eqn{\gamma}. A double.
#' @param dBeta True value of \eqn{\beta}. A double.
#' @return Returns a list with the simulated paths for \eqn{y} and \eqn{\sigma^2}
#' @export
# b) Simulate from the model
gjrgarch.sim = function(iT,dOmega,dAlpha,dGamma,dBeta) {
  vSigma2 = numeric(iT)
  vY = numeric(iT)
  vSigma2[1] = dOmega / ( 1 - dAlpha - 0.5*dGamma - dBeta)
  vY[1] = rnorm(1,0,1) * sqrt(vSigma2[1]) 
  for (t in 2:iT) {
    vSigma2[t] = dOmega + (dAlpha + dGamma * ifelse(vY[t-1] < 0, 1, 0) )*vY[t-1]^2 + dBeta * vSigma2[t-1]
    vY[t] =  rnorm(1,0,1) * sqrt(vSigma2[t]) 
  }
  SimPaths = list(vY = vY,vSigma2 = vSigma2)
  return(SimPaths)
}

#------------------------------------------------------------------------#
#  Problem 2)   Parameter estimation in the GJR-GARCH(1,1) model         #
#------------------------------------------------------------------------#
#' Mapping function for constrained optimization
#'
#' Map a working parameter \eqn{\tilde{\theta}} to natural parameter \eqn{\theta}.
#' For \eqn{[0, \infty]} parameter constraints, we use the exponential transformation
#' \deqn{\theta = exp(\tilde{\theta})} 
#' For \eqn{[L, U]} constraints, we use the modified logistic transformation
#' \deqn{\theta = L + \frac{(U - L)}{1 + exp(-\tilde{\theta})}}
#' @param dParTilde Parameter value to map
#' @param dLbound Lower bound for restricted parameter space 
#' @param dUbound Upper bound for restricted parameter space 
#'
#' @return Double containg the natural parameter value.
#'
#' @examples
#' dParTilde              = 0.5                
#' dLbound           = 0                                                 
#' dUbound           = Inf                                     
#' Map2Natural(dParTilde, dLbound, dUbound)
#' @export
Map2Natural = function(dParTilde, dLbound, dUbound) {      #Mapping to constrained parameter space
  if (dLbound == 0 && dUbound == Inf) {
    dPar = exp(dParTilde)
  } else {
    dPar = dLbound + (dUbound - dLbound)/(1 + exp(-dParTilde))
  }
  return(dPar)
}

#' Mapping function for constrained optimization
#'
#' Map a natural parameter \eqn{\theta} to working parameter  \eqn{\tilde{\theta}}. 
#' For \eqn{[0, \infty]} constraints, we use the log transformation
#' \deqn{\tilde{\theta} = log(\theta)}
#' For \eqn{[L, U]} constraints, we use the inverse of the modified logistic transformation
#' \deqn{\tilde{\theta} = log\left(\frac{\theta - L}{U - \theta}\right)}
#'
#' @param dPar Parameter value to map
#' @param dLbound Lower bound for restricted parameter space 
#' @param dUbound Upper bound for restricted parameter space 
#'
#' @return Double containg the working parameter value.
#'
#' @examples
#' dPar              = 0.5                
#' dLbound           = 0                                                 
#' dUbound           = Inf                                     
#' Map2Working(dPar, dLbound, dUbound) 
#' @export
Map2Working = function(dPar, dLbound, dUbound) {
  if (dLbound == 0 && dUbound == Inf) {
    dParTilde = log(dPar)
  } else {
    dParTilde = log((dPar - dLbound)/(dUbound - dPar))
  }
  return(dParTilde)
}

#' Mapping of model parameters
#' 
#' This function maps all model parameters to the restricted parameter space. The function calls the documented Map2Natural function.
#' @param vPar Vector of input "tilde" parameters (double)
#' @param dLowerLimit Lower bound for restricted parameter space.
#' @param dUpperLimit Upper bound for restricted parameter space.
#' @return Returns a numeric vector of "natural" mapped parameter values 
#' @export
Map2NaturalAll = function(vPar, dLowerLimit = 1e-4, dUpperLimit = 1-1e-4) {
  dOmega = Map2Natural(vPar[1], 0, Inf)
  dAlpha = Map2Natural(vPar[2], dLowerLimit, dUpperLimit)
  dGamma = Map2Natural(vPar[3], dLowerLimit, 2.0 * (dUpperLimit - dAlpha))
  dBeta =  Map2Natural(vPar[4], dLowerLimit, dUpperLimit - dAlpha - 0.5 * dGamma)
  return( c(dOmega, dAlpha, dGamma, dBeta) )
}

#' Mapping of model parameters
#' 
#' This function maps all model parameters from the restricted parameter space to working parameters.
#' The function calls the documented Map2Working function from the ExamPackage.
#' @param vPar Vector of input "natural" parameters (double)
#' @param dLowerLimit Lower bound for restricted parameter space.
#' @param dUpperLimit Upper bound for restricted parameter space.
#' @return Returns a numeric vector of "tilde" working parameter values
#' @export
Map2WorkingAll = function(vPar, dLowerLimit = 1e-4, dUpperLimit = 1-1e-4) {
  dOmegaTilde = Map2Working(vPar[1], 0, Inf)
  dAlphaTilde = Map2Working(vPar[2], dLowerLimit, dUpperLimit)
  dGammaTilde = Map2Working(vPar[3], dLowerLimit, 2.0 * (dUpperLimit - vPar[2]))
  dBetaTilde =  Map2Working(vPar[4], dLowerLimit, dUpperLimit - vPar[2] - 0.5 * vPar[3])
  return( c(dOmegaTilde, dAlphaTilde, dGammaTilde, dBetaTilde) )
}

#' Find starting values
#' 
#' This function solves for a reasonable vector of starting values given the simulated data.
#' @param vY Vector of simulated data from the GJR-GARCH(1,1) model
#' @return Returns a numeric vector of starting values
#' @export
StartingValues = function(vY) {
  dAlpha = 0.05
  dGamma = 0.07
  dBeta = 0.90
  dOmega = var(vY) * (1.0 - dAlpha - 0.5 * dGamma - dBeta)#Using the expression for the unconditional variance
  return(  c(dOmega,dAlpha,dGamma,dBeta) )
}

#' Negative average log-likelihood
#' 
#' This function calculates the negative average log-likelihood for the GJR-GARCH(1,1) model
#' @param vPar Vector of input "natural" parameters
#' @param vY Vector of simulated data from the GJR-GARCH(1,1) model
#' @return Returns a list. The first element is a double containing the average negative log-likelihood. The second element is the filtered variance process.
#' @export
NegAvgLL = function(vY, vPar) {
  dOmega = vPar[1]
  dAlpha = vPar[2]
  dGamma = vPar[3]
  dBeta = vPar[4]
  vSigma2 = numeric(length(vY))#Initialize 
  vSigma2[1] = dOmega / ( 1 - dAlpha - 0.5*dGamma - dBeta)#First draw
  dLL = dnorm(vY[1], 0, sqrt(vSigma2[1]), log = TRUE)
  for (t in 2:length(vY)) {
    vSigma2[t] = dOmega + (dAlpha + dGamma * ifelse(vY[t-1] < 0, 1, 0) )*vY[t-1]^2 + dBeta * vSigma2[t-1]
    dLL = dLL + dnorm(vY[t], 0, sqrt(vSigma2[t]), log = TRUE)
  }
  return(list(dLLK = - dLL / length(vY), vSigma2 = vSigma2))
}

#' Likelihood link function
#' 
#' This function reparameterizes the average negative log-likelihood function
#' @param vParTilde Vector of input "working" parameters (double)
#' @param vY Vector of simulated data from the GJR-GARCH(1,1) model
#' @return Returns a double containing the average negative log-likelihood.
#' @export
NegAvgLLlink = function(vY, vParTilde) {
  vPar = Map2NaturalAll(vParTilde)
  return( NegAvgLL(vY = vY, vPar)$dLLK )
}

#' GJG-GARCH(1,1) estimation
#' 
#' This function estimates the parameters of the GJG-GARCH(1,1) model.
#' @param vY Vector of simulated data from the GJR-GARCH(1,1) model
#' @return Returns a list. The first element is the vector of estimated parameter values. The second element is the filtered variance process.
#' @export
gjrgarch.est = function(vY) {
  vPar0 = StartingValues(vY = vY)   #Starting values
  vPar0Tilde = Map2WorkingAll(vPar0) #Log-transformed starting values
  vParTilde = optim(vPar0Tilde, NegAvgLLlink, vY = vY, method = "BFGS")$par #Estimated parameter vector
  vParEst = Map2NaturalAll(vParTilde) #Transform back estimated parameters
  VarianceProcessEst = NegAvgLL(vY = vY, vParEst)$vSigma2 #Estimated (filtered) variance process
  return(list(vParEst = vParEst, vVarFiltered = VarianceProcessEst)) 
}

#' GJG-GARCH(1,1) estimation using C++
#' 
#' This function estimates the parameters of the GJG-GARCH(1,1) model using C++ functions
#' @param vY Vector of simulated data from the GJR-GARCH(1,1) model
#' @return Returns a list. The first element is the vector of estimated parameter values. The second element is the filtered variance process.
#' @export
Cgjrgarch.est = function(vY) {
  vPar0 = CStartingValues(vY = vY)   #Starting values
  vPar0Tilde = CMap2WorkingAll(vPar0) #Log-transformed starting values
  vParTilde = optim(vPar0Tilde, CNegAvgLLlink, vY = vY, method = "BFGS")$par #Estimated parameter vector
  vParEst = CMap2NaturalAll(vParTilde) #Transform back estimated parameters
  VarianceProcessEst = CNegAvgLL(vY = vY, vParEst)$vSigma2 #Estimated (filtered) variance process
  return(list(vParEst = vParEst, vVarFiltered = VarianceProcessEst)) 
}

#' Plots. An example of a function that we might want to pass to C++.
#' 
#' @param vY Vector to plot.
#' @return Returns NULL as the function is just used as a plotting tool for C++.
#' @export
plotR = function(vY, vSigma2) {
  plot(vY, xlab = "t", ylab = "Log return", main = "Simulated log returns", type = "l")
  plot(vSigma2, xlab = "t", ylab = "Log return variance", main = "Simulated variance process", type = "l")
  abline(h = dOmega / (1-dAlpha-0.5*dGamma-dBeta),col = 'red')
  return(NULL) 
}





