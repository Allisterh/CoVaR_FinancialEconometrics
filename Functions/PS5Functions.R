#Needed to make C++ functions work in the package
#' @useDynLib PS5Package
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------------------------------#
#  Problem 1)   Simulating log returns from Dyn-Poi        distribution  #
#------------------------------------------------------------------------#
#' Simulate data from the Dynamic Poisson model
#' 
#' This function simulates data from the dynamic poisson model
#' @param iT Number of time periods to simulate. An integer.
#' @param dAlpha True value of \eqn{\alpha}. A double.
#' @param dPhi True value of \eqn{\gamma}. A double.
#' @return Returns a list with the simulated path for \eqn{y}
#' @export
# b) Simulate from the model
DynPoissonSim = function(iT, dAlpha, dPhi) {
  vY = numeric(iT)
  vY[1] = rpois(1, dPhi / (1-dAlpha) )
  for (i in 2:iT) {
    vY[i] = rpois(1, dPhi + dAlpha * vY[i - 1])
  }
  return(vY)
}

#------------------------------------------------------------------------#
#  Problem 2)   Parameter estimation in the Dyn-Poi model                #
#------------------------------------------------------------------------#
#' Negative average log-likelihood
#' 
#' This function calculates the negative average log-likelihood for the Dyn-Poi model
#' @param vPar Vector of input "natural" parameters. First parameter is \eqn{\phi}. Second parameter is \eqn{\alpha}.
#' @param vY Vector of simulated data from the Dyn-Poi model
#' @return Returns a double containing the average negative log-likelihood.
#' @export
NegAvgll = function(vY, vPar) {
  dPhi = vPar[1]
  dAlpha = vPar[2]
  LL = dpois(vY[1], dPhi / (1-dAlpha), log = TRUE)
  for ( t in 2:length(vY) ) {
    LL = LL + dpois( vY[t], dPhi + dAlpha * vY[t-1], log = TRUE)
  }
  return( -LL / length(vY) )
}

#' Mapping function for constrained optimization
#'
#' Map a working parameter \eqn{\tilde{\theta}} to natural parameter \eqn{\theta}.
#' For \eqn{[0, \infty]} parameter constraints, we use the exponential transformation
#' \deqn{\theta = exp(\tilde{\theta})} 
#' For \eqn{[L, U]} constraints, we use the modified logistic transformation
#' \deqn{\theta = L + \frac{(U - L)}{1 + exp(-\tilde{\theta})}}
#' @param dPar Parameter value to map
#' @param dLbound Lower bound for restricted parameter space 
#' @param dUbound Upper bound for restricted parameter space 
#'
#' @return Double containg the natural parameter value.
#'
#' @examples
#' dPar              = 0.5                
#' dLbound           = 0                                                 
#' dUbound           = Inf                                     
#' Map2Natural(dPar, dLbound, dUbound)
#' @export
Map2Natural = function(dPar, dLbound, dUbound) {      #Mapping to constrained parameter space
  if (dLbound == 0 && dUbound == Inf) {
    return(exp(dPar))
  } else {
    return( dLbound + (dUbound - dLbound)/(1 + exp(-dPar)))
  }
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
#' @return Double containg  the working parameter value.
#'
#' @examples
#' dPar              = 0.5                
#' dLbound           = 0                                                 
#' dUbound           = Inf                                     
#' Map2Working(dPar, dLbound, dUbound) 
#' @export
Map2Working = function(dPar, dLbound, dUbound) {
  if (dLbound == 0 && dUbound == Inf) {
    return( log(dPar) )
  } else {
    return( log((dPar - dLbound)/(dUbound - dPar)) )
  }
}

#' Likelihood link function
#' 
#' This function reparameterizes the average negative log-likelihood function
#' @param vParTilde Vector of input "working" parameters. First parameter is \eqn{\phi}. Second parameter is \eqn{\alpha}.
#' @param vY Vector of simulated data from model
#' @return Returns a double containing the average negative log-likelihood.
#' @export
NegAvgLLlink = function(vParTilde, vY) {
  #Mapping alpha in (0,1) and phi in (0,Inf)
  dPhi = Map2Natural(vParTilde[1],0,Inf)
  dAlpha = Map2Natural(vParTilde[2], 0, 1) 
  return( NegAvgll(vY, c(dPhi, dAlpha) ) )
}

