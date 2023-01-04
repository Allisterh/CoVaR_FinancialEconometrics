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
  return(dTilde)
}


