// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------#
//  Problem 1)   Simulating log returns from GJR-GARCH(1,1) distribution  #
//------------------------------------------------------------------------#
//' Simulate data from the GJR-GARCH(1,1) model
//' 
//' This function simulates data from the GJR-GARCH(1,1) model using a C++ function
//' 
//' @param iT Number of time periods to simulate. An integer.
//' @param dOmega True value of \eqn{\omega}. A double
//' @param dAlpha True value of \eqn{\alpha}. A double.
//' @param dGamma True value of \eqn{\gamma}. A double.
//' @param dBeta True value of \eqn{\beta}. A double.
//' @return Returns a list with the simulated paths for \eqn{y} and \eqn{\sigma^2}
//' @export
// [[Rcpp::export]]
List Cgjrgarchsim(int iT, double dOmega, double dAlpha, double dGamma, double dBeta) {
  arma::vec vSigma2 = zeros(iT), vY = zeros(iT);
  int I;
  vSigma2(0) = dOmega / ( 1 - dAlpha - 0.5*dGamma - dBeta);
  vY(0) = Rf_rnorm(0,1) * sqrt(vSigma2(0)); 
                                
  for (int t = 1; t < iT; t++) {
    I = (vY(t-1) < 0) ? 1 : 0;
    vSigma2(t) = dOmega + (dAlpha + dGamma * I ) * pow(vY(t-1), 2) + dBeta * vSigma2(t-1);
    vY(t) =  Rf_rnorm(0,1) * sqrt(vSigma2(t));
  }
  
  return Rcpp::List::create(Rcpp::Named("vY") = vY,
                                  Named("vSigma2") = vSigma2);
}

//' Mapping function for constrained optimization
//'
//' Map a working parameter \eqn{\tilde{\theta}} to natural parameter \eqn{\theta}.
//' For \eqn{[0, \infty]} parameter constraints, we use the exponential transformation
//' \deqn{\theta = exp(\tilde{\theta})} 
//' For \eqn{[L, U]} constraints, we use the modified logistic transformation
//' \deqn{\theta = L + \frac{(U - L)}{1 + exp(-\tilde{\theta})}}
//' @param dParTilde Parameter value to map
//' @param dLbound Lower bound for restricted parameter space 
//' @param dUbound Upper bound for restricted parameter space 
//'
//' @return Double containg the natural parameter value.
//'
//' @examples
//' dParTilde         = 0.5;                
//' dLbound           = 0;                                                 
//' dUbound           = Inf;                                     
//' CMap2Natural(dParTilde, dLbound, dUbound)
//' @export
// [[Rcpp::export]]
double CMap2Natural(double dParTilde, double dLbound, double dUbound) {      //Mapping to constrained parameter space
  double dPar;
  if (dLbound == 0 && dUbound == R_PosInf) {
   dPar = exp(dParTilde);
  } else {
    dPar = dLbound + (dUbound - dLbound)/(1 + exp(-dParTilde));
  }
  return dPar;
}


//' Mapping function for constrained optimization
//'
//' Map a natural parameter \eqn{\theta} to working parameter  \eqn{\tilde{\theta}}. 
//' For \eqn{[0, \infty]} constraints, we use the log transformation
//' \deqn{\tilde{\theta} = log(\theta)}
//' For \eqn{[L, U]} constraints, we use the inverse of the modified logistic transformation
//' \deqn{\tilde{\theta} = log\left(\frac{\theta - L}{U - \theta}\right)}
//' @param dPar Parameter value to map
//' @param dLbound Lower bound for restricted parameter space 
//' @param dUbound Upper bound for restricted parameter space 
//'
//' @return Double containg the working parameter value.
//'
//' @examples
//' dPar              = 0.5;                
//' dLbound           = 0;                                                 
//' dUbound           = Inf;                                     
//' CMap2Working(dPar, dLbound, dUbound)
//' @export
// [[Rcpp::export]]
double CMap2Working(double dPar, double dLbound, double dUbound) {      //Mapping to constrained parameter space
  double dParTilde;
  if (dLbound == 0 && dUbound == R_PosInf) {
    dParTilde = log(dPar);
  } else {
    dParTilde = log((dPar - dLbound)/(dUbound - dPar));
  }
  return dParTilde;
}

//------------------------------------------------------------------------#
//  Problem 2)   Parameter estimation in the GJR-GARCH(1,1) model         #
//------------------------------------------------------------------------#
//' Mapping of model parameters
//' 
//' This function maps all model parameters to the restricted parameter space. The function calls the documented Map2Natural function.
//' @param vPar Vector of input "tilde" parameters (double)
//' @param dLowerLimit Lower bound for restricted parameter space.
//' @param dUpperLimit Upper bound for restricted parameter space.
//' @return Returns a numeric vector of "natural" mapped parameter values 
//' @export
// [[Rcpp::export]]
arma::vec CMap2NaturalAll(arma::vec vPar, double dLowerLimit = 1e-4, double dUpperLimit = 1-1e-4) {
  double dOmega = CMap2Natural(vPar(0), 0, R_PosInf);
  double dAlpha = CMap2Natural(vPar(1), dLowerLimit, dUpperLimit);
  double dGamma = CMap2Natural(vPar(2), dLowerLimit, 2.0 * (dUpperLimit - dAlpha));
  double dBeta =  CMap2Natural(vPar(3), dLowerLimit, dUpperLimit - dAlpha - 0.5 * dGamma);
  arma::vec vParOut = zeros(vPar.size());
  vParOut(0) = dOmega;
  vParOut(1) = dAlpha;
  vParOut(2) = dGamma;
  vParOut(3) = dBeta;
  return vParOut;
}

//' Mapping of model parameters
//' 
//' This function maps all model parameters from the restricted parameter space to working parameters.
//' The function calls the documented Map2Working function from the ExamPackage.
//' @param vPar Vector of input "natural" parameters (double)
//' @param dLowerLimit Lower bound for restricted parameter space.
//' @param dUpperLimit Upper bound for restricted parameter space.
//' @return Returns a numeric vector of "tilde" working parameter values
//' @export
// [[Rcpp::export]]
arma::vec CMap2WorkingAll(arma::vec vPar, double dLowerLimit = 1e-4, double dUpperLimit = 1-1e-4) {
  double dOmega = CMap2Working(vPar(0), 0, R_PosInf);
  double dAlpha = CMap2Working(vPar(1), dLowerLimit, dUpperLimit);
  double dGamma = CMap2Working(vPar(2), dLowerLimit, 2.0 * (dUpperLimit - dAlpha));
  double dBeta =  CMap2Working(vPar(3), dLowerLimit, dUpperLimit - dAlpha - 0.5 * dGamma);
  arma::vec vParOut = zeros(vPar.size());
  vParOut(0) = dOmega;
  vParOut(1) = dAlpha;
  vParOut(2) = dGamma;
  vParOut(3) = dBeta;
  return vParOut;
}

//' Find starting values
//' 
//' This function solves for a reasonable vector of starting values given the simulated data.
//' @param vY Vector of simulated data from the GJR-GARCH(1,1) model
//' @return Returns a numeric vector of starting values
//' @export
// [[Rcpp::export]]
arma::vec CStartingValues(arma::vec vY) {
  double dAlpha = 0.05, dGamma = 0.07, dBeta = 0.90, dOmega = var(vY) * (1.0 - dAlpha - 0.5 * dGamma - dBeta); //Using the expression for the unconditional variance
  arma::vec vPar0 = zeros(4);
  vPar0(0) = dOmega;
  vPar0(1) = dAlpha;
  vPar0(2) = dGamma;
  vPar0(3) = dBeta;
  return vPar0;
}

//' Negative average log-likelihood
//' 
//' This function calculates the negative average log-likelihood for the GJR-GARCH(1,1) model
//' @param vPar Vector of input "natural" parameters (double)
//' @param vY Vector of simulated data from the GJR-GARCH(1,1) model
//' @return Returns a list. The first element is a double containing the average negative log-likelihood. The second element is the filtered variance process.
//' @export
// [[Rcpp::export]]
List CNegAvgLL(arma::vec vY, arma::vec vPar) {
  
  double dOmega = vPar(0), dAlpha = vPar(1), dGamma = vPar(2), dBeta = vPar(3), dLL;
  int I, iT = vY.size();
  arma::vec vSigma2 = zeros(iT);//Initialize length of variance process vector
  
  vSigma2(0) = dOmega / ( 1 - dAlpha - 0.5*dGamma - dBeta); //First draw
  dLL = R::dnorm( vY(0), 0, sqrt(vSigma2(0)), TRUE);//Log-normal density 
  
  for (int t = 1; t < iT; t++) {
    I = (vY(t-1) < 0) ? 1 : 0;
    vSigma2(t) = dOmega + (dAlpha + dGamma * I ) * pow(vY(t-1), 2) + dBeta * vSigma2(t-1);
    dLL += R::dnorm( vY(t), 0, sqrt(vSigma2(t)), TRUE);
  }
  
  return Rcpp::List::create(Rcpp::Named("dLLK") = - dLL / iT,
                                  Named("vSigma2") = vSigma2);
}

//' Likelihood link function
//' 
//' This function reparameterizes the average negative log-likelihood function
//' @param vParTilde Vector of input "working" parameters (double)
//' @param vY Vector of simulated data from the GJR-GARCH(1,1) model
//' @return Returns a list. The first element is a double containing the average negative log-likelihood. The second element is the filtered variance process.
//' @export
// [[Rcpp::export]]
double CNegAvgLLlink(arma::vec vY, arma::vec vParTilde) {
  arma::vec vPar = CMap2NaturalAll(vParTilde);
  List lLLK = CNegAvgLL(vY = vY, vPar = vPar);
  return  lLLK["dLLK"];
}



