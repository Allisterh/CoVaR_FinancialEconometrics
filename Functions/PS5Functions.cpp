// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//Problem 1: Simulating data

//' Simulate data from the dynamic poisson model
//' 
//' Function simulates data from the dynamic poisson model for a given set of parameter values
//' 
//' @param dPhi The parameter \eqn{\phi}. A double.
//' @param dAlpha The parameter \eqn{\alpha}. A double.
//' @param iT Number of simulated periods. An integer.
//' @return A vector of a simulated time series from the model.
//' @export
// [[Rcpp::export]]
  arma::vec CDynPoissonSim(int iT, double dPhi, double dAlpha) {

  arma::vec vY = zeros(iT);
  
  //First obs
  vY(0) = Rf_rpois(dPhi / (1 - dAlpha) );
  
  //Subsequent draws
  for (int t = 1; t < iT; t++) {
    vY(t) = Rf_rpois( dPhi + dAlpha * vY(t - 1) );
  }
  
  return vY;
}

//Problem 3: 

//' Negative average log-likelihood for dynamic poisson model
//' 
//' Function calculates the average log-likelihood for the dynamic poisson model
//' 
//' @param vY Is a numeric vector of simulated data from the model
//' @param vParTilde Vector of input "working" parameters. First parameter is \eqn{\phi}. Second parameter is \eqn{\alpha}.
//' @return Returns a double containing the average negative log-likelihood.
//' @export
// [[Rcpp::export]]
double CNegAvgll(arma::vec vY, arma::vec vPar) {
  
  double dPhi = vPar(0);
  double dAlpha = vPar(1);

  int iT = vY.size();
  
  double dLL = Rf_dpois( vY(0), dPhi / (1-dAlpha), 1);
  
  for (int t = 1; t < iT; t++) {
    dLL += Rf_dpois( vY(t),dPhi + dAlpha * vY(t-1), 1 );
  }  
  
  double NegAvgLL = -dLL / iT;
  
  return NegAvgLL;
}