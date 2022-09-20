// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
List GarchSim(int iT, double dOmega, double dAlpha, double dBeta) {
  
  vec vY(iT);      // observations. Note: vector declaration with vec
  vec vSigma2(iT); // conditional variances
  
  // initialize at the unconditional value
  vSigma2(0) = dOmega/(1.0 - dAlpha - dBeta);
  
  // sample the first observations
  vY(0) = pow(vSigma2(0), 0.5) * Rf_rnorm(0.0, 1.0); // Note: The use of pow(double, double) for powers and the function Rf rnorm(double, double) for simulation from a Gaussian distribution
  
  for (int t = 1; t<iT; t++) {
    vSigma2(t) =  dOmega + dAlpha * pow(vY(t - 1), 2.0) + dBeta * vSigma2(t - 1);
    vY(t) = pow(vSigma2(t), 0.5) * Rf_rnorm(0.0, 1.0);
  }
  
  List lOut;                  // Note: The List declaration
  lOut["vSigma2"] = vSigma2;  // Note: The list elements assignment
  lOut["vY"] = vY;
  
  return lOut;
}
