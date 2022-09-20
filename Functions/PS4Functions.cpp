// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
//############################################################################
//             11: Simulating compund Poisson variable (Monte Carlo)         #
//############################################################################
//' Monte Carlo simulation
//' 
//' Function simulates data from the dynamic poisson model for a given set of parameter values
//' 
//' @param dLambda Poisson parameter \eqn{\lambda}. A double.
//' @param dMu Log-normal dis. parameter \eqn{\mu}. A double.
//' @param dSigma2 dMu Log-normal dis. parameter \eqn{\sigma^2}. A double.
//' @param iB Number of Monte-Carlo replications. An integer.
//' @return The vector of simulated total payments (vector of compound Poisson variables)
//' @export
// [[Rcpp::export]]
arma::vec MonteCarloC(double dLambda, double dMu, double dSigma2, int iB) { 
  arma::vec vTotalAmount = zeros(iB), vX, dN;
  for (int b = 0; b < iB; b++) {
    dN = Rcpp::rpois(1,dLambda);
    vX =  Rcpp::rlnorm(dN(0,0), dMu, sqrt(dSigma2));
    vTotalAmount(b) = accu(vX);
  }
  return vTotalAmount;
}
