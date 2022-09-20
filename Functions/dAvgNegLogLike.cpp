// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// 3) Write a function in C++ which returns the negative average log–likelihood for a series of observations

//[[Rcpp::export]]
double Neg_Avg_llkC(vec vY, double dPhi, double dAlpha) {

  int iT = vY.size();

  vec vLambda(iT); //intensities
  //initialize lambda
  vLambda(0) = dPhi/(1.0 - dAlpha);

  //initialize likelihood
  double dAvgLogLike = Rf_dpois(vY(0), vLambda(0), 1);

  for (int t = 1; t < iT; t++) {
    vLambda(t) = dPhi + dAlpha * vY(t - 1);
    dAvgLogLike = dAvgLogLike + Rf_dpois(vY(t), vLambda(t), 1);
  }

  return -dAvgLogLike/(iT * 1.0);

}


