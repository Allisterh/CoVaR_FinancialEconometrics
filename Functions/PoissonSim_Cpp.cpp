// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

//[[Rcpp::export]]
Rcpp::List PoissonSim_Cpp(int iT, double dPhi, double dAlpha) {
  vec vY(iT);
  vec vLambda(iT);
  vLambda(0) =  dPhi/(1-dAlpha);
  vY(0) = Rf_rpois(vLambda(0));
  for (int t = 1; t < iT; t++) {
    vLambda(t) = dPhi + dAlpha * vLambda(t-1);
    vY(t) = Rf_rpois(vLambda(t));
  }
  List lOut;
  lOut["vLambda"] = vLambda;
  lOut["vY"] = vY;
  return lOut;
} 