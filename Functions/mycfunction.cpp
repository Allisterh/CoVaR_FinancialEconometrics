// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
mat FunC(mat mX, mat mY) {
  mat mZ = mX * mY;
  mat mZInv = mZ.i();
  return mZInv;
}

