// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat matprodC(arma::mat m1, arma::mat m2) {
  if(m1.n_cols != m2.n_rows)
    stop("Incompatible matrix dimensions");
  
  return m1 * m2;
}
