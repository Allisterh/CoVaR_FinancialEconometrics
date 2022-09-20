#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cumsumC(NumericVector x) {
  int n = x.size();
  NumericVector out(n);
  
  for(int i = 0; i < n; i++) {
    out[i] = out[i - 1] + x[i];
  }
  return out;
}

/*** R
vX <- 1:10
cumsumC(vX)
cumsum(vX)
vY <- runif(1e6)
microbenchmark(cumsum(vY), cumsumC(vY), unit = "milliseconds")
 */
