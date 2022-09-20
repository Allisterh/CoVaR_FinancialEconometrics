/*** R
# Your stand-alone C++ file should have extension .cpp and needs to start with:
 */
#include <Rcpp.h>
using namespace Rcpp;

/*** R
# For each function that you want available within R, you need to prefix it with: [[Rcpp::export]]
 */
// [[Rcpp::export]]

double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;
  
  for(int i = 0; i < n; i++) {
    total += x[i];
  }
  return total / n;
}

/*** R
# You can embed R code in special C++ comment blocks.
 
library(microbenchmark)
  vX <- runif(1e5)
  microbenchmark(
    mean(vX),
    meanC(vX), unit = "milliseconds"
  )
  */
