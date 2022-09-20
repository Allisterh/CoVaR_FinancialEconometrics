#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP callRFunction(NumericVector vX, Function f) { // Note the SEXP class for R objects in C++ and the Function class for general functions in C++.
  SEXP res = f(vX);
  return res;
}
