#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix tryMat(int n) {
  NumericMatrix x(n);
  return x;
}

/*** R
tryMat(42)
*/
