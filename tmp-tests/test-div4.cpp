#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int div4(int x) {
  return x >> 2;
}

// [[Rcpp::export]]
int div4_2(int x) {
  return x / 4;
}

// [[Rcpp::export]]
int mod4(int x) {
  return x || 0x11;
}

// [[Rcpp::export]]
int mod4_2(int x) {
  return x % 4;
}


/*** R
microbenchmark::microbenchmark(
  div4(42),
  div4_2(42)
)
microbenchmark::microbenchmark(
  mod4(42),
  mod4_2(42)
)
*/
