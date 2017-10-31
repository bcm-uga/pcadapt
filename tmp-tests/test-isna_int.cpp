#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
bool is_na_int(int x) {
  return R_IsNA(x);
}

// [[Rcpp::export]]
bool is_na_int2(int x) {
  return IntegerVector::is_na(x);
}

// [[Rcpp::export]]
bool is_na_int3(int x) {
  return ISNA(x);
}


/*** R
is_na_int(NA_integer_)
is_na_int2(NA_integer_)
is_na_int3(NA_integer_)
*/
