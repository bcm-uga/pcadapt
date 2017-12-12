#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix getCode(int NA_VAL = 3) {
  
  IntegerVector num = IntegerVector::create(0, NA_VAL, 1, 2);
  IntegerMatrix code(4, 256);
  
  int i, k, k2;
  int coeff = 1;
  for (i = 0; i < 4; i++) {
    for (k = 0; k < 256; k++) {
      k2 = k / coeff;
      code(i, k) = num[k2 % 4];
    }
    coeff *= 4;
  }
  
  return code;
}

/*** R
all.equal(getCode(), pcadapt:::getCode())
all.equal(getCode(5), pcadapt:::getCode(5))
*/
