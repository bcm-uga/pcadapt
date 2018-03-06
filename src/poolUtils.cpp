#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix poolCov(NumericMatrix& m) {
  int nPOP = m.ncol();
  int nSNP = m.nrow();
  NumericMatrix cov(nPOP, nPOP);
  for (int i = 0; i < nPOP; i++) {
    for (int j = 0; j < nPOP; j++) {
      for (int k = 0; k < nSNP; k++) {
        cov(i, j) += m(k, i) * m(k, j);
      }
    }
  }
  return cov;
}

