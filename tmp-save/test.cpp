#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector prodGx(const NumericMatrix &G,
                     const NumericVector &x,
                     const NumericVector &p) {
  int nSNP = G.nrow();
  int nIND = G.ncol();
  NumericVector res(nSNP);
  for (int k = 0; k < nSNP; k++) {
    double value = 0;
    int nbmv = 0;
    for (int j = 0; j < nIND; j++) {
      if ((!NumericVector::is_na(G(k, j)))) {
        value += (G(k, j) - 2 * p[k]) * x[j] / (sqrt(2 * p[k] * (1 - p[k])));
      } else {
        nbmv++;
      }
    } 
    res[k] = value / (nIND - nbmv);
  }
  return(res);
}

// [[Rcpp::export]]
NumericVector prodtGx(const NumericMatrix &G,
                      const NumericVector &x,
                      const NumericVector &p) {
  int nSNP = G.nrow();
  int nIND = G.ncol();
  NumericVector res(nIND);
  for (int j = 0; j < nIND; j++) {
    double value = 0;
    int nbmv = 0;
    for (int k = 0; k < nSNP; k++) {
      if ((!NumericVector::is_na(G(k, j)))) {
        value += (G(k, j) - 2 * p[k]) * x[k] / (sqrt(2 * p[k] * (1 - p[k])));  
      } else {
        nbmv++;
      }
    }
    res[j] = value / (nSNP - nbmv);
  }
  return(res);
}