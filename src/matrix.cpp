#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector cmpt_af_matrix(const NumericMatrix &G) {
  int nIND = G.nrow();
  int nSNP = G.ncol();
  NumericVector af(nSNP);
  for (int j = 0; j < nSNP; j++) {
    int n_available = 0;
    for (int i = 0; i < nIND; i++) {
      if (!NumericVector::is_na(G(i, j))) {
        af[j] += G(i, j);
        n_available++;
      }
    }
    af[j] /= (2 * n_available);
  }
  return(af);
}

//[[Rcpp::export]]
NumericVector prodGx_matrix(const NumericMatrix &G,
                            const NumericVector &x, 
                            const NumericVector &m, 
                            const NumericVector &s,
                            const LogicalVector &pass) {
  
  // Input vector of length p
  // Output vector of length n
  int nIND = G.nrow();
  int nSNP = G.ncol();
  NumericVector y(nIND);
  for (int j = 0; j < nSNP; j++) {
    if (pass[j]) {
      for (int i = 0; i < nIND; i++) {
        if (!NumericVector::is_na(G(i, j))) {
          y[i] += (G(i, j) - 2 * m[j]) * x[j] / s[j];
        }
      }
    }  
  }
  return y;
}

//[[Rcpp::export]]
NumericVector prodtGx_matrix(const NumericMatrix &G,
                             const NumericVector &x, 
                             const NumericVector &m, 
                             const NumericVector &s,
                             const LogicalVector &pass) {
  
  // Input vector of length n
  // Output vector of length p
  int nIND = G.nrow();
  int nSNP = G.ncol();
  NumericVector y(nSNP);
  for (int j = 0; j < nSNP; j++) {
    if (pass[j]) {
      for (int i = 0; i < nIND; i++) {
        if (!NumericVector::is_na(G(i, j))) {
          y[j] += (G(i, j) - 2 * m[j]) * x[i] / s[j];
        }
      }
    }  
  }
  return y;
}