#include <pcadapt/bed-acc.h>
#include <pcadapt/mat-acc.h>

template <class C>
ListOf<NumericVector> nb_nona(C macc, 
                              const LogicalVector &pass, 
                              int sum_pass) {
  
  int n = macc.nrow();
  int m = macc.ncol();
  
  IntegerVector n_nona(m, n);
  IntegerVector m_nona(n, sum_pass);
  //IntegerVector m_nona(n, m);
  int i, j;
  
  // WARNING: do not use std::size_t because of `n - 4`
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      if (macc(i, j) == 3) {
        
        n_nona[j]--;
        //m_nona[i]--;  
        if (pass[j]) {
          m_nona[i]--;  
        }
        
      }
    }
  }
  
  return List::create(m_nona, n_nona);
}

// [[Rcpp::export]]
ListOf<NumericVector> nb_nona(SEXP obj,
                              const NumericMatrix &lookup_scale,
                              const IntegerMatrix &lookup_byte,
                              const LogicalVector &pass,
                              int sum_pass) {
  
  if (Rf_isMatrix(obj)) {
    matAcc macc(obj, lookup_scale);
    return nb_nona(macc, pass, sum_pass);
  } else {
    XPtr<bed> xpMat(obj);
    bedAcc macc(xpMat, lookup_scale, lookup_byte);
    return nb_nona(macc, pass, sum_pass);
  }
}