#include <pcadapt/bed-acc.h>
#include <pcadapt/mat-acc.h>

template <class C>
ListOf<NumericVector> nb_nona(C macc) {
  
  int n = macc.nrow();
  int m = macc.ncol();
  
  IntegerVector n_nona(m, n);
  IntegerVector m_nona(n, m);
  int i, j;
  
  
  // WARNING: do not use std::size_t because of `n - 4`
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      if (macc(i, j) == 3) {
        n_nona[j]--;
        m_nona[i]--;
      }
    }
  }
  
  return List::create(m_nona, n_nona);
}

// [[Rcpp::export]]
ListOf<NumericVector> nb_nona(SEXP obj,
                            const NumericMatrix& lookup_scale,
                            const IntegerMatrix& lookup_byte) {
  
  if (Rf_isMatrix(obj)) {
    matAcc macc(obj, lookup_scale);
    return nb_nona(macc);
  } else {
    XPtr<bed> xpMat(obj);
    bedAcc macc(xpMat, lookup_scale, lookup_byte);
    return nb_nona(macc);
  }
}