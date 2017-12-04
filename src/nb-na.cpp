/******************************************************************************/

#include <pcadapt/bed-acc.h>
#include <pcadapt/mat-acc.h>

/******************************************************************************/

template <class C>
ListOf<NumericVector> nb_nona(C macc) {
  
  int n = macc.nrow();
  int p = macc.ncol();
  int i, j;
  
  IntegerVector n_nona(p, n);
  IntegerVector p_nona(n, p);
  
  for (j = 0; j < p; j++) {
    for (i = 0; i < n; i++) {
      if (macc(i, j) == 3) {
        n_nona[j]--;
        p_nona[i]--; 
      }
    }
  }
  
  return List::create(_["p"] = p_nona, _["n"] = n_nona);
}

/******************************************************************************/

// [[Rcpp::export]]
ListOf<NumericVector> nb_nona(SEXP obj,
                              const NumericMatrix& lookup_scale,
                              const IntegerMatrix& lookup_byte,
                              const IntegerVector& ind_col) {

  // Need access NA as 3
  if (Rf_isMatrix(obj)) {
    matAcc macc(obj, lookup_scale, ind_col);
    return nb_nona(macc);
  } else {
    XPtr<bed> xpMat(obj);
    bedAcc macc(xpMat, lookup_scale, lookup_byte, ind_col);
    return nb_nona(macc);
  }
}

/******************************************************************************/