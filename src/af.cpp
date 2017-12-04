/******************************************************************************/

#include <pcadapt/bed-acc.h>
#include <pcadapt/mat-acc.h>

/******************************************************************************/

template <class C>
NumericVector AF(C macc) {
  
  size_t n = macc.nrow();
  size_t p = macc.ncol();
  size_t i, j, n_available;
  
  double x;
  NumericVector af(p);
  
  for (j = 0; j < p; j++) {
    n_available = n; // Counts the number of available values for SNP j
    for (i = 0; i < n; i++) {
      x = macc(i, j);
      if (x == 3) { // Checking a 3 is much faster that checking a NA
        n_available--;
      } else {
        af[j] += x;
      }
    }
    af[j] /= 2 * n_available;
  }
  
  return af;
}

/******************************************************************************/

// Dispatch function for af
// [[Rcpp::export]]
NumericVector get_af(SEXP obj,
                     const NumericMatrix& lookup_scale,
                     const IntegerMatrix& lookup_byte,
                     const IntegerVector& ind_col) {
  
  if (Rf_isMatrix(obj)) {
    matAcc macc(obj, lookup_scale, ind_col);
    return AF(macc);
  } else {
    XPtr<bed> xpMat(obj);
    bedAcc macc(xpMat, lookup_scale, lookup_byte, ind_col);
    return AF(macc);
  }
  
}

/******************************************************************************/