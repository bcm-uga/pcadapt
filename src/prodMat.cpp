#include <pcadapt/bed-acc.h>
#include <pcadapt/mat-acc.h>

/******************************************************************************/

template <class C>
NumericVector pMatVec4(C macc, const NumericVector &x, const IntegerVector &pass) {
  
  int n = macc.nrow();
  int m = macc.ncol();
  
  NumericVector res(n);
  int i, j;
  
  // WARNING: do not use std::size_t because of `m - 4`
  for (j = 0; j <= m - 4; j += 4) { // unrolling optimization
    for (i = 0; i < n; i++) {
      res[i] += (x[j] * macc(i, j) * pass[j] + x[j+1] * macc(i, j+1) * pass[j+1]) +
        (x[j+2] * macc(i, j+2) * pass[j+2] + x[j+3] * macc(i, j+3) * pass[j+3]);
    } // The parentheses are somehow important.
  }
  for (; j < m; j++) {
    for (i = 0; i < n; i++) {
      res[i] += x[j] * macc(i, j) * pass[j];
    }
  }
  
  return res;
}

// [[Rcpp::export]]
NumericVector pMatVec4(SEXP obj, const NumericVector &x,
                       const NumericMatrix &lookup_scale,
                       const IntegerMatrix &lookup_byte,
                       const IntegerVector &pass) {
  
  if (Rf_isMatrix(obj)) {
    matAcc macc(obj, lookup_scale);
    return pMatVec4(macc, x, pass);
  } else {
    XPtr<bed> xpMat(obj);
    bedAcc macc(xpMat, lookup_scale, lookup_byte);
    return pMatVec4(macc, x, pass);
  }
}

/******************************************************************************/

template <class C>
NumericVector cpMatVec4(C macc, const NumericVector &x, const IntegerVector &pass) {
  
  int n = macc.nrow();
  int m = macc.ncol();
  
  NumericVector res(m);
  double tmp;
  int i, j;
  
  // WARNING: do not use std::size_t because of `n - 4`
  for (j = 0; j < m; j++) {
    if (pass[j] == 1) {
      tmp = 0;
      for (i = 0; i <= n - 4; i += 4) { // unrolling optimization
        tmp += (macc(i, j) * x[i] + macc(i+1, j) * x[i+1]) +
          (macc(i+2, j) * x[i+2] + macc(i+3, j) * x[i+3]);
      }
      for (; i < n; i++) {
        tmp += macc(i, j) * x[i];
      }
      res[j] = tmp;
    }
  }
  return res;
}

// [[Rcpp::export]]
NumericVector cpMatVec4(SEXP obj, const NumericVector &x,
                        const NumericMatrix &lookup_scale,
                        const IntegerMatrix &lookup_byte,
                        const IntegerVector &pass) {
  
  if (Rf_isMatrix(obj)) {
    matAcc macc(obj, lookup_scale);
    return cpMatVec4(macc, x, pass);
  } else {
    XPtr<bed> xpMat(obj);
    bedAcc macc(xpMat, lookup_scale, lookup_byte);
    return cpMatVec4(macc, x, pass);
  }
}
