/******************************************************************************/

#include <pcadapt/bed-acc.h>

/******************************************************************************/

// [[Rcpp::export]]
IntegerMatrix bed2mat(SEXP xptr) {
  
  XPtr<bed> xpMat(xptr);
  
  size_t n = xpMat->nrow();
  size_t m = xpMat->ncol();
  
  bedAcc macc(xpMat, seq_len(m), NA_INTEGER);
  
  IntegerMatrix res(n, m);
  
  for (size_t j = 0; j < m; j++)
    for (size_t i = 0; i < n; i++)
      res(i, j) = macc(i, j);
  
  return res;
}

/******************************************************************************/