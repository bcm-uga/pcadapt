/******************************************************************************/

#include <pcadapt/bed-acc.h>
#include <pcadapt/mat-acc.h>

/******************************************************************************/

template <class C, class C2>
double total_var_scaled(C macc, C2 macc_scaled) {
  
  int n = macc.nrow();
  int m = macc.ncol();
  
  double res = 0;
  
  for (int j = 0; j < m; j++) {
    
    double col_norm = 0;
    int nona = 0;
    
    for (int i = 0; i < n; i++) {
      if (macc(i, j) != 3) {
        double g_ij = macc_scaled(i, j);
        col_norm += g_ij * g_ij;
        nona++;
      }
    }
    
    res += n * col_norm / nona;
  }
  
  return res;
}

/******************************************************************************/

// [[Rcpp::export]]
double total_var_scaled(SEXP obj,          // af should be ALL allele frequencies
                   const IntegerVector& ind_col,
                   const NumericVector& af,
                   double ploidy) {
  
  if (Rf_isMatrix(obj)) {
    matAcc macc(obj, ind_col);
    matAccScaled macc_scaled(obj, ind_col, af, ploidy, 0);
    return total_var_scaled(macc, macc_scaled);
  } else {
    XPtr<bed> xpMat(obj);
    bedAcc macc(xpMat, ind_col);
    bedAccScaled macc_scaled(xpMat, ind_col, af, ploidy, 0);
    return total_var_scaled(macc, macc_scaled);
  }
}

/******************************************************************************/