#include <pcadapt/bed-acc.h>
#include <pcadapt/mat-acc.h>

template <class C>
NumericMatrix multLinReg(C macc, 
                         const NumericMatrix &u,
                         const NumericVector &d,
                         NumericMatrix &v) {
  
  size_t n = macc.nrow();
  size_t p = macc.ncol();
  size_t K = u.ncol();
  
  NumericMatrix Z(p, K);
  
  for (size_t j = 0; j < p; j++) {
    
    for (size_t k = 0; k < K; k++) {
      for (size_t i = 0; i < n; i++) {
        if (macc(i, j) != 3) {
          v(j, k) += macc(i, j) * u(i, k) / d[k];   
        }
      }
    }
    
    double residual = 0;
    int n_available = 0;
    NumericVector sum_squared_u(K);  
    for (size_t i = 0; i < n; i++) {
      double y = 0;
      for (size_t k = 0; k < K; k++) {
        y += u(i, k) * d[k] * v(j, k) ; // Y = UDV
      }
      if (macc(i, j) != 3) {
        residual += (y - macc(i, j)) * (y - macc(i, j)); // macc(i, j) is normalized
        n_available++;
        for (int k = 0; k < K; k++) {
          // sum_squared_u = 1 if SNP j has no missing value
          // sum_squared_u < 1 if SNP j has at least one missing value
          sum_squared_u[k] += u(i, k) * u(i, k); 
        }
      }
    }
    
    /* t-score */
    for (size_t k = 0; k < K; k++) {
      if (residual > 0 && n_available > K) {
        Z(j, k) = v(j, k) * d[k] / sqrt(residual / (n_available - K));
      }
      if (sum_squared_u[k] > 0) {
        // this should never happen
        Z(j, k) /= sqrt(sum_squared_u[k]);
      }
    }
  }
  
  return Z;
  
}

// Dispatch function for multLinReg
// [[Rcpp::export]]
NumericMatrix multLinReg(SEXP obj,
                         const NumericMatrix &lookup_scale,
                         const IntegerMatrix &lookup_byte,
                         const IntegerVector& ind_col,
                         const NumericMatrix &u,
                         const NumericVector &d,
                         NumericMatrix &v) {
  
  if (Rf_isMatrix(obj)) {
    matAcc macc(obj, lookup_scale, ind_col);
    return multLinReg(macc, u, d, v);
  } else {
    XPtr<bed> xpMat(obj);
    bedAcc macc(xpMat, lookup_scale, lookup_byte, ind_col);
    return multLinReg(macc, u, d, v);
  }
  
}