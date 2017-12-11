/******************************************************************************/

// [[Rcpp::depends(RcppArmadillo)]]
#include <pcadapt/bed-acc.h>
#include <pcadapt/mat-acc.h>

/******************************************************************************/

// // [[Rcpp::export]]
// void test(arma::uvec x) {
//   Rcout << x.subvec(0, 2) << std::endl;
// }

arma::mat invCrossprodArma(const arma::mat& x) {
  return inv_sympd(x.t() * x);
}

/******************************************************************************/

template <class C>
NumericMatrix multLinReg(C macc, 
                         const arma::mat& u) {
  
  size_t n = macc.nrow();
  size_t p = macc.ncol();
  int K = u.n_cols;
  
  arma::mat K_inv(K, K);
  arma::colvec betas(K);
  NumericMatrix Z(p, K);
  arma::colvec x(n);
  arma::uvec ind_row(n);
  int k, nona;
  
  for (size_t j = 0; j < p; j++) {
    
    nona = 0;
    for (size_t i = 0; i < n; i++) {
      x[i] = macc(i, j);
      if (x[i] != 3) ind_row[nona++] = i;
    }
    arma::uvec ind_row2 = ind_row.subvec(0, nona - 1); 
    arma::colvec x2 = x.elem(ind_row2);
    
    arma::mat u2 = u.rows(ind_row2);
    
    K_inv = invCrossprodArma(u2);
    betas = K_inv * (u2.t() * x2);
    arma::colvec eps = x2 - u2 * betas;
    double tmp = dot(eps, eps) / (nona - K);
    
    for (k = 0; k < K; k++) {
      Z(j, k) = betas(k) / sqrt(K_inv(k, k) * tmp);
    }
  }
  
  return Z;
}

/******************************************************************************/

// Dispatch function for multLinReg
// [[Rcpp::export]]
NumericMatrix multLinReg(SEXP obj,
                         const NumericMatrix& lookup_scale,
                         const IntegerMatrix& lookup_byte,
                         const IntegerVector& ind_col,
                         const arma::mat& u) {
  
  if (Rf_isMatrix(obj)) {
    matAcc macc(obj, lookup_scale, ind_col);
    return multLinReg(macc, u);
  } else {
    XPtr<bed> xpMat(obj);
    bedAcc macc(xpMat, lookup_scale, lookup_byte, ind_col);
    return multLinReg(macc, u);
  }
}

/******************************************************************************/