// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <pcadapt/bed-acc.h>
#include <pcadapt/mat-acc.h>

// template <class C>
// NumericMatrix multLinReg(C macc, 
//                          const NumericMatrix &u,
//                          const NumericVector &d,
//                          NumericMatrix &v) {
//   
//   size_t n = macc.nrow();
//   size_t p = macc.ncol();
//   size_t K = u.ncol();
//   
//   NumericMatrix Z(p, K);
//   
//   for (size_t j = 0; j < p; j++) {
//     
//     for (size_t k = 0; k < K; k++) {
//       for (size_t i = 0; i < n; i++) {
//         if (macc(i, j) != 3) {
//           v(j, k) += macc(i, j) * u(i, k) / d[k];   
//         }
//       }
//     }
//     
//     double residual = 0;
//     int n_available = 0;
//     NumericVector sum_squared_u(K);  
//     for (size_t i = 0; i < n; i++) {
//       double y = 0;
//       for (size_t k = 0; k < K; k++) {
//         y += u(i, k) * d[k] * v(j, k) ; // Y = UDV
//       }
//       if (macc(i, j) != 3) {
//         residual += (y - macc(i, j)) * (y - macc(i, j)); // macc(i, j) is normalized
//         n_available++;
//         for (int k = 0; k < K; k++) {
//           // sum_squared_u = 1 if SNP j has no missing value
//           // sum_squared_u < 1 if SNP j has at least one missing value
//           sum_squared_u[k] += u(i, k) * u(i, k); 
//         }
//       }
//     }
//     
//     /* t-score */
//     for (size_t k = 0; k < K; k++) {
//       if (residual > 0 && n_available > K) {
//         Z(j, k) = v(j, k) * d[k] / sqrt(residual / (n_available - K));
//       }
//       if (sum_squared_u[k] > 0) {
//         // this should never happen
//         Z(j, k) /= sqrt(sum_squared_u[k]);
//       }
//     }
//   }
//   
//   return Z;
//   
// }

arma::mat invCrossprodArma(const arma::mat& x) {
  return inv_sympd(x.t() * x);
}

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

// [[Rcpp::export]]
void test(arma::uvec x) {
  Rcout << x.subvec(0, 2) << std::endl;
}

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