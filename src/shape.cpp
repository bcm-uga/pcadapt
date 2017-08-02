#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat fJ_cpp(int n){
  arma::mat zz(n, n);
  zz.ones();
  zz /= n;
  arma::mat H(n, n);
  H.eye();
  H -= zz;
  return(H);
}

// [[Rcpp::export]]
arma::mat fcnt_cpp(arma::mat &a){
  int nrow = a.n_rows;
  arma::mat aa = fJ_cpp(nrow);
  return(aa * a);
}

//' @export
//'
// [[Rcpp::export]]
arma::mat pca_rotation(arma::mat &a, arma::mat &b){
  arma::mat fcnt_a = fcnt_cpp(a);
  arma::mat fcnt_b = fcnt_cpp(b);
  arma::mat x = fcnt_a.t() * fcnt_b;
  arma::mat u;
  arma::vec s;
  arma::mat v;
  svd(u, s, v, x);
  arma::mat R;
  R = v * u.t();
  arma::cx_vec eigval_v;
  arma::cx_mat eigvec_v;
  arma::cx_vec eigval_u;
  arma::cx_mat eigvec_u;
  eig_gen(eigval_v, eigvec_v, v);
  eig_gen(eigval_u, eigvec_u, u);
  arma::cx_double tmp_v = prod(eigval_v);
  arma::cx_double tmp_u = prod(eigval_u);
  double chk1 = real(tmp_v);
  double chk2 = real(tmp_u);
  if ((chk1 < 0) && (chk2 > 0)) {
    for (int i = 0; i < v.n_rows; i++){
      v(i, v.n_cols - 1) *= (-1);
    }
    R = v * u.t();
  }
  if ((chk2 < 0) && (chk1 > 0)) {
    for (int i = 0; i < u.n_rows; i++){
      u(i, u.n_cols - 1) *= (-1);
    }
    R = v * u.t();
  }
  return(R);
}