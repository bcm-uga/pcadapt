#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

#define NA 9

using namespace Rcpp;


//' Covariance for loaded genotype data
//' 
//' \code{rsvd_cpp} computes the randomized SVD.
//' 
//' @param A a numeric matrix.
//' @param k an integer.
//' 
//' @return The returned value is a Rcpp::List.
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List rsvd_cpp(arma::mat &A, int k){
  
  // int m = A.n_rows; 
  int n = A.n_cols;
  
  int p = 10;
  int q = 1;
  
  // Oversampling parameter
  int l = k + p;
  if (l > n){
    l = n;
  }
  // int nu = k;
  // int nv = k;
  
  arma::mat O = arma::randu(n, l);
  arma::mat Y = A * O;
  arma::mat Q;
  arma::mat R;
  arma::mat Z;
  
  if (q > 0){
    for (int i = 0; i < q; i++){
      arma::qr_econ(Q, R, Y);
      Y = Q;
      Z = A.t() * Y;
      arma::qr_econ(Q, R, Z);
      Z = Q;
      Y = A * Z;
    }
  }
  arma::qr_econ(Q, R, Y);
  arma::mat B = Q.t() * A;
  return Rcpp::List::create(Rcpp::Named("Q") = Q ,
                            Rcpp::Named("B") = B);
}