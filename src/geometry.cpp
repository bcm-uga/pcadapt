#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

//' Cartesian coordinates to barycentric coordinates
//' 
//' \code{cart2bary_cpp} returns the barycentric coordinates.
//' 
//' @param X a matrix.
//' @param P a matrix.
//' 
//' @return The returned value is a numeric matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::mat cart2bary_cpp(const arma::mat X, const arma::mat P){
  int M = P.n_rows;
  int N = P.n_cols;
  if (X.n_cols != N){
    stop("Simplex X must have same number of columns as point matrix P");
  }
  if (X.n_rows != (N + 1)){
    stop("Simplex X must have N columns and N+1 rows");
  }
  arma::colvec J(N, arma::fill::ones);
  arma::mat X1 = X.rows(0, N - 1) - (J * X.row(N));
  
  return(X1);
}