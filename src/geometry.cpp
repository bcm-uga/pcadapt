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
  int M = P.n_rows; // number of points
  int N = P.n_cols; // number of coordinates
  if (X.n_cols != N){
    // Centroids must have the same number of coordinates
    stop("Simplex X must have same number of columns as point matrix P");
  }
  if (X.n_rows != (N + 1)){
    stop("Simplex X must have N columns and N+1 rows");
  }
  arma::colvec J_N(N, arma::fill::ones);
  arma::colvec J_M(M, arma::fill::ones);
  arma::mat X1 = X.rows(0, N - 1) - (J_N * X.row(N));
  arma::mat A = (P - J_M * X.row(N)) ;
  arma::mat Beta(M, N + 1, arma::fill::zeros);
  arma::mat tBeta  = arma::solve(X1.t(), A.t()); 
  Beta.cols(0, N - 1) = tBeta.t();
  double value = 0.0;
  for (int i = 0; i < M; i++){
    value = sum(Beta.row(i));
    Beta(i, N) = 1.0 - value; 
  }
  return(Beta);
}