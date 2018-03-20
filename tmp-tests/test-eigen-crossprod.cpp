// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::export]]
MatrixXd crossprodEigen(const Map<MatrixXd> x) {
  return x.transpose() * x;
}

// [[Rcpp::export]]
arma::mat crossprodArma(const arma::mat& x) {
  return x.t() * x;
}

// [[Rcpp::export]]
MatrixXd invCrossprodEigen(const Map<MatrixXd> x) {
  return crossprodEigen(x).inverse();
}

// // [[Rcpp::export]]
// MatrixXd invCrossprodEigen2(const Map<MatrixXd> x) {
//   MatrixXd I = MatrixXd::Identity(x.cols(), x.cols());
//   return crossprodEigen(x).llt().solve(I);
// }

// [[Rcpp::export]]
arma::mat invCrossprodArma(const arma::mat& x) {
  return inv_sympd(crossprodArma(x));
}

// [[Rcpp::export]]
arma::mat invCrossprodArmaRow(const arma::mat& x,
                           const arma::Row<uint32_t>& ind) {
  return inv_sympd(crossprodArma(x.rows(ind)));
}


/*** R
mat <- matrix(0, 1000, 20); mat[] <- rnorm(length(mat))

all.equal(K <- crossprodEigen(mat), K2 <- crossprod(mat))
all.equal(crossprodArma(mat), crossprod(mat))
microbenchmark::microbenchmark(
  crossprodEigen(mat),
  crossprodArma(mat),
  crossprod(mat)
)

all.equal(invCrossprodEigen(mat), solve(K2))
all.equal(invCrossprodArma(mat), solve(K2))
microbenchmark::microbenchmark(
  invCrossprodEigen(mat),
  invCrossprodArma(mat),
  invCrossprodArmaRow(mat, 0:999),
  invCrossprodArmaRow(mat, 0:50),
  solve(crossprod(mat))
)
*/
