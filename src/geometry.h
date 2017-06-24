#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

arma::mat cart2bary_cpp(const arma::mat X, const arma::mat P);