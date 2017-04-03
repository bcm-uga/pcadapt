#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

arma::mat fJ_cpp(int n);

arma::mat fcnt_cpp(arma::mat &a);

arma::mat pca_rotation(arma::mat &a, arma::mat &b);