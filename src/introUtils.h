#include <RcppArmadillo.h>
#include "toolbox.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

IntegerVector get_window(int i, 
                         const arma::vec &map, 
                         const double window_size,
                         const int side);

arma::mat cmpt_global_pca(const arma::mat &geno, 
                          const arma::mat &V, 
                          const arma::vec &sigma);

arma::mat cmpt_local_pca(const arma::mat &geno, const arma::mat &V, 
                         const arma::vec &sigma, const int beg, const int end);

void updt_local_scores(arma::mat &u, 
                       const arma::mat &geno, 
                       const arma::mat &V, 
                       const arma::vec &sigma, 
                       const int beg_old, 
                       const int end_old,
                       const int beg_new,
                       const int end_new);