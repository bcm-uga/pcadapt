#include <RcppArmadillo.h>
#include "procrustes.h"
#include "registration.h"
#include "toolbox.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

arma::mat cmpt_global_pca(const arma::mat &geno, 
                          const arma::mat &V, 
                          const arma::vec &sigma);

arma::mat cmpt_local_pca(const arma::mat &geno, const arma::mat &V, 
                         const arma::vec &sigma, const int beg, const int end);

void updt_local_scores_deprecated(arma::mat &u, const arma::mat &geno, 
                                  const arma::mat &V, const arma::vec &sigma, 
                                  const int beg, const int end);

void updt_local_scores(arma::mat &u, 
                         const arma::mat &geno, 
                         const arma::mat &V, 
                         const arma::vec &sigma, 
                         const int beg_old, 
                         const int end_old,
                         const int beg_new,
                         const int end_new);