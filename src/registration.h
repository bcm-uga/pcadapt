#include <RcppArmadillo.h>
#include "toolbox.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

Rcpp::List cmpt_centroids(arma::mat u, const arma::vec lab, const int pop1, 
                          const int pop2);

void cmpt_transformation(arma::mat &uloc, 
                         arma::mat &uglob, 
                         const arma::vec &lab, 
                         const int ancstrl1, 
                         const int ancstrl2, 
                         arma::vec &s, 
                         arma::vec &dloc, 
                         arma::vec &dglob, 
                         arma::mat &R);

arma::mat rescale_local_pca(arma::mat &u, arma::vec &s, arma::vec &dep_loc, 
                            arma::vec &dep_glob, arma::mat &R);