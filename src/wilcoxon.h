#include <RcppArmadillo.h>
#include "toolbox.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

arma::vec get_rank(const arma::vec &v_temp);

double cmpt_wilcoxon_stat(arma::mat &usc,
                          arma::mat &uglob, 
                          int direction, 
                          arma::vec &lab, 
                          int adm, 
                          int axis);