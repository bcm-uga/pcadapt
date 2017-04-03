#include <RcppArmadillo.h>
#include "toolbox.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

//' Wilcoxon rank
//' 
//' \code{get_rank} returns the ordering index.
//' 
//' @param v_temp a numeric vector.
//' 
//' @return The returned value is a vector of integers.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::vec get_rank(const arma::vec &v_temp){
  int n = v_temp.n_elem;
  arma::vec v_sort(n);
  for (int i = 0; i < n; i++){
    v_sort[i] = v_temp[i];
  }
  arma::uvec idx = sort_index(v_sort);
  arma::vec rank(n);
  rank.zeros();
  for (int i = 0; i < n; i++){
    int tmp = idx[i];
    rank[tmp] = i + 1;
  }
  return(rank);
}

//' Wilcoxon statistics
//' 
//' \code{cmpt_wilcoxon_stat} computes the Wilcoxon statistics.
//' 
//' @param usc a matrix of rescaled local scores.
//' @param uglob a matrix of global scores.
//' @param direction an integer.
//' @param lab a vector of integers.
//' @param adm an integer.
//' @param axis a numeric value.
//' 
//' @return The returned value is a numeric value.
//' 
//' @export
//' 
// [[Rcpp::export]]
double cmpt_wilcoxon_stat(arma::mat &usc,
                          arma::mat &uglob, 
                          int direction, 
                          arma::vec &lab, 
                          int adm, 
                          int axis){
  int nIND = uglob.n_rows;
  int nAND = get_nb_ind(lab, adm);
  arma::vec tmp(nIND);
  arma::vec Z(nIND);
  Z.zeros();
  double diff;
  double m = 0;
  for (int j = 0; j < nIND; j++){
    if (lab[j] == adm){
      m += uglob(j, axis) / nAND;
    }
  }
  
  for (int j = 0; j < nIND; j++){
    diff = usc(j, axis) - m;
    tmp[j] = fabs(diff);
    if ((direction == 1) && (diff > 0)){
      Z[j] = 1;
    } else if ((direction == (-1)) && (diff < 0)){
      Z[j] = 1;
    }
  }
  arma::vec tmp_sort = get_rank(tmp);
  double W = 0;
  for (int j = 0; j < nIND; j++){
    if (lab[j] == adm){
      W += (double) Z[j] * tmp_sort[j];
    }
  }
  return(W);
}