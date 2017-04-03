#include <RcppArmadillo.h>
#include "introUtils.h"
#include "procrustes.h"
#include "registration.h"
#include "toolbox.h"
#include "wilcoxon.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

//' Axis of projection
//' 
//' \code{get_axis} returns the axis onto which projection should be performed.
//' 
//' @param uglob a matrix of global scores.
//' @param lab a vector of integers.
//' @param anc1 an integer.
//' @param anc2 an integer.
//' 
//' @return The returned value is a numeric vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::vec get_axis(arma::mat &uglob, const arma::vec &lab, const int anc1, 
                   const int anc2){
  Rcpp::List res = cmpt_centroids(uglob, lab, anc1, anc2);
  arma::vec m1 = res[0];
  arma::vec m2 = res[1];
  return(m2 - m1);
}

//' Directional statistics
//' 
//' \code{cmpt_directional_stat} computes the displacement of admixed 
//' individuals such that positive (resp. negative) displacement corresponds to 
//' a displacement towards ancestral population 2 (resp. 1).
//' 
//' @param usc a matrix of rescaled local scores.
//' @param uglob a matrix of global scores.
//' @param lab a vector of integers.
//' @param adm an integer.
//' @param ax a numeric vector.
//' 
//' @return The returned value is a numeric value.
//' 
//' @export
//' 
// [[Rcpp::export]]
double cmpt_directional_stat(arma::mat &usc,
                             arma::mat &uglob, 
                             const arma::vec &lab, 
                             const int adm, 
                             arma::vec &ax){
  int nIND = uglob.n_rows; 
  double stat = 0;
  for (int j = 0; j < nIND; j++){
    if (lab[j] == adm){
      stat += arma::dot(usc.row(j) - uglob.row(j), ax);
    }
  }
  return(stat);
}

//' Wilcoxon statistics
//' 
//' \code{cmpt_all_stat} computes the statistics.
//' 
//' @param geno a genotype matrix.
//' @param V a loading matrix.
//' @param sigma a vector of singular values.
//' @param window_size an integer.
//' @param direction an integer.
//' @param lab a vector of integers.
//' @param ancstrl1 an integer.
//' @param ancstrl2 an integer.
//' @param adm an integer.
//' @param axis a numeric vector.
//' 
//' @return The returned value is a numeric vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::vec cmpt_all_stat(const arma::mat &geno, 
                        const arma::mat &V, 
                        const arma::vec &sigma, 
                        const int window_size,  
                        const int direction, 
                        const arma::vec lab, 
                        const int ancstrl1,
                        const int ancstrl2,
                        const int adm, 
                        const arma::vec axis){
  int nSNP = geno.n_cols;
  int nIND = geno.n_rows;
  int K = V.n_cols;
  arma::vec s(K, arma::fill::zeros);
  arma::vec dglob(K, arma::fill::zeros);
  arma::vec dloc(K, arma::fill::zeros);
  arma::mat usc(nIND, K, arma::fill::zeros);
  arma::vec stat(nSNP, arma::fill::zeros);
  arma::vec ax(K, arma::fill::zeros);
  arma::mat R(K, K, arma::fill::zeros); // ROTATION CORRECTION
  R.eye();
  
  arma::mat uglob = cmpt_global_pca(geno, V, sigma);
  arma::mat uloc = cmpt_local_pca(geno, V, sigma, 0, window_size);
  ax = get_axis(uglob, lab, ancstrl1, ancstrl2);
  for (int k = 0; k < K; k++){
    ax[k] *= axis[k];   
  }
  cmpt_transformation(uloc, uglob, lab, ancstrl1, ancstrl2, s, dloc, dglob, R);
  usc = rescale_local_pca(uloc, s, dglob, dloc, R);
  for (int i = 1; i < (nSNP - window_size); i++){
    updt_local_scores(uloc, geno, V, sigma, i, i + window_size);
    cmpt_transformation(uloc, uglob, lab, ancstrl1, ancstrl2, s, dloc, dglob, R);
    usc = rescale_local_pca(uloc, s, dloc, dglob, R);
    stat[i] = cmpt_directional_stat(usc, uglob, lab, adm, ax);
    //stat[i] = cmpt_wilcoxon_stat(usc, uglob, direction, lab, adm, axis);
  }
  stat[0] = stat[1];
  for (int i = (nSNP - window_size); i < nSNP; i++){
    stat[i] = stat[nSNP - window_size - 1];
  }
  return(stat);
}


//' \code{cmpt_new_win} computes the statistics.
//' 
//' @param i an integer.
//' @param map a vector containing the genetic positions in Morgans.
//' @param window_size a numeric value specifying the window size en Morgans.
//' 
//' @return The returned value is a numeric vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::vec cmpt_new_win(int i, const arma::vec &map, const double window_size){
  int n = map.n_elem;
  double half_window = window_size / 2.0;
  int idx_left = i;
  int idx_right = i;
  while ((map[i] - map[idx_left] < half_window) && (idx_left > 0)){
    idx_left -= 1;
  }
  while ((map[idx_right] - map[i] < half_window) && (idx_right < n - 1)){
    idx_right += 1;
  }
  arma::vec lr(2, arma::fill::zeros);
  lr[0] = idx_left;
  lr[1] = idx_right;
  return(lr);
}
