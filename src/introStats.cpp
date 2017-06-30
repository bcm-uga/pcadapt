#include <RcppArmadillo.h>
#include "introUtils.h"
#include "procrustes.h"
#include "registration.h"
#include "toolbox.h"
#include "wilcoxon.h"
#include "geometry.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

//' Axis of projection
//' 
//' \code{get_axis} returns the axis onto which projection should be performed.
//' 
//' @param u a matrix of global scores.
//' @param labels a vector of integers.
//' @param pop1 an integer.
//' @param pop2 an integer.
//' 
//' @return The returned value is a numeric vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::vec get_axis(arma::mat &u, 
                   const arma::vec &labels, 
                   const int pop1, 
                   const int pop2){
  Rcpp::List res = cmpt_centroids(u, labels, pop1, pop2);
  arma::vec m1 = res[0];
  arma::vec m2 = res[1];
  return(m2 - m1);
}

//' Directional statistics
//' 
//' \code{cmpt_directional_stat} computes the displacement of admixed 
//' individuals such that a positive (resp. negative) displacement corresponds 
//' to a displacement towards ancestral population 2 (resp. 1).
//' 
//' @param usc a matrix of rescaled local scores.
//' @param uglob a matrix of global scores.
//' @param labels a vector of integers.
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
                             const arma::vec &labels, 
                             const int adm, 
                             arma::vec &ax){
  int nIND = uglob.n_rows; 
  double stat = 0;
  for (int j = 0; j < nIND; j++){
    if (labels[j] == adm){
      stat += arma::dot(usc.row(j) - uglob.row(j), ax);
    }
  }
  return(stat);
}


//' \code{get_window}
//' 
//' @param i an integer.
//' @param map a vector containing the genetic positions in Morgans.
//' @param window_size a numeric value specifying the window size en Morgans.
//' @param side an integer specifying whether the window should be aligned on 
//' the left, middle or right.
//' 
//' @return The returned value is a numeric vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
IntegerVector get_window(int i, 
                         const arma::vec &map, 
                         const double window_size,
                         const int side){
  int n = map.n_elem;
  double half_window = window_size / 2.0;
  int idx_left = i;
  int idx_right = i;
  
  if (side == -1){
    while ((map[i] - map[idx_left] < 2 * half_window) && (idx_left > 0)){
      idx_left -= 1;
    }
  } else if (side == 0){
    while ((map[i] - map[idx_left] < half_window) && (idx_left > 0)){
      idx_left -= 1;
    }
    while ((map[idx_right] - map[i] < half_window) && (idx_right < n - 1)){
      idx_right += 1;
    } 
  } else if (side == 1){
    while ((map[idx_right] - map[i] < 2 * half_window) && (idx_right < n - 1)){
      idx_right += 1;
    } 
  }
  
  IntegerVector lr(2);
  lr[0] = idx_left;
  lr[1] = idx_right;
  
  return(lr);
}


//' Introgression statistics
//' 
//' \code{cmpt_stat_introgr} computes the statistics.
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
//' @param map a numeric vector containing the genetic positions.
//' @param with_map an integer specifying whether the genetic positions have
//' been provided.
//' @param side an integer specifying whether the window should be aligned on 
//' the left, middle or right.
//' 
//' @return The returned value is a numeric vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::vec cmpt_stat_introgr(const arma::mat &geno, 
                            const arma::mat &V, 
                            const arma::vec &sigma, 
                            const int window_size,  
                            const int direction, 
                            const arma::vec lab, 
                            const int ancstrl1,
                            const int ancstrl2,
                            const int adm, 
                            const arma::vec axis,
                            const arma::vec map,
                            const int with_map,
                            const int side){
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
  
  IntegerVector idx_old(2);
  IntegerVector idx_new(2);
  
  int hws = window_size / 2;
  int loop_beg;
  int loop_end;
  
  if (side == -1 and with_map == 0){
    idx_old = get_window(window_size, map, window_size, side);    
    loop_beg = window_size + 1;
    loop_end = nSNP;
  } else if (side == 0 and with_map == 0){
    idx_old = get_window(hws, map, window_size, side);    
    loop_beg = hws + 1;
    loop_end = nSNP - hws;
  } else if (side == 1 and with_map == 0){
    idx_old = get_window(0, map, window_size, side);    
    loop_beg = 1;
    loop_end = window_size;
  } else if (with_map == 1){
    idx_old = get_window(0, map, window_size, side);
    loop_beg = idx_old[1] / 2;
    IntegerVector idx_last;
    idx_last = get_window(nSNP - 1, map, window_size, side);
    loop_end = idx_last[0] + (idx_last[1] - idx_last[0]) / 2; 
  }
  
  arma::mat uloc = cmpt_local_pca(geno, V, sigma, idx_old[0], idx_old[1]);
  
  ax = get_axis(uglob, lab, adm, ancstrl2); // from adm to ancstrl2
  for (int k = 0; k < K; k++){
    ax[k] *= axis[k];   
  }
  cmpt_transformation(uloc, uglob, lab, ancstrl1, ancstrl2, s, dloc, dglob, R);
  usc = rescale_local_pca(uloc, s, dglob, dloc, R);
  for (int i = loop_beg; i < loop_end; i++){
    idx_new = get_window(i, map, window_size, side);
    updt_local_scores(uloc, geno, V, sigma, idx_old[0], idx_old[1],
                      idx_new[0], idx_new[1]);
    cmpt_transformation(uloc, uglob, lab, ancstrl1, ancstrl2, s, dloc, dglob, R);
    usc = rescale_local_pca(uloc, s, dloc, dglob, R);
    stat[i] = cmpt_directional_stat(usc, uglob, lab, adm, ax);
    idx_old[0] = idx_new[0];
    idx_old[1] = idx_new[1];
  }
  for (int i = 0; i < loop_beg; i++){
    stat[i] = NA_REAL;
  }
  for (int i = loop_end; i < nSNP; i++){
    stat[i] = NA_REAL;
  }
  stat[0] = 0;
  stat[nSNP - 1] = 0;
  return(stat);
}


//' Introgression statistics
//' 
//' \code{cmpt_stat_introgr_bary} computes the statistics.
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
//' @param map a numeric vector containing the genetic positions.
//' @param with_map an integer specifying whether the genetic positions have
//' been provided.
//' @param side an integer specifying whether the window should be aligned on 
//' the left, middle or right.
//' 
//' @return The returned value is a numeric vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::vec cmpt_stat_introgr_bary(const arma::mat &geno, 
                            const arma::mat &V, 
                            const arma::vec &sigma, 
                            const int window_size,  
                            const int direction, 
                            const arma::vec lab, 
                            const int ancstrl1,
                            const int ancstrl2,
                            const int adm, 
                            const arma::vec axis,
                            const arma::vec map,
                            const int with_map,
                            const int side){
  int nSNP = geno.n_cols;
  int nIND = geno.n_rows;
  int K = V.n_cols;
  arma::mat usc(nIND, K, arma::fill::zeros);
  arma::vec stat(nSNP, arma::fill::zeros);  
  arma::mat uglob = cmpt_global_pca(geno, V, sigma);
  
  IntegerVector idx_old(2);
  IntegerVector idx_new(2);
  
  int hws = window_size / 2;
  int loop_beg;
  int loop_end;
  
  if (side == -1 and with_map == 0){
    idx_old = get_window(window_size, map, window_size, side);    
    loop_beg = window_size + 1;
    loop_end = nSNP;
  } else if (side == 0 and with_map == 0){
    idx_old = get_window(hws, map, window_size, side);    
    loop_beg = hws + 1;
    loop_end = nSNP - hws;
  } else if (side == 1 and with_map == 0){
    idx_old = get_window(0, map, window_size, side);    
    loop_beg = 1;
    loop_end = window_size;
  } else if (with_map == 1){
    idx_old = get_window(0, map, window_size, side);
    loop_beg = idx_old[1] / 2;
    IntegerVector idx_last;
    idx_last = get_window(nSNP - 1, map, window_size, side);
    loop_end = idx_last[0] + (idx_last[1] - idx_last[0]) / 2; 
  }
  
  arma::mat uloc = cmpt_local_pca(geno, V, sigma, idx_old[0], idx_old[1]);
  arma::mat mglob = cmpt_centroids_bary(uglob, lab, ancstrl1, ancstrl2);
  arma::mat mloc(nIND, 1, arma::fill::zeros);
  mloc.col(0) = uglob.col(0);
  arma::mat res_glob = cart2bary_cpp(mglob, mloc);
  int number_of_admixed = 0;
  number_of_admixed = get_nb_ind(lab, adm);
  arma::mat res;
  for (int i = loop_beg; i < loop_end; i++){
    idx_new = get_window(i, map, window_size, side);
    updt_local_scores(uloc, geno, V, sigma, idx_old[0], idx_old[1], idx_new[0], idx_new[1]);
    arma::mat mglob = cmpt_centroids_bary(uloc, lab, ancstrl1, ancstrl2);
    mloc.col(0) = uloc.col(0);
    res = cart2bary_cpp(mglob, mloc);
    for (int j = 0; j < nIND; j++){
      if (lab[j] == adm){
        //stat[i] += (res(j, 1) - res_glob(j, 1)) / number_of_admixed;
        stat[i] += res(j, 1) / number_of_admixed;
      }
    }
    idx_old[0] = idx_new[0];
    idx_old[1] = idx_new[1];
  }
  double mean = sum(stat) / (loop_end - loop_beg);
  for (int i = 0; i < loop_beg; i++){
    stat[i] = NA_REAL;
  }
  for (int i = loop_end; i < nSNP; i++){
    stat[i] = NA_REAL;
  }
  stat[0] = mean;
  stat[nSNP - 1] = mean;
  return(stat);
}



