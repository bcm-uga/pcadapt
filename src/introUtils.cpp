#include <RcppArmadillo.h>
#include "procrustes.h"
#include "registration.h"
#include "toolbox.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

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

//' Global Principal Component Analysis
//' 
//' \code{cmpt_global_pca} computes the scores using all genetic markers.
//' 
//' @param geno a genotype matrix.
//' @param V a loading matrix.
//' @param sigma a vector of singular values.
//' 
//' @return The returned value is a matrix of scores.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::mat cmpt_global_pca(const arma::mat &geno, const arma::mat &V,
                          const arma::vec &sigma){
  int nIND = geno.n_rows;
  int nSNP = V.n_rows;
  int K = V.n_cols;
  arma::mat u(nIND, K, arma::fill::zeros);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      for (int i = 0; i < nSNP; i++){
        u(j, k) += geno(j, i) * V(i, k) / sigma[k];
      }
    }
  }
  return(u);
}

//' Local Principal Component Analysis
//' 
//' \code{cmpt_local_pca} computes the scores using a subset of genetic markers.
//' 
//' @param geno a genotype matrix.
//' @param V a loading matrix.
//' @param sigma a vector of singular values.
//' @param beg an integer specifying the first marker to be included.
//' @param end an integer specifying the first marker to be excluded.
//' 
//' @return The returned value is a matrix of scores.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::mat cmpt_local_pca(const arma::mat &geno, const arma::mat &V, 
                         const arma::vec &sigma, const int beg, const int end){
  // [beg, end) 
  int nIND = geno.n_rows;
  int nSNP = V.n_rows;
  int K = V.n_cols;
  arma::mat u(nIND, K, arma::fill::zeros);
  double cst = (double) nSNP / (end - beg);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      for (int i = beg; i < end; i++){
        u(j, k) += geno(j, i) * V(i, k) * cst / sigma[k];
      }
    }
  }
  return(u);
}

//' Update local Principal Component Analysis
//' 
//' \code{updt_local_scores} computes the scores using a subset of genetic 
//' markers.
//' 
//' @param u a score matrix.
//' @param geno a genotype matrix.
//' @param V a loading matrix.
//' @param sigma a vector of singular values.
//' @param beg an integer specifying the first marker to be included.
//' @param end an integer specifying the first marker to be excluded.
//' 
//' @return The returned value is a score matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
void updt_local_scores_deprecated(arma::mat &u, const arma::mat &geno, const arma::mat &V, 
                       const arma::vec &sigma, const int beg, const int end){
  int nIND = geno.n_rows; 
  int nSNP = geno.n_cols;
  int K = u.n_cols;
  double cst = (double) nSNP / (end - beg);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      u(j, k) -= (geno.at(j, beg - 1) * V(beg - 1, k)) * cst / sigma[k];
      u(j, k) += (geno.at(j, end) * V(end, k)) * cst / sigma[k];
    }
  }
}

//' Update local Principal Component Analysis
//' 
//' \code{updt_local_scores} computes the scores using a subset of genetic 
//' markers.
//' 
//' @param u a score matrix.
//' @param geno a genotype matrix.
//' @param V a loading matrix.
//' @param sigma a vector of singular values.
//' @param beg_old an integer specifying the first marker to be included.
//' @param end_old an integer specifying the first marker to be excluded.
//' @param beg_new an integer specifying the first marker to be included.
//' @param end_new an integer specifying the first marker to be excluded.
//' 
//' @return The returned value is a score matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
void updt_local_scores(arma::mat &u, 
                       const arma::mat &geno, 
                       const arma::mat &V, 
                       const arma::vec &sigma, 
                       const int beg_old, 
                       const int end_old,
                       const int beg_new,
                       const int end_new){
  int nIND = geno.n_rows; 
  int nSNP = geno.n_cols;
  int K = u.n_cols;
  double cst_old = (double) nSNP / (end_old - beg_old);
  u /= cst_old;
  
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      if (beg_new > beg_old){
        for (int p = beg_old; p < beg_new; p++){
          u(j, k) -= (geno.at(j, p) * V(p, k)) / sigma[k];
        }
      }
      if (end_new > end_old){
        for (int q = end_old; q < end_new; q++){
          u(j, k) += (geno.at(j, q) * V(q, k)) / sigma[k];
        }
      }
    }
  }
  
  double cst_new = (double) nSNP / (end_new - beg_new);
  u *= cst_new;
}

