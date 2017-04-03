#include <RcppArmadillo.h>
#include "procrustes.h"
#include "registration.h"
#include "toolbox.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

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
void updt_local_scores(arma::mat &u, const arma::mat &geno, const arma::mat &V, 
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
