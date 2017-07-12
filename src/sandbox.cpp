#include <RcppArmadillo.h>
#include "introUtils.h"
#include "toolbox.h"
#include "geometry.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

//' Centroids
//' 
//' \code{scores_centroids_cpp} returns the centroids of each population.
//' 
//' @param scores a numeric matrix.
//' @param pop a string vector.
//' @param popUnique a string vector.
//' 
//' @return The returned value is a numeric matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::mat scores_centroids_cpp(const arma::mat &scores,
                               const StringVector &pop,
                               const StringVector &popUnique){
  int nb_pop = popUnique.size();
  int nb_ind = scores.n_rows;
  int K = scores.n_cols;
  
  arma::mat centroids(nb_pop, K, arma::fill::zeros);
  for (int i = 0; i < nb_pop; i++){
    int tmp = 0; // counts the number of individuals in the i-th population
    for (int j = 0; j < nb_ind; j++){
      if (pop[j] == popUnique[i]){
        tmp++;
        for (int k = 0; k < K; k++){
          centroids(i, k) += scores(j, k);    
        }
      }
    }
    for (int k = 0; k < K; k++){
      centroids(i, k) /= (double) tmp;
    }
  }
  return(centroids);
}

//' Simplex
//' 
//' \code{centroids_to_simplex_cpp}
//' 
//' @param centroids a numeric matrix.
//' @param popUnique a string vector.
//' @param admixed a character vector.
//' 
//' @return The returned value is a numeric matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::mat centroids_to_simplex_cpp(const arma::mat &centroids,
                                   const StringVector &popUnique,
                                   const CharacterVector &admixed){
  
  int nb_pop = popUnique.size();
  int K = centroids.n_cols;
  std::string str_admixed = Rcpp::as<std::string> (admixed);
  
  arma::mat simplex(nb_pop - 1, K);
  int idx = 0;
  for (int i = 0; i < nb_pop; i++){
    if (popUnique[i] != str_admixed){
      for (int k = 0; k < K; k++){
        simplex(idx, k) = centroids(i, k);
      }
      idx++;
    }  
  }
  return(simplex);
}


//' Barycentric coordinates to ancestry
//' 
//' \code{centroids_to_simplex_cpp}
//' 
//' @param scores a numeric matrix.
//' @param pop a string vector.
//' @param popUnique a string vector.
//' @param admixed a character vector.
//' 
//' @return The returned value is a numeric vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::vec bary_to_ancestry_cpp(const arma::mat &scores,
                               const StringVector &pop,
                               const StringVector &popUnique,
                               const CharacterVector &admixed){
  arma::mat centroids = scores_centroids_cpp(scores,
                                             pop,
                                             popUnique);
  arma::mat simplex = centroids_to_simplex_cpp(centroids,
                                               popUnique,
                                               admixed);
  int nb_admixed = get_nb_ind(pop, admixed);
  int K = scores.n_cols;
  std::string str_admixed = Rcpp::as<std::string> (admixed);

  arma::mat scores_admixed(nb_admixed, K, arma::fill::zeros);
  int idx = 0;
  for (int i = 0; i < scores.n_rows; i++){
    if (pop[i] == str_admixed){
      for (int k = 0; k < K; k++){
        scores_admixed(idx, k) = scores(i, k);
      }
      idx++;
    }
  }
  
  arma::mat bary_coord = cart2bary_cpp(simplex, scores_admixed);
  
  arma::vec ancestry(bary_coord.n_cols, arma::fill::zeros);
  for (int k = 0; k < bary_coord.n_cols; k++){
    ancestry[k] = mean(bary_coord.col(k));
  }
  
  return(ancestry);
}


//' Simplex
//' 
//' \code{centroids_to_simplex_cpp}
//' 
//' @param sgeno a scaled genotype matrix.
//' @param d a numeric vector.
//' @param v a numeric matrix.
//' @param pop a string vector.
//' @param popUnique a string vector.
//' @param admixed a character vector.
//' @param window_size a numeric value.
//' @param map a numeric vector.
//' 
//' @return The returned value is a numeric matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::mat slidingWindows(const arma::mat &sgeno,
                         const arma::vec &d,
                         const arma::mat &v,
                         const StringVector &pop,
                         const StringVector &popUnique,
                         const CharacterVector &admixed,
                         const double window_size,
                         const NumericVector &map){
  
  int nSNP = sgeno.n_cols;
  int nIND = sgeno.n_rows;
  int nPOP = popUnique.size();
  int hws = (int) window_size / 2;
  int K = std::max(1, nPOP - 2);
  
  IntegerVector ix_o = get_window(hws, map, window_size);    
  IntegerVector ix_n(2);
  
  int loop_start = hws + 1;
  int loop_end = nSNP - hws;
  arma::mat u = cmpt_local_pca(sgeno, v, d, ix_o[0], ix_o[1]);
  
  arma::mat stat(nSNP, nPOP - 1);
  stat.fill(NA_REAL);
  arma::mat tmp(nIND, K, arma::fill::zeros);
  
  for (int i = loop_start; i < loop_end; i++){
    ix_n = get_window(i, map, window_size);  
    updt_local_scores(u, sgeno, v, d, ix_o[0], ix_o[1], ix_n[0], ix_n[1]);
    for (int k = 0; k < K; k++){
      tmp.col(k) = u.col(k);
    }
    stat.row(i) = bary_to_ancestry_cpp(tmp, pop, popUnique, admixed).t();
    ix_o[0] = ix_n[0];
    ix_o[1] = ix_n[1];
  }
  
  return(stat);
  
}
