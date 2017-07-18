#include <RcppArmadillo.h>
#include "introUtils.h"
#include "toolbox.h"
#include "geometry.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

//' Get population size
//' 
//' \code{get_pop_size} 
//' 
//' @param pop a string vector.
//' @param popUnique a string vector.
//' 
//' @return The returned value is a numeric matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
IntegerVector get_pop_size(const StringVector &pop,
                           const StringVector &popUnique){
  
  int nb_ind = pop.size();
  int nb_pop = popUnique.size();
  IntegerVector popSize(nb_pop);
  
  for (int i = 0; i < nb_pop; i++){
    for (int j = 0; j < nb_ind; j++){
      if (pop[j] == popUnique[i]){
        popSize[i] += 1;  
      }
    }
  }
  return(popSize);
}

// [[Rcpp::export]]
void updt_centroids_cpp(arma::mat &centroids,
                        const arma::mat &scores,
                        const StringVector &pop,
                        const StringVector &popUnique,
                        const IntegerVector &popSize,
                        int K){
  int nb_pop = centroids.n_rows;
  int nb_ind = scores.n_rows;
  
  for (int i = 0; i < nb_pop; i++){
    for (int k = 0; k < K; k++){
      centroids(i, k) = 0.0;
      for (int j = 0; j < nb_ind; j++){
        if (pop[j] == popUnique[i]){
          centroids(i, k) += (double) scores(j, k) / popSize[i];    
        }
      }
    }
  }
}

// [[Rcpp::export]]
void updt_simplex_cpp(arma::mat &simplex,
                      const arma::mat &centroids,
                      const StringVector &popUnique,
                      const CharacterVector &admixed){
  
  std::string str_admixed = Rcpp::as<std::string> (admixed);
  
  int idx = 0;
  for (int i = 0; i < popUnique.size(); i++){
    if (popUnique[i] != str_admixed){
      for (int k = 0; k < centroids.n_cols; k++){
        simplex(idx, k) = centroids(i, k);
      }
      idx++;
    }  
  }
}

////////////////////////////////////////////////////////////////////////////////

//' Introgression statistics
//' 
//' \code{slidingWindows_fast} 
//' 
//' @param sgeno a scaled genotype matrix.
//' @param d a numeric vector.
//' @param v a numeric matrix.
//' @param pop a string vector.
//' @param popUnique a string vector.
//' @param admixed a character vector.
//' @param window_size a numeric value.
//' @param map a numeric vector.
//' @param with_map an integer.
//' 
//' @return The returned value is a numeric matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::mat slidingWindows_fast(const arma::mat &sgeno, 
                              const arma::vec &d, 
                              const arma::mat &v,
                              const StringVector &pop,
                              const StringVector &popUnique,
                              const CharacterVector &admixed,
                              const int window_size,  
                              const arma::vec map,
                              const int with_map){
  int nSNP = sgeno.n_cols;
  int nIND = sgeno.n_rows;
  int nPOP = popUnique.size();
  int K = std::max(1, nPOP - 2);
  std::string str_admixed = Rcpp::as<std::string> (admixed);
  
  int number_of_admixed = get_nb_ind(pop, admixed);
  IntegerVector popSize = get_pop_size(pop, popUnique);
  
  arma::mat stat(nSNP, nPOP - 1, arma::fill::zeros); // the simplex must have nPOP - 1 points
  stat.fill(NA_REAL);
  
  IntegerVector ix_o(2);
  IntegerVector ix_n(2);
  
  int hws = window_size / 2;
  int loop_start = 0;
  int loop_end = 0;
  
  if (with_map == 1){
    loop_start = 0;
    loop_end = nSNP - 1;
  } else {
    loop_start = hws + 1;
    loop_end = nSNP - hws;  
  }
  
  arma::mat u = cmpt_local_pca(sgeno, v, d, ix_o[0], ix_o[1]);
  arma::mat uK(nIND, K, arma::fill::zeros);
  arma::mat tmp(popUnique.size(), K, arma::fill::zeros);
  arma::mat simplex(popUnique.size() - 1, K, arma::fill::zeros);
  double value = 0.0;
  for (int i = loop_start; i < loop_end; i++){
    ix_n = get_window(i, map, window_size);
    updt_local_scores(u, sgeno, v, d, ix_o[0], ix_o[1], ix_n[0], ix_n[1]);
    updt_centroids_cpp(tmp, u, pop, popUnique, popSize, K);
    updt_simplex_cpp(simplex, tmp, popUnique, admixed);
    for (int k = 0; k < K; k++){
      uK.col(k) = u.col(k);  
    }
    arma::mat res = cart2bary_cpp(simplex, uK);
    
    for (int k = 0; k < nPOP - 1; k++){
      value = 0.0;
      for (int j = 0; j < nIND; j++){
        if (pop[j] == str_admixed){
          value += res(j, k) / number_of_admixed;  
        }
      }
      stat(i, k) = value;  
    }
    
    ix_o[0] = ix_n[0];
    ix_o[1] = ix_n[1];
    
  }
  
  return(stat);
}