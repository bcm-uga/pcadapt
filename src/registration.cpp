#include <RcppArmadillo.h>
#include "toolbox.h"

// [[Rcpp::depends("RcppArmadillo")]]

#define NA 9

using namespace Rcpp;

//' Ancestral populations centroids
//' 
//' \code{cmpt_centroids} returns the average scores for each population.
//' 
//' @param u a matrix of scores.
//' @param lab a vector of integers.
//' @param pop1 an integer.
//' @param pop2 an integer.
//' 
//' @return The returned value is a list.
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List cmpt_centroids(arma::mat u, const arma::vec lab, const int pop1, 
                          const int pop2){
  int nIND = u.n_rows;
  int K = u.n_cols;
  arma::vec m1(K, arma::fill::zeros);
  arma::vec m2(K, arma::fill::zeros);
  int c1 = get_nb_ind(lab, pop1);
  int c2 = get_nb_ind(lab, pop2);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      if (lab[j] == pop1){
        m1[k] += u(j, k) / c1;
      } else if (lab[j] == pop2){
        m2[k] += u(j, k) / c2;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("m1") = m1,
                            Rcpp::Named("m2") = m2);
}

//' Match global centroids and local centroids
//' 
//' \code{cmpt_transformation} computes the local centroid, the global centroid 
//' and the scaling factor.
//' 
//' @param uloc a matrix of local scores.
//' @param uglob a matrix of global scores.
//' @param lab a vector of integers.
//' @param ancstrl1 an integer.
//' @param ancstrl2 an integer.
//' @param s a numeric vector.
//' @param dloc a numeric vector.
//' @param dglob a numeric vector.
//' @param R a numeric matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
void cmpt_transformation(arma::mat &uloc, 
                         arma::mat &uglob, 
                         const arma::vec &lab, 
                         const int ancstrl1, 
                         const int ancstrl2, 
                         arma::vec &s, 
                         arma::vec &dloc, 
                         arma::vec &dglob, 
                         arma::mat &R){
  Rcpp::List mglob = cmpt_centroids(uglob, lab, ancstrl1, ancstrl2);
  Rcpp::List mloc = cmpt_centroids(uloc, lab, ancstrl1, ancstrl2);
  arma::vec mglob1 = mglob[0];
  arma::vec mglob2 = mglob[1];
  arma::vec mloc1 = mloc[0];
  arma::vec mloc2 = mloc[1];
  dglob = (mglob1 + mglob2);
  dglob /= 2.0;
  dloc =  (mloc1 + mloc2);
  dloc /= 2.0;
  int K = s.n_elem;
  for (int k = 0; k < K; k++){
    s[k] = fabs(mglob1[k] - mglob2[k]) / fabs(mloc1[k] - mloc2[k]);
  }
}

//' Rescale local scores
//' 
//' \code{rescale_local_pca} returns the rescaled local scores.
//' 
//' @param u a matrix of scores.
//' @param s a numeric vector.
//' @param dep_loc a numeric vector.
//' @param dep_glob a numeric vector.
//' @param R a numeric matrix.
//' 
//' @return The returned value is a list.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::mat rescale_local_pca(arma::mat &u, arma::vec &s, arma::vec &dep_loc, 
                            arma::vec &dep_glob, arma::mat &R){
  int nIND = u.n_rows;
  int K = u.n_cols;
  arma::mat usc(nIND, K);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      usc(j, k) = u(j, k) * s[k];
      usc(j, k) = usc(j, k) + dep_glob[k] - dep_loc[k];
    }
  }
  usc *= R;
  return(usc);
}