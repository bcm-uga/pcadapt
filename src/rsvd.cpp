#include <RcppArmadillo.h>
#include "toolbox.h"

// [[Rcpp::depends("RcppArmadillo")]]

#define NA 9

using namespace Rcpp;


//' Randomized SVD
//' 
//' \code{rsvd_cpp} computes the randomized SVD.
//' 
//' @param A a numeric matrix.
//' @param k an integer.
//' 
//' @return The returned value is a Rcpp::List.
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List rsvd_cpp(arma::mat &A, int k){
  
  // int m = A.n_rows; 
  int n = A.n_cols;
  
  int p = 10;
  int q = 1;
  
  // Oversampling parameter
  int l = k + p;
  if (l > n){
    l = n;
  }
  // int nu = k;
  // int nv = k;
  
  arma::mat O = arma::randu(n, l);
  arma::mat Y = A * O;
  arma::mat Q;
  arma::mat R;
  arma::mat Z;
  
  if (q > 0){
    for (int i = 0; i < q; i++){
      arma::qr_econ(Q, R, Y);
      Y = Q;
      Z = A.t() * Y;
      arma::qr_econ(Q, R, Z);
      Z = Q;
      Y = A * Z;
    }
  }
  arma::qr_econ(Q, R, Y);
  arma::mat B = Q.t() * A;
  return Rcpp::List::create(Rcpp::Named("Q") = Q ,
                            Rcpp::Named("B") = B);
}

//' A * x
//' 
//' \code{AX} computes the product A * x.
//' 
//' @param filename a character string.
//' @param x a numeric vector.
//' @param nSNP an integer.
//' @param nIND an integer.
//' 
//' @return The returned value is a numeric vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::vec AX(std::string filename, arma::vec x, int nSNP, int nIND){
  FILE *xfile;
  xfile = fopen(filename.c_str(), "r");
  arma::vec resprod(nSNP, arma::fill::zeros);
  float value;
  for (int i = 0; i < nSNP; i++){
    for (int j = 0; j < nIND; j++){
      if (fscanf(xfile, "%g", &value) != EOF){
        resprod[i] += (double) value * x[j];
      }
    }
  }
  fclose(xfile);
  return(resprod);
}

//' A.t * x
//' 
//' \code{tAX} computes the product A.t * x.
//' 
//' @param filename a character string.
//' @param x a numeric vector.
//' @param nSNP an integer.
//' @param nIND an integer.
//' 
//' @return The returned value is a numeric vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::vec tAX(std::string filename, arma::vec x, int nSNP, int nIND){
  FILE *xfile;
  xfile = fopen(filename.c_str(), "r");
  arma::vec resprod(nIND, arma::fill::zeros);
  float value;
  for (int i = 0; i < nSNP; i++){
    for (int j = 0; j < nIND; j++){
      if (fscanf(xfile, "%g", &value) != EOF){
        resprod[j] += (double) value * x[i];
      }
    }
  }
  fclose(xfile);
  return(resprod);
}
