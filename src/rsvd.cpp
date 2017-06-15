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

//' @export
//' 
// [[Rcpp::export]]
arma::mat AX_fun(std::string filename, int l, int nSNP, int nIND){
  FILE *xfile;
  xfile = fopen(filename.c_str(), "r");
  Rprintf("Reading file %s...\n", filename.c_str());
  arma::mat resprod(l, nIND, arma::fill::zeros);
  
  for (int i = 0; i < nSNP; i++){
    arma::vec geno(nIND, arma::fill::zeros);
    arma::vec random_vec(l);
    random_vec = runif(l);
    float value;
    for (int j = 0; j < nIND; j++){
      if (fscanf(xfile, "%g", &value) != EOF){
        geno[j] = (double) value;
      }
      for (int p = 0; p < l; p++){
        resprod(p, j) += geno[j] * random_vec[p];
      }
    }
  }
  return(resprod);
}

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
