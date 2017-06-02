#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

#define NA 9

using namespace Rcpp;

//' Fast computation of rolling products.
//' 
//' \code{test} 
//' 
//' @param x a genotype matrix.
//' 
//' @return The returned value is a numerical vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
double logit_fun(double x){
  return(1 / (1 + exp(-x)));
}

//' Fast computation of rolling products.
//' 
//' \code{roll_prod} 
//' 
//' @param x a genotype matrix.
//' @param regcoeff a numerical vector.
//' @param window_size an integer.
//' 
//' @return The returned value is a numerical vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::vec roll_prod(const arma::mat &x, const arma::vec &regcoeff, 
                    int window_size){
  int nSNP = x.n_cols;
  int nIND = x.n_rows;
  arma::vec u(nIND, arma::fill::zeros);
  arma::vec stat(nSNP, arma::fill::zeros);
  double res = 0;
  int loop_beg = window_size / 2;
  int loop_end = nSNP - window_size;
  
  for (int i = 0; i < window_size; i++){
    for (int j = 0; j < nIND; j++){
      u[j] += regcoeff[i] * x(j, i);  
    }
  }
  
  for (int j = 0; j < nIND; j++){
    res += (u[j] / nIND);
  }
  stat[loop_beg] = logit_fun(res);
  
  for (int i = 1; i < loop_end; i++){
    res = 0;
    for (int j = 0; j < nIND; j++){
      u[j] -= regcoeff[i - 1] * x(j, i - 1);
      u[j] += regcoeff[i + window_size] * x(j, i + window_size);
      res += u[j] / nIND;
    } 
    stat[i + loop_beg] =  logit_fun(res);
  }
  for (int i = 0; i < loop_beg; i++){
    stat[i] = NA_REAL;
  }
  for (int i = (loop_end + window_size / 2); i < nSNP; i++){
    stat[i] = NA_REAL;
  }
  return(stat);
}
