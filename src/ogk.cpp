#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector colMedian_cpp(arma::mat &x){
  int ncol = x.n_cols;
  NumericVector out(ncol);
  for (int j = 0; j < ncol; j++){ 
    out[j] = arma::median(x.col(j));
  }
  return out;
}

// [[Rcpp::export]]
double Erho_cpp(double b){
  double pb;
  double db;
  double res;
  pb = R::pnorm5(b, 0.0, 1.0, 1, 0);
  db = R::dnorm4(b, 0.0, 1.0, 0);
  res = 2 * ((1 - b * b) * pb - b * db + b * b) - 1;
  return(res);
}

// [[Rcpp::export]]
double Es2_cpp(double c){
  double q;
  double tmp;
  double res;
  q = R::qnorm5(0.75, 0.0, 1.0, 1, 0);
  tmp = c * q;
  res = Erho_cpp(tmp);
  return(res);
}

// [[Rcpp::export]]
Rcpp::List scaleTau2_matrix(arma::mat &x){
  int n = x.n_rows;
  int p = x.n_cols;
  double c1 = 4.5;
  double c2 = 3.0;
  double nEs2;
  NumericVector medx(p);
  NumericVector sigma0(p);
  medx = colMedian_cpp(x);
  NumericVector xvec(n);
  NumericVector rho(n);
  NumericVector w(n);
  NumericVector muvec(p);
  double mu;
  double w_tot;
  double sum_rho;
  double tmp;
  for (int j = 0; j < p; j++){
    mu = 0;
    w_tot = 0;
    sum_rho = 0;
    for (int i = 0; i < n; i++) {
      xvec[i] = std::abs(x(i, j) - medx[j]);
    }
    sigma0[j] = Rcpp::median(xvec);
    for (int i = 0; i < n; i++){
      xvec[i] /= (sigma0[j] * c1);
      w[i] = 1 - (xvec[i] * xvec[i]);
      w[i] = (std::abs(w[i]) + w[i]) * (std::abs(w[i]) + w[i]) / 4;
      mu += x(i, j) * w[i];
      w_tot += w[i]; 
    }
    muvec[j] = mu / w_tot;
    for (int i = 0; i < n; i++){
      tmp = ((x(i, j) - muvec[j]) / sigma0[j]) * ((x(i, j) - muvec[j]) / sigma0[j]);
      if (tmp > (c2 * c2)){
        rho[i] = c2 * c2;  
      } else {
        rho[i] = tmp;
      }
      sum_rho += rho[i];
    }
    nEs2 = n * Es2_cpp(c2);
    sigma0[j] *= sqrt(sum_rho / nEs2);
  }
  return Rcpp::List::create(Rcpp::Named("mu") = muvec,
                            Rcpp::Named("sigma0") = sigma0);
}

// [[Rcpp::export]]
NumericVector scaleTau2_vector(arma::vec &x){
  int n = x.size();
  double c1 = 4.5;
  double c2 = 3.0;
  double nEs2;
  double medx;
  double sigma0;
  medx = arma::median(x);
  NumericVector xvec(n);
  NumericVector rho(n);
  NumericVector w(n);
  double mu;
  double w_tot;
  double sum_rho;
  double tmp;
  mu = 0;
  w_tot = 0;
  sum_rho = 0;
  for (int i = 0; i < n; i++) {
    xvec[i] = std::abs(x[i] - medx);
  }
  sigma0 = Rcpp::median(xvec);
  for (int i = 0; i < n; i++){
    xvec[i] /= (sigma0 * c1);
    w[i] = 1 - (xvec[i] * xvec[i]);
    w[i] = (std::abs(w[i]) + w[i]) * (std::abs(w[i]) + w[i]) / 4;
    mu += x[i] * w[i];
    w_tot += w[i]; 
  }
  mu /= w_tot;
  for (int i = 0; i < n; i++){
    tmp = ((x[i] - mu) / sigma0) * ((x[i] - mu) / sigma0);
    if (tmp > (c2 * c2)){
      rho[i] = c2 * c2;  
    } else {
      rho[i] = tmp;
    }
    sum_rho += rho[i];
  }
  nEs2 = n * Es2_cpp(c2);
  sigma0 *= sqrt(sum_rho / nEs2);
  NumericVector musigma(2);
  musigma[0] = mu;
  musigma[1] = sigma0;
  return(musigma);
}

// [[Rcpp::export]]
double covGK_cpp(arma::vec x, arma::vec y){
  arma::vec sum_xy = x + y;
  arma::vec diff_xy = x - y;
  NumericVector ms_sum;
  NumericVector ms_diff;
  ms_sum = scaleTau2_vector(sum_xy);
  ms_diff = scaleTau2_vector(diff_xy);
  return(((ms_sum[1] * ms_sum[1]) - (ms_diff[1] * ms_diff[1])) / 4);
}

// [[Rcpp::export]]
Rcpp::List ogk_step(arma::mat &x){
  int p = x.n_cols;
  Rcpp::List ms;
  ms = scaleTau2_matrix(x);
  NumericVector sigma0 = ms[1];
  arma::mat U(p, p);
  for (int i = 0; i < p; i++){
    for (int j = i; j < p; j++){
      if (i == j){
        U(i, j) = 1.0;  
      } else {
        U(i, j) = covGK_cpp(x.col(i) / sigma0[i], x.col(j) / sigma0[j]);
        U(j, i) = U(i, j);
      }
    }
  }
  arma::vec eigval;
  arma::mat eigvec;
  arma::mat V;
  arma::eig_sym(eigval, eigvec, U, "std");
  arma::mat D = arma::diagmat(as<arma::vec>(ms[1]));
  V = x * (arma::solve(D, eigvec));
  ms = scaleTau2_matrix(V);
  arma::mat L = arma::diagmat(as<arma::vec>(ms[1]) % as<arma::vec>(ms[1]));
  arma::mat A = D * eigvec;
  return Rcpp::List::create(Rcpp::Named("cov") = A * L * A.t(),
                            Rcpp::Named("center") = A * as<arma::vec>(ms[0]),
                            Rcpp::Named("V") = V,
                            Rcpp::Named("A") = A);
}

// [[Rcpp::export]]
arma::vec getDistance_cpp(arma::mat &x, arma::rowvec center, arma::mat cov){
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i = 0; i < n; i++){
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}

//' Robust estimates for location and scatter
//'
//' \code{covRob_cpp} implements the Orthogonalized Gnanadesikan-Kettenring estimator
//' of Maronna and Zamar.
//'
//' @param x a numeric matrix.
//' 
//' @examples
//' ## see also ?pcadapt for examples
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::List covRob_cpp(arma::mat& x){
  int n = x.n_rows;
  int p = x.n_cols;
  double df = (double) p;
  
  /* First iteration */
  Rcpp::List msva;   
  msva = ogk_step(x);
  arma::mat V = as<arma::mat>(msva[2]);
  arma::mat A = as<arma::mat>(msva[3]);
  
  /* Second iteration */
  msva = ogk_step(V);
  arma::mat cov = as<arma::mat>(msva[0]);
  arma::vec center = as<arma::mat>(msva[1]);
  cov = A * cov * A.t();
  center = A * center;
  arma::mat Z = as<arma::mat>(msva[2]);
  
  Rcpp::List musigma;
  musigma = scaleTau2_matrix(Z);
  NumericVector mu = musigma[0];
  NumericVector sigma0 = musigma[1];
  NumericVector d(n);
  double medd;
  for (int j = 0; j < p; j++){
    for (int i = 0; i < n; i++){
      Z(i, j) -= mu[j];
      Z(i, j) /= sigma0[j];
      d[i] += Z(i, j) * Z(i, j); 
    }
  }
  medd = Rcpp::median(d);
  double tmp = R::qchisq(0.5, df, 1, 0);
  double cdelta = medd / tmp;
  double beta = 0.9;
  double quantile = R::qchisq(beta, df, 1, 0);
  double qq;
  qq = quantile * cdelta;
  IntegerVector wt(n);
  double sum_wt = 0;
  arma::rowvec wcenter(p);
  wcenter.zeros();
  arma::mat wcov(p, p);
  wcov.zeros();
  for (int i = 0; i < n; i++){
    if (d[i] < qq){
      wt[i] = 1;
      sum_wt += 1;
      for (int j = 0; j < p; j++){
        wcenter[j] += x.at(i, j);
      }
    } 
  }
  for (int j = 0; j < p; j++){
    wcenter[j] /= sum_wt;
  }
  for (int i = 0; i < p; i ++){
    for (int j = i; j < p; j++){
      for (int k = 0; k < n; k++){
        if (wt[k] == 1){
          wcov.at(i, j) += (x.at(k, i) - wcenter[i]) * (x.at(k, j) - wcenter[j]) / sum_wt;
        } 
      }
      wcov.at(j, i) = wcov.at(i, j); 
    }
  }
  arma::vec wdist;
  wdist = getDistance_cpp(x, wcenter, wcov);
  return Rcpp::List::create(Rcpp::Named("cov") = wcov,
                            Rcpp::Named("center") = wcenter,
                            Rcpp::Named("dist") = wdist);
}
