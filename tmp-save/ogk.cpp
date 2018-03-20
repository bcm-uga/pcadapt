#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;               	      // 'maps' rather than copies 
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;

using namespace Rcpp;
using namespace Eigen;

NumericVector colMedian_rcpp(NumericMatrix &x) {
  int ncol = x.ncol();
  NumericVector out(ncol);
  for (int j = 0; j < ncol; j++) { 
    NumericVector tmp(x(_, j));
    out[j] = median(tmp);
  }
  return out;
}


double Erho_rcpp(double b){
  double pb;
  double db;
  double res;
  pb = R::pnorm5(b, 0.0, 1.0, 1, 0);
  db = R::dnorm4(b, 0.0, 1.0, 0);
  res = 2 * ((1 - b * b) * pb - b * db + b * b) - 1;
  return(res);
}

double Es2_rcpp(double c){
  double q;
  double tmp;
  double res;
  q = R::qnorm5(0.75, 0.0, 1.0, 1, 0);
  tmp = c * q;
  res = Erho_rcpp(tmp);
  return(res);
}

Rcpp::List scaleTau2_matrix_rcpp(NumericMatrix &x){
  int n = x.nrow();
  int p = x.ncol();
  double c1 = 4.5;
  double c2 = 3.0;
  double nEs2;
  NumericVector medx(p);
  NumericVector sigma0(p);
  medx = colMedian_rcpp(x);
  NumericVector xvec(n);
  NumericVector rho(n);
  NumericVector w(n);
  NumericVector muvec(p);
  double mu;
  double w_tot;
  double sum_rho;
  double tmp;
  for (int j = 0; j < p; j++) {
    mu = 0;
    w_tot = 0;
    sum_rho = 0;
    for (int i = 0; i < n; i++) {
      xvec[i] = std::abs(x(i, j) - medx[j]);
    }
    sigma0[j] = Rcpp::median(xvec);
    for (int i = 0; i < n; i++) {
      xvec[i] /= (sigma0[j] * c1);
      w[i] = 1 - (xvec[i] * xvec[i]);
      w[i] = (std::abs(w[i]) + w[i]) * (std::abs(w[i]) + w[i]) / 4;
      mu += x(i, j) * w[i];
      w_tot += w[i]; 
    }
    muvec[j] = mu / w_tot;
    for (int i = 0; i < n; i++) {
      tmp = ((x(i, j) - muvec[j]) / sigma0[j]) * ((x(i, j) - muvec[j]) / sigma0[j]);
      if (tmp > (c2 * c2)){
        rho[i] = c2 * c2;  
      } else {
        rho[i] = tmp;
      }
      sum_rho += rho[i];
    }
    nEs2 = n * Es2_rcpp(c2);
    sigma0[j] *= sqrt(sum_rho / nEs2);
  }
  return Rcpp::List::create(Rcpp::Named("mu") = muvec,
                            Rcpp::Named("sigma0") = sigma0);
}

NumericVector scaleTau2_vector_rcpp(NumericVector &x) {
  int n = x.size();
  double c1 = 4.5;
  double c2 = 3.0;
  double nEs2;
  double medx;
  double sigma0;
  medx = median(x);
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
  sigma0 = median(xvec);
  for (int i = 0; i < n; i++) {
    xvec[i] /= (sigma0 * c1);
    w[i] = 1 - (xvec[i] * xvec[i]);
    w[i] = (std::abs(w[i]) + w[i]) * (std::abs(w[i]) + w[i]) / 4;
    mu += x[i] * w[i];
    w_tot += w[i]; 
  }
  mu /= w_tot;
  for (int i = 0; i < n; i++) {
    tmp = ((x[i] - mu) / sigma0) * ((x[i] - mu) / sigma0);
    if (tmp > (c2 * c2)) {
      rho[i] = c2 * c2;  
    } else {
      rho[i] = tmp;
    }
    sum_rho += rho[i];
  }
  nEs2 = n * Es2_rcpp(c2);
  sigma0 *= sqrt(sum_rho / nEs2);
  NumericVector musigma(2);
  musigma[0] = mu;
  musigma[1] = sigma0;
  return(musigma);
}

double covGK_rcpp(NumericVector &x, NumericVector &y){
  NumericVector sum_xy = x + y;
  NumericVector diff_xy = x - y;
  NumericVector ms_sum;
  NumericVector ms_diff;
  ms_sum = scaleTau2_vector_rcpp(sum_xy);
  ms_diff = scaleTau2_vector_rcpp(diff_xy);
  return(((ms_sum[1] * ms_sum[1]) - (ms_diff[1] * ms_diff[1])) / 4);
}

Eigen::MatrixXd getEigenValues(NumericMatrix M) {
  Eigen::Map<Eigen::MatrixXd> Mm(Rcpp::as< Eigen::Map<Eigen::MatrixXd> > (M));
  SelfAdjointEigenSolver<Eigen::MatrixXd> es(Mm);
  return es.eigenvectors();
}

Rcpp::List ogk_step_rcpp(NumericMatrix &x) {
  int p = x.ncol();
  Rcpp::List ms;
  ms = scaleTau2_matrix_rcpp(x);
  NumericVector sigma0 = ms[1];
  NumericMatrix U(p, p);
  for (int i = 0; i < p; i++) {
    for (int j = i; j < p; j++) {
      if (i == j) {
        U(i, j) = 1.0;  
      } else {
        NumericVector tmp_i = x.column(i) / sigma0[i];
        NumericVector tmp_j = x.column(j) / sigma0[j];
        U(i, j) = covGK_rcpp(tmp_i , tmp_j);
        U(j, i) = U(i, j);
      }
    }
  }
  NumericVector eigval;
  Eigen::MatrixXd eigvec;
  eigvec = getEigenValues(U);
  Eigen::Map<Eigen::MatrixXd> eigenX(Rcpp::as< Eigen::Map<Eigen::MatrixXd> > (x));
  Eigen::Map<Eigen::VectorXd> ms_vec(Rcpp::as< Eigen::Map<Eigen::VectorXd> > (sigma0));
  Eigen::MatrixXd A = ms_vec.asDiagonal() * eigvec;
  Eigen::MatrixXd V = eigenX * ms_vec.asDiagonal().inverse() * eigvec;
  
  return Rcpp::List::create(Rcpp::Named("V") = V,
                            Rcpp::Named("A") = A);
  
}

Eigen::VectorXd getDistance_rcpp(NumericMatrix &x,
                                 NumericVector &center,
                                 NumericMatrix &cov) {
  
  int n = x.nrow();
  int p = x.ncol();
  Eigen::MatrixXd x_cen(n, p);
  Eigen::Map<Eigen::MatrixXd> covEigen(Rcpp::as< Eigen::Map<Eigen::MatrixXd> > (cov));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      x_cen(i, j) = x(i, j) - center[j];
    }
  }
  
  Eigen::VectorXd res = (x_cen * covEigen.inverse()).cwiseProduct(x_cen).rowwise().sum();
  
  return res;
  
}

// [[Rcpp::export]]
Rcpp::List covRob_rcpp(NumericMatrix &x){
  int n = x.nrow();
  int p = x.ncol();
  double df = (double) p;
  
  /* First iteration */
  Rcpp::List msva;   
  msva = ogk_step_rcpp(x);
  
  NumericMatrix V = msva[0];
  NumericMatrix A = msva[1];
  
  /* Second iteration */
  msva = ogk_step_rcpp(V);
  NumericMatrix Z = msva[0];
  
  Rcpp::List musigma;
  musigma = scaleTau2_matrix_rcpp(Z);
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
  
  NumericVector wcenter(p);
  NumericMatrix wcov(p, p);
  NumericVector wdist(n);
  
  for (int i = 0; i < n; i++) {
    if (d[i] < qq){
      wt[i] = 1;
      sum_wt += 1;
      for (int j = 0; j < p; j++) {
        wcenter[j] += x.at(i, j);
      }
    } 
  }
  for (int j = 0; j < p; j++) {
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
  
  wdist = getDistance_rcpp(x, wcenter, wcov);
  return Rcpp::List::create(Rcpp::Named("cov") = wcov,
                            Rcpp::Named("center") = wcenter,
                            Rcpp::Named("dist") = wdist);
}
