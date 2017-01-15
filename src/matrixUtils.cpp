#include <Rcpp.h>
using namespace Rcpp;

#define NA 9

//' Product of a matrix with its transpose
//' 
//' \code{tAA_cpp} computes the product of a real-valued matrix with its transpose. 
//' 
//' @param x a matrix.
//' @param nrow an integer specifying the number of rows of \code{x}.
//' @param ncol an integer specifying the number of columns of \code{x}.
//' 
//' @return The returned value is the result of the product of \code{x} with its transpose.
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericMatrix tAA_cpp(NumericMatrix x, int nrow, int ncol){
  NumericMatrix xcov(ncol,ncol);
  for (int i = 0; i < ncol; i++){
    for (int j = i; j < ncol; j++){
      for (int k = 0; k < nrow; k++){
        if ((x(k, i) != NA) && (x(k, j) != NA) && (!NumericVector::is_na(x(k, i))) && (!NumericVector::is_na(x(k, j)))){
          xcov(i, j) += x(k, i) * x(k, j);
        }
      }
      xcov(j, i) = xcov(i, j);
    }
  } 
  return(xcov);
}

int get_rows_matrix(double *xs, int nIND, int ploidy, double min_maf, double *maf_i){
  double mean = 0;
  double var = 0;
  double af = 0;
  int na = 0;
  
  for (int j = 0; j < nIND; j++){
    if ((xs[j] != NA) && (!NumericVector::is_na(xs[j]))){
      mean += xs[j]; 
    } else {
      na += 1;
    }
  }
  
  if (na >= nIND){
    Rcpp::stop("Detected SNP with missing values only, please remove it before proceeding.");  
  } else {
    mean /= (nIND - na);
    if (ploidy == 2){
      af = mean / 2.0;
      var = 2.0 * af * (1 - af);
    } else {
      af = mean;
      var = af * (1 - af);
    }
    if (af > 0.5){
      double tmp_af = af; 
      af = 1.0 - tmp_af;
    }
  }
  
  if (af >= min_maf){
    for (int j = 0; j < nIND; j++){
      if ((xs[j] != NA) && (!NumericVector::is_na(xs[j]))){
        if (var > 0){
          xs[j] -= mean;
          xs[j] /= sqrt(var);
        }
      } else {
        xs[j] = 0;
      }
    }
  } else {
    for (int j = 0; j < nIND; j++){
      xs[j] = 0;  
    }
  }
  
  *maf_i = af;
  return(na);
}

//' Covariance for loaded genotype data
//' 
//' \code{cmpt_cov_file} computes the covariance matrix of a genotype matrix.
//' 
//' @param input a genotype matrix.
//' @param min_maf a value between \code{0} and \code{0.45} specifying the threshold of minor allele frequencies above which p-values are computed.
//' @param ploidy an integer specifying the ploidy of the individuals.
//' 
//' @return The returned value is a Rcpp::List containing the covariance matrix, the number of individuals and the number of genetic markers present in the data.
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List cmpt_cov_matrix(NumericMatrix input, double min_maf, int ploidy){
  int nIND = input.ncol();
  int nSNP = input.nrow();
  NumericMatrix xs(nSNP, nIND);
  NumericMatrix xcov(nIND, nIND);
  double mean;
  double var;
  double af;
  int na;
  for (int i = 0; i < nSNP; i++){
    mean = 0;
    var = 0;
    af = 0;
    na = 0;
    for (int j = 0; j < nIND; j++){
      if ((input(i, j) != NA) && (!NumericVector::is_na(input(i, j)))){
        mean += input(i, j); 
      } else {
        na += 1;
      }
    }
    
    if (na >= nIND){
      Rcpp::stop("Detected SNP with missing values only, please remove it before proceeding.");  
    } else {
      mean /= (nIND - na);
      if (ploidy == 2){
        af = mean / 2.0;
        var = 2.0 * af * (1 - af);
      } else {
        af = mean;
        var = af * (1 - af);
      }
      if (af > 0.5){
        double tmp_af = af; 
        af = 1.0 - tmp_af;
      }
    }
    
    if (af >= min_maf){
      for (int j = 0; j < nIND; j++){
        if ((input(i, j) != NA) && (!NumericVector::is_na(input(i, j)))){
          if (var > 0){
            xs(i, j) = input(i, j) - mean;
            xs(i, j) /= sqrt(var);
          }
        } else {
          xs(i, j) = 0;
        }
      }
    } else {
      for (int j = 0; j < nIND; j++){
        xs(i, j) = 0;  
      }
    }
  }
  xcov = tAA_cpp(xs, nSNP, nIND);
  return Rcpp::List::create(Rcpp::Named("xcov") = xcov,
                            Rcpp::Named("nIND") = nIND,
                            Rcpp::Named("nSNP") = nSNP);
}

//' Linear regression
//' 
//' \code{lrfunc_matrix} performs the multiple linear regression of the genotype matrix on the scores.
//' 
//' @param Geno a genotype matrix.
//' @param scores a matrix containing the scores.
//' @param nIND an integer specifying the number of individuals present in the data.
//' @param nSNP an integer specifying the number of genetic markers present in the data.
//' @param K an integer specifying the number of principal components to retain.
//' @param ploidy an integer specifying the ploidy of the individuals.
//' @param min_maf a value between \code{0} and \code{0.45} specifying the threshold of minor allele frequencies above which p-values are computed.
//' 
//' @return The returned value is a Rcpp::List containing the multiple linear regression z-scores, the minor allele frequencies and the number of missing values for each genetic marker.
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List lrfunc_matrix(NumericMatrix Geno, NumericMatrix scores, int nIND, int nSNP, int K, int ploidy, double min_maf){
  double maf_i;
  double *residuals = new double[nSNP]();
  double *Y = new double[nIND]();
  double *GenoRowScale = new double[nIND]();
  double *sum_scores_sq = new double[K]();
  int *check_na = new int[nIND]();
  NumericMatrix Z(nSNP, K);
  IntegerVector missing(nSNP);
  NumericVector maf(nSNP);
  
  /* Linear regression */
  for (int i = 0; i < nSNP; i++){
    for (int j = 0; j < nIND; j ++){
      GenoRowScale[j] = Geno(i, j);
    }
    missing[i] = get_rows_matrix(GenoRowScale, nIND, ploidy, min_maf, &maf_i);
    maf[i] = maf_i;
    residuals[i] = 0;
    for (int k = 0; k < K; k++){
      for (int j = 0; j < nIND; j++){
        Z(i, k) += GenoRowScale[j] * scores(j, k);
      }
    }
    for (int j = 0; j < nIND; j++){
      for (int k = 0; k < K; k++){
        Y[j] += Z(i, k) * scores(j, k);
      }
    }
    for (int j = 0; j < nIND; j++){
      residuals[i] += (GenoRowScale[j] - Y[j]) * (GenoRowScale[j] - Y[j]);
      Y[j] = 0;
    }
    if ((nIND - K - missing[i]) <= 0){
      residuals[i] = 0.0;
    } else {
      residuals[i] /= nIND - K - missing[i];
    }
  }
  
  /* Correcting for missing values */
  for (int i = 0; i < nSNP; i++){
    for (int k = 0; k < K; k++){
      sum_scores_sq[k] = 0;
    }
    for (int k = 0; k < K; k++){
      for (int j = 0; j < nIND; j++){
        if (k == 0){
          if ((Geno(i, j) != NA) && !NumericVector::is_na(Geno(i, j))){
            sum_scores_sq[k] += scores(j, k) * scores(j, k);
            check_na[j] = 1;
          } else {
            check_na[j] = 0;
          }
        } else {
            if (check_na[j] == 1){
              sum_scores_sq[k] += scores(j, k) * scores(j, k);
            }
        } 
      }
      if (residuals[i] == 0.0){
        Z(i, k) = NA_REAL;
      } else {
        Z(i, k) /= sqrt(residuals[i]); 
      }
      if (sum_scores_sq[k] > 0){
        Z(i, k) /= sqrt(sum_scores_sq[k]);
      }
    }
  }
  delete[] residuals;
  delete[] Y;
  delete[] GenoRowScale;
  delete[] sum_scores_sq;
  delete[] check_na;
  return Rcpp::List::create(Rcpp::Named("zscores") = Z,
                            Rcpp::Named("maf") = maf,
                            Rcpp::Named("missing") = missing);
}
