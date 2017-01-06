#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_row_cpp(NumericVector x, int n, int ploidy, double min_maf){
  NumericVector xs(n);
  double mean = 0;
  double var = 0;
  double maf = 0;
  int na = 0;
  
  for (int j = 0; j < n; j++){
    if (!NumericVector::is_na(x[j])){
      mean += x[j]; 
    } else {
      na += 1;
    }
  }
  
  if (na >= n){
    mean = NA_REAL;
    maf = NA_REAL;
  } else {
    mean /= (n - na);
    if (ploidy == 2){
      maf = mean / 2.0;
      var = 2.0 * maf * (1 - maf);
    } else {
      maf = mean;
      var = maf * (1 - maf);
    }
    if (maf > 0.5){
      double tmp_maf = maf; 
      maf = 1.0 - tmp_maf;
    }
  }
  
  if (maf >= min_maf){
    for (int j = 0; j < n; j++){
      if (!NumericVector::is_na(x[j])){
        if (var > 0){
          xs[j] = (x[j] - mean) / sqrt(var);
        }
      } else {
        xs[j] = 0;
      }
    }
  } else {
    for (int j = 0; j < n; j++){
      xs[j] = 0;  
    }
  }
  return(xs);
}

//' Compute graph matrix with the heat kernel function.
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix lrfunc_cpp(NumericMatrix scores, NumericMatrix Geno, int K, int nIND, int nSNP, int ploidy, double min_maf){
  
  NumericMatrix Z(nSNP, K);
  NumericMatrix U(nIND, K);
  NumericMatrix t_U(K, nIND);
  NumericVector residuals(nSNP);
  NumericVector Y(nIND);
  NumericVector na(nIND);
  NumericVector GenoRow(nIND);
  NumericVector GenoRowScale(nIND);
  
  for (int a = 0; a < K; a++){
    for (int b = 0; b < nIND; b++){
      U(b, a) = scores(b, a);
      t_U(a, b) = scores(b, a);
    }
  }
  
  // Linear regression //
  for (int i = 0; i < nSNP; i++){
    for (int j = 0; j < nIND; j ++){
      GenoRow[j] = Geno(i, j);
    }
    GenoRowScale = get_row_cpp(GenoRow, nIND, ploidy, min_maf);
    residuals[i] = 0;
    for (int k = 0; k < K; k++){
      for (int j = 0; j < nIND; j++){
        Z(i, k) += GenoRowScale[j]*U(j, k);
      }
    }
    for (int j = 0; j < nIND; j++){
      for (int k = 0; k < K; k++){
        Y[j] += Z(i, k) * t_U(k, j);
      }
    }
    for (int j = 0; j < nIND; j++){
      residuals[i] += (GenoRowScale[j] - Y[j])*(GenoRowScale[j] - Y[j]);
      Y[j] = 0;
    }
    if (nIND - K <= 0){
      residuals[i] = 0.0;
    } else {
      residuals[i] /= nIND - K;
    }
  }
  
  // Correcting for missing values //
  NumericVector sum_scores_sq(K);
  NumericVector check_na(nIND);
  for (int i = 0; i < nSNP; i++){
    for (int k = 0; k < K; k++){
      sum_scores_sq[k] = 0;
    }
    for (int k = 0; k < K; k++){
      for (int j = 0; j < nIND; j++){
        if (k == 0){
          if (!NumericVector::is_na(Geno(i, j))){
            sum_scores_sq[k] += U(j, k) * U(j, k);
            check_na[j] = 1;
          } else {
            check_na[j] = 0;
          }
        } else {
            if (check_na[j] == 1){
              sum_scores_sq[k] += U(j, k) * U(j, k);
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
  return Z;
}



