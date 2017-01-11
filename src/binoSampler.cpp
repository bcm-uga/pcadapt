#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericMatrix sample_geno_cpp(NumericMatrix freq, double ploidy, NumericVector sample_size, int nIND){
  int nPOOL = freq.nrow();
  int nSNP = freq.ncol();
  NumericMatrix m = no_init(nSNP, nIND);
  NumericVector v = no_init(nSNP);
  NumericVector prob = no_init(nSNP);
  int current_position = 0;
  for (int k = 0; k < nPOOL; k++){
    for (int fill = 0; fill < nSNP; fill++){
      prob[fill] = freq(k, fill);
    }
    for (int j = 0; j < sample_size[k]; j++){
      std::transform(prob.begin(), prob.end(), v.begin(), [=](double p){return R::rbinom(ploidy, p);}); 
      for (int i = 0; i < nSNP; i++){
        m(i, current_position) = v[i];
      }
      current_position++;
    }
  }
  return(m);
}

