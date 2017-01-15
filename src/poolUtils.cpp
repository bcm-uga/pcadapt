#include <Rcpp.h>
using namespace Rcpp;

#define NA 9

//' Sample genotype matrix from pooled samples
//' 
//' \code{sample_geno_cpp} sample genotypes based on observed allelic frequencies.
//' 
//' @param freq a matrix containing allele frequencies.
//' @param ploidy an integer specifying the ploidy of the sampled individuals.
//' @param sample_size a vector specifying the number of individuals to be sampled for each pool.
//' 
//' @return The returned value is a numeric vector of length 2.
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericMatrix sample_geno_cpp(NumericMatrix freq, double ploidy, IntegerVector sample_size){
  int nIND = Rcpp::sum(sample_size);
  int nPOOL = freq.nrow();
  int nSNP = freq.ncol();
  IntegerVector na(nSNP);
  NumericMatrix m = no_init(nSNP, nIND);
  NumericVector v = no_init(nSNP);
  NumericVector prob = no_init(nSNP);
  int current_position = 0;
  for (int k = 0; k < nPOOL; k++){
    for (int fill = 0; fill < nSNP; fill++){
      if ((freq(k, fill) != NA) && !NumericVector::is_na(freq(k, fill))){
        prob[fill] = freq(k, fill);
        na[fill] = 0;
      } else {
        prob[fill] = 0.0;
        na[fill] = 1;
      }
    }
    for (int j = 0; j < sample_size[k]; j++){
      std::transform(prob.begin(), prob.end(), v.begin(), [=](double p){return R::rbinom(ploidy, p);}); 
      for (int i = 0; i < nSNP; i++){
        if (na[i] != 1){
          m(i, current_position) = v[i];
        } else {
          m(i, current_position) = NA;
        }
      }
      current_position++;
    }
  }
  return(m);
}