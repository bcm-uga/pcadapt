/******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

#define NA 9

/******************************************************************************/

//' Sample genotype matrix from pooled samples
//' 
//' \code{sample_geno_matrix} sample genotypes based on observed allelic frequencies.
//' 
//' @param freq a matrix containing allele frequencies.
//' @param ploidy an integer specifying the ploidy of the sampled individuals.
//' @param sample_size a vector specifying the number of individuals to be 
//'   sampled for each pool.
//' 
//' @return The returned value is a numeric vector of length 2.
//' 
// [[Rcpp::export]]
NumericMatrix sample_geno_matrix(const NumericMatrix& freq, 
                                 double ploidy, 
                                 const IntegerVector& sample_size){
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
      std::transform(prob.begin(), prob.end(), v.begin(), 
                     [=](double p){return R::rbinom(ploidy, p);}); 
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

/******************************************************************************/

//' Sample genotype matrix from pooled samples
//' 
//' \code{sample_geno_file} sample genotypes based on observed allelic frequencies.
//' 
//' @param input a character string specifying the name of the file containing 
//'   the allele frequencies.
//' @param output a character string specifying the name of the output file.
//' @param ploidy an integer specifying the ploidy of the sampled individuals.
//' @param sample_size a vector specifying the number of individuals to be 
//'   sampled for each pool.
//' 
//' @return The returned value is a numeric vector of length 2.
//' 
// [[Rcpp::export]]
NumericVector sample_geno_file(std::string input, 
                               std::string output, 
                               double ploidy, 
                               const IntegerVector& sample_size){
  FILE *file_in;
  file_in = fopen(input.c_str(), "r");
  FILE *file_out;
  file_out = fopen(output.c_str(), "w");
  NumericVector file_size = get_size_cpp(input);
  int nSNP = file_size[1];
  int nPOOL = file_size[0];
  IntegerVector na(nSNP);
  NumericVector v = no_init(nSNP);
  NumericVector prob = no_init(nSNP);
  float value;
  
  for (int k = 0; k < nPOOL; k++){
    for (int fill = 0; fill < nSNP; fill++){
      if (fscanf(file_in, "%g", &value) != EOF){
        if ((value != 9.0) || !NumericVector::is_na(value)){
          prob[fill] = (double) value;
          na[fill] = 0;
        } else {
          prob[fill] = 0.0;
          na[fill] = 1;
        }
      }
    }
    for (int j = 0; j < sample_size[k]; j++){
      std::transform(prob.begin(), prob.end(), v.begin(), 
                     [=](double p){return R::rbinom(ploidy, p);}); 
      for (int i = 0; i < nSNP; i++){
        int tmp = v[i];
        if (i < (nSNP - 1)){
          if (na[i] != 1){
            fprintf(file_out, "%d ", tmp);
          } else {
            fprintf(file_out, "%d ", 9);
          }
        } else if (i == (nSNP - 1)){
          if (na[i] != 1){
            fprintf(file_out, "%d", tmp);
          } else {
            fprintf(file_out, "%d", 9);
          }  
        }
      }
      fprintf(file_out, "\n");
    }
  }
  fclose(file_in);
  fclose(file_out);
  return(file_size);
}

/******************************************************************************/