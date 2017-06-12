#include <RcppArmadillo.h>
#include "toolbox.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

#define NA 9

void add_to_cov_cpp(arma::mat &cov, arma::mat &genoblock){
  double tmp = 0;
  for (int i = 0; i < cov.n_rows; i++){
    for (int j = i; j < cov.n_rows; j++){
      if (i != j){
        tmp = arma::dot(genoblock.col(i), genoblock.col(j));
        cov(i, j) += tmp;
        cov(j, i) += tmp;
      } else {
        tmp = arma::dot(genoblock.col(i), genoblock.col(j));
        cov(i, j) += tmp;
      }
    }
  }
}

//' Covariance for loaded genotype data
//' 
//' \code{cmpt_cov_cpp} computes the covariance matrix of a genotype matrix.
//' 
//' @param filename a character string specifying the name of the file to be 
//' processed with \code{pcadapt}.
//' @param xmatrix a genotype matrix.
//' @param min_maf a value between \code{0} and \code{0.45} specifying the 
//' threshold of minor allele frequencies above which p-values are computed.
//' @param ploidy an integer specifying the ploidy of the individuals.
//' @param type an integer specifying the input type.
//' @param blocksize an integer specifying the number of rows for each block.
//' 
//' @return The returned value is a Rcpp::List containing the covariance matrix, 
//' the number of individuals and the number of genetic markers present in the 
//' data.
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List cmpt_cov_cpp(std::string filename, arma::mat &xmatrix, 
                        double min_maf, int ploidy, int type, 
                        int blocksize = 120){
  int nSNP = 0;
  int nIND = 0;
  FILE *xfile;
  if (type == 0){
    xfile = fopen(filename.c_str(), "r");
    Rprintf("Reading file %s...\n", filename.c_str());
    NumericVector file_size = get_size_cpp(filename);
    nSNP = file_size[0];
    nIND = file_size[1];
  } else if (type == 1){
    nIND = xmatrix.n_cols;
    nSNP = xmatrix.n_rows;
  }
  Rprintf("Number of SNPs: %i\n", nSNP);
  Rprintf("Number of individuals: %i\n", nIND);  
  arma::mat scale_geno(nSNP, nIND);
  int unused_na = 0;
  double unused_maf = 0;
  int b;
  arma::mat xcov(nIND, nIND, arma::fill::zeros);
  arma::vec unused_missing(nIND, arma::fill::zeros);
  
  for (int i = 0; i < nSNP; i += blocksize){
    if (nSNP - i < blocksize){
      b = nSNP - i;
    } else {
      b = blocksize;
    }
    arma::mat geno(b, nIND, arma::fill::zeros);
    if (type == 0){
      unused_na = scale_rows(xfile, xmatrix, geno, nIND, ploidy, min_maf, 0, b, 
                             &unused_maf, unused_missing, 0);
    } else if (type == 1){
      unused_na = scale_rows(xfile, xmatrix, geno, nIND, ploidy, min_maf, i, 
                             (i + b), &unused_maf, unused_missing, 1); 
    }
    add_to_cov_cpp(xcov, geno);
  }
  
  if (type == 0){
    fclose(xfile);
  }
  
  return Rcpp::List::create(Rcpp::Named("xcov") = xcov,
                            Rcpp::Named("nIND") = nIND,
                            Rcpp::Named("nSNP") = nSNP);
}
