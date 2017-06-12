#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

#define NA 9

NumericVector get_size_cpp(std::string filename);

int get_nb_ind(const arma::vec &lab, const int anc);

NumericVector cmpt_minor_af(arma::mat &xmatrix, int ploidy);

arma::mat scale_geno(arma::mat &xmatrix, int ploidy, arma::vec maf, 
                     arma::vec keep_or_not);

int scale_rows(FILE *xfile, arma::mat &xmatrix, arma::mat &xs, int nIND, 
               int ploidy, double min_maf, int begin, int end, double *maf_i, 
               arma::vec &missing_i, int type);

void add_to_cov_cpp(arma::mat &cov, arma::mat &genoblock);

Rcpp::List cmpt_cov_cpp(std::string filename, arma::mat &xmatrix, 
                        double min_maf, int ploidy, int type, int blocksize);

arma::mat cmpt_loadings(std::string filename, arma::mat &xmatrix, 
                        arma::mat &scores, int nIND, int nSNP, int K, 
                        int ploidy, double min_maf, arma::vec &sigma, int type);

Rcpp::List lrfunc_cpp(std::string filename, arma::mat &xmatrix, 
                      arma::mat &scores, int nIND, int nSNP, int K, int ploidy, 
                      double min_maf, arma::vec &sigma, int type);

NumericVector sample_geno_file(std::string input, std::string output, 
                               double ploidy, IntegerVector sample_size);

