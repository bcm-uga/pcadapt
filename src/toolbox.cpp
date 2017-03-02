#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

#define NA 9

//' File size
//' 
//' \code{get_size_cpp} returns the number of genetic markers and the number of individuals present in the data.
//' 
//' @param filename a character string specifying the name of the file to be processed with \code{pcadapt}.
//' 
//' @return The returned value is a numeric vector of length 2.
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericVector get_size_cpp(std::string filename){
  NumericVector file_size(2);
  FILE *input;
  if ((input = fopen(filename.c_str(), "r")) == NULL){
    Rprintf("Error, invalid input file.\n");
  }
  int currentchar;
  int nrow = 0;
  int ncol = 0;
  currentchar = fgetc(input);
  while(currentchar != EOF){
    if (nrow == 0 && currentchar != '\n' && currentchar != '\r' && currentchar != ' '){
      ncol ++;
    }
    if (currentchar == '\n'){
      nrow ++;
    }
    currentchar = fgetc(input);
  }
  fclose(input);
  file_size[0] = nrow;
  file_size[1] = ncol;
  return file_size;
}

//' Minor allele frequencies
//' 
//' \code{cmpt_minor_af} computes the minor allele frequencies.
//' 
//' @param xmatrix a genotype matrix.
//' @param ploidy an integer specifying the ploidy of the individuals.
//' 
//' @return The returned value is a numeric matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericVector cmpt_minor_af(arma::mat &xmatrix, int ploidy){
  double mean = 0;
  double var = 0;
  double af = 0;
  int na = 0;
  int nSNP = xmatrix.n_rows;
  int nIND = xmatrix.n_cols;
  NumericVector maf(nSNP);
  for (int i = 0; i < nSNP; i++){
    var = 0;
    mean = 0;
    na = 0;
    for (int j = 0; j < nIND; j++){
      if ((xmatrix.at(i, j) != NA) && (!NumericVector::is_na(xmatrix.at(i, j)))){
        mean += xmatrix.at(i, j); 
      } else {
        na += 1;
      }
    }
    if (na >= nIND){
      af = NA_REAL;
      //Rcpp::stop("Detected SNP with missing values only, please remove it before proceeding."); 
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
    maf[i] = af;
  }
  return(maf);
}

//' Scale genotype matrices
//' 
//' \code{scale_geno} scales the genotype matrix.
//' 
//' @param xmatrix a genotype matrix.
//' @param ploidy an integer specifying the ploidy of the individuals.
//' @param maf a vector of minor allele frequencies.
//' @param keep_or_not a vector of integers.
//' 
//' @return The returned value is a numeric matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::mat scale_geno(arma::mat &xmatrix, int ploidy, arma::vec maf, arma::vec keep_or_not){
  int nSNP = xmatrix.n_rows;
  int nIND = xmatrix.n_cols;
  int nSNP_kept = arma::sum(keep_or_not);
  arma::mat geno(nSNP_kept, nIND);
  int current_position = 0;
  double mean = 0.0;
  double sd = 1.0;
  for (int i = 0; i < nSNP; i++){
    mean = 0.0;
    sd = 1.0;
    if (keep_or_not[i] == 1){
      mean = arma::mean(xmatrix.row(i));
      mean /= ploidy;
      sd = sqrt(ploidy * maf[i] * (1 - maf[i]));
      geno.row(current_position) = xmatrix.row(i);
      geno.row(current_position) -= mean;
      geno.row(current_position) /= sd;
      current_position++;
    }
  }
  return(geno);
}

int scale_rows(FILE *xfile, arma::mat &xmatrix, arma::mat &xs, int nIND, int ploidy, double min_maf, int begin, int end, double *maf_i, arma::vec &missing_i, int type){
  double mean = 0;
  double var = 0;
  double af = 0;
  int na = 0;
  int blocksize = end - begin;
  
  for (int i = 0; i < blocksize; i++){
    var = 0;
    mean = 0;
    na = 0;
    
    if (type == 0){
      float value;
      for (int j = 0; j < nIND; j++){
        if (fscanf(xfile, "%g", &value) != EOF){
          xs(i, j) = (double) value;
          if (value != NA){
            mean += (double) value;
          } else {
            missing_i[j] = 1.0;
            na += 1;
          }
        }
      }
    } else if (type == 1){
      for (int j = 0; j < nIND; j++){
        xs(i, j) = xmatrix((begin + i), j);
        if ((xs(i, j) != NA) && (!NumericVector::is_na(xs(i, j)))){
          mean += xs(i, j); 
        } else {
          missing_i[j] = 1.0;
          na += 1;
        }
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
        if ((xs(i, j) != NA) && (!NumericVector::is_na(xs(i, j)))){
          if (var > 0){
            xs(i, j) -= mean;
            xs(i, j) /= sqrt(var);
          }
        } else {
          xs(i, j) = 0.0; 
        }
      }
    } else {
      for (int j = 0; j < nIND; j++){
        xs(i, j) = 0.0;  
      }
    }
  }
  
  if (blocksize == 1){
    *maf_i = af;
  } else {
    *maf_i = 0;
  }
  
  return(na);
}

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
//' @param filename a character string specifying the name of the file to be processed with \code{pcadapt}.
//' @param xmatrix a genotype matrix.
//' @param min_maf a value between \code{0} and \code{0.45} specifying the threshold of minor allele frequencies above which p-values are computed.
//' @param ploidy an integer specifying the ploidy of the individuals.
//' @param type an integer specifying the input type.
//' 
//' @return The returned value is a Rcpp::List containing the covariance matrix, the number of individuals and the number of genetic markers present in the data.
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List cmpt_cov_cpp(std::string filename, arma::mat &xmatrix, double min_maf, int ploidy, int type){
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
  int blocksize = 14;
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
      unused_na = scale_rows(xfile, xmatrix, geno, nIND, ploidy, min_maf, 0, b, &unused_maf, unused_missing, 0);
    } else if (type == 1){
      unused_na = scale_rows(xfile, xmatrix, geno, nIND, ploidy, min_maf, i, (i + b), &unused_maf, unused_missing, 1); 
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

//' Compute the loadings
//' 
//' \code{cmpt_loadings} returns the loadings.
//' 
//' @param filename a character string specifying the name of the file to be processed with \code{pcadapt}.
//' @param xmatrix a genotype matrix.
//' @param scores a matrix containing the scores.
//' @param nIND an integer specifying the number of individuals present in the data.
//' @param nSNP an integer specifying the number of genetic markers present in the data.
//' @param K an integer specifying the number of principal components to retain.
//' @param ploidy an integer specifying the ploidy of the individuals.
//' @param min_maf a value between \code{0} and \code{0.45} specifying the threshold of minor allele frequencies above which p-values are computed.
//' @param sigma a numeric vector.
//' @param type an integer specifying the input type.
//' 
//' @return The returned value is a matrix containing the loadings.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::mat cmpt_loadings(std::string filename, arma::mat &xmatrix, arma::mat &scores, int nIND, int nSNP, int K, int ploidy, double min_maf, arma::vec &sigma, int type){
  FILE *xfile;
  
  if (type == 0){
    xfile = fopen(filename.c_str(), "r");
  } 
  
  double maf_i;
  double unused_missing = 0.0;
  
  arma::mat GenoRowScale(1, nIND, arma::fill::zeros);
  arma::mat V(nSNP, K, arma::fill::zeros);
  arma::vec unused_check_na(nIND, arma::fill::zeros);
  
  for (int i = 0; i < nSNP; i++){
    if (type == 0){
      unused_missing = scale_rows(xfile, xmatrix, GenoRowScale, nIND, ploidy, min_maf, 0, 1, &maf_i, unused_check_na, 0);  
    } else if (type == 1){
      unused_missing = scale_rows(xfile, xmatrix, GenoRowScale, nIND, ploidy, min_maf, i, i + 1, &maf_i, unused_check_na, 1);  
    }
    
    V.row(i) = GenoRowScale * scores;
    for (int k = 0; k < K; k++){
      V(i, k) *= sqrt(nSNP / sigma[k]);
    }
  }
  
  if (type == 0){
    fclose(xfile);
  } 
  
  return(V);
}

//' Linear regression
//' 
//' \code{lrfunc_cpp} performs the multiple linear regression of the genotype matrix on the scores.
//' 
//' @param filename a character string specifying the name of the file to be processed with \code{pcadapt}.
//' @param xmatrix a genotype matrix.
//' @param scores a matrix containing the scores.
//' @param nIND an integer specifying the number of individuals present in the data.
//' @param nSNP an integer specifying the number of genetic markers present in the data.
//' @param K an integer specifying the number of principal components to retain.
//' @param ploidy an integer specifying the ploidy of the individuals.
//' @param min_maf a value between \code{0} and \code{0.45} specifying the threshold of minor allele frequencies above which p-values are computed.
//' @param sigma a numeric vector.
//' @param type an integer specifying the input type.
//' 
//' @return The returned value is a Rcpp::List containing the multiple linear regression z-scores, the minor allele frequencies and the number of missing values for each genetic marker.
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List lrfunc_cpp(std::string filename, arma::mat &xmatrix, arma::mat &scores, int nIND, int nSNP, int K, int ploidy, double min_maf, arma::vec &sigma, int type){
  FILE *xfile;
  
  if (type == 0){
    xfile = fopen(filename.c_str(), "r");
  } 
  
  double maf_i;
  double residual;
  arma::mat Y(1, nIND, arma::fill::zeros);
  arma::mat GenoRowScale(1, nIND, arma::fill::zeros);
  arma::vec sum_scores_sq(K, arma::fill::zeros);
  arma::mat Z(nSNP, K, arma::fill::zeros);
  arma::mat V(nSNP, K, arma::fill::zeros);
  
  NumericVector missing(nSNP);
  NumericVector maf(nSNP);
  arma::vec check_na(nIND, arma::fill::zeros);
  
  for (int i = 0; i < nSNP; i++){
    if (type == 0){
      missing[i] = scale_rows(xfile, xmatrix, GenoRowScale, nIND, ploidy, min_maf, 0, 1, &maf_i, check_na, 0);  
    } else if (type == 1){
      missing[i] = scale_rows(xfile, xmatrix, GenoRowScale, nIND, ploidy, min_maf, i, i + 1, &maf_i, check_na, 1);  
    }
    maf[i] = maf_i;
    Z.row(i) = GenoRowScale * scores;
    for (int k = 0; k < K; k++){
      V(i, k) = Z(i, k) * sqrt(nSNP / sigma[k]);
    }
    Y = Z.row(i) * scores.t();
    residual = dot(GenoRowScale - Y, (GenoRowScale - Y).t());
    Y.zeros();
    if ((nIND - K - missing[i]) <= 0){
      residual = 0.0;
    } else {
      residual /= nIND - K - missing[i];
    }
    
    /* Correcting for missing values */
    sum_scores_sq.zeros();
    for (int k = 0; k < K; k++){
      for (int j = 0; j < nIND; j++){
        if (check_na[j] != 1.0){
          sum_scores_sq[k] +=  (double) scores(j ,k) * scores(j, k);
        }
      }
      if (residual == 0.0){
        Z(i, k) = 0.0;
      } else {
        Z(i ,k) /= sqrt(residual);
      }
      if (sum_scores_sq[k] > 0){
        Z(i ,k) /= sqrt(sum_scores_sq[k]);
      }
    }
    check_na.zeros();
  }
  
  if (type == 0){
    fclose(xfile);
  } 
  
  return Rcpp::List::create(Rcpp::Named("zscores") = Z,
                            Rcpp::Named("loadings") = V,
                            Rcpp::Named("maf") = maf,
                            Rcpp::Named("missing") = missing);
}

//' Sample genotype matrix from pooled samples
//' 
//' \code{sample_geno_file} sample genotypes based on observed allelic frequencies.
//' 
//' @param input a character string specifying the name of the file containing the allele frequencies.
//' @param output a character string specifying the name of the output file.
//' @param ploidy an integer specifying the ploidy of the sampled individuals.
//' @param sample_size a vector specifying the number of individuals to be sampled for each pool.
//' 
//' @return The returned value is a numeric vector of length 2.
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericVector sample_geno_file(std::string input, std::string output, double ploidy, IntegerVector sample_size){
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
      std::transform(prob.begin(), prob.end(), v.begin(), [=](double p){return R::rbinom(ploidy, p);}); 
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
