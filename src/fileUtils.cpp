#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

#define NA 9

//' File size
//' 
//' \code{get_size_file} returns the number of genetic markers and the number of individuals present in the data.
//' 
//' @param path a character string specifying the name of the file to be processed with \code{pcadapt}.
//' 
//' @return The returned value is a numeric vector of length 2.
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericVector get_size_file(std::string path){
  NumericVector file_size(2);
  FILE *input;
  int nbbn = 0;
  int nbsp = 0;
  if ((input = fopen(path.c_str(), "r")) == NULL){
    Rprintf("Error, invalid input file.\n");
  }
  int currentchar;
  int prevchar;
  currentchar = fgetc(input);
  while(currentchar != EOF){
    if (currentchar == 10){
      nbbn++;
      if (prevchar != 32 && prevchar != '\t'){
        nbsp++;
      }
    }
    if ((currentchar == 32 || prevchar == '\t') && (prevchar != 32 || prevchar != '\t')){
      nbsp++;
    }
    prevchar = currentchar;
    currentchar = fgetc(input);
  }
  fclose(input);
  file_size[0] = nbbn;
  file_size[1] = nbsp/nbbn;
  return file_size;
}

void tr(double *x, int nrow, int ncol){
  double *tmp = new double[nrow * ncol];
  for (int i = 0; i < ncol * nrow; i++){
    tmp[i] = x[i];
  }
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      x[j * nrow + i] = tmp[ i * ncol + j];
    }
  }
  delete[] tmp;
}

void tAA_arma(double *A, double *tAA, int nrow, int ncol){
  double *tA = new double[nrow * ncol];
  for (int i = 0; i < nrow * ncol; i++){
    tA[i] = A[i];
  }
  arma::blas_int nrtA = ncol, nctB = ncol, nrtB = nrow;
  arma::blas_int ldc = nrtA, lda = nrtA, ldb = nrtB;
  double alpha = 1, beta = 0;
  tr(A, nrow, ncol);
  arma::dgemm_("N", "N", &nrtA, &nctB, &nrtB, &alpha, tA, &lda, A, &ldb, &beta, tAA, &ldc);
  tr(A, ncol, nrow);
  delete[] tA;
}

void add_to_cov(double *xcov, double *scratchcov, int nIND, double *geno, int blocksize){
  tAA_arma(geno, scratchcov, blocksize, nIND);
  for (int i = 0; i < nIND*nIND; i++){
    xcov[i] += scratchcov[i];
  }
}

int get_rows_file(double *xs, FILE *input, int nIND, int ploidy, double min_maf, int blocksize, double *maf_i){
  double mean = 0;
  double var = 0;
  double af = 0;
  int na = 0;
  float value;
  
  for (int i = 0; i < blocksize; i++){
    var = 0;
    mean = 0;
    na = 0;
    for (int j = 0; j < nIND; j++){
      if (fscanf(input, "%g", &value) != EOF){
        xs[j + i * nIND] = (double) value;
        if (value != NA){
          mean += (double) value;
        } else {
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
        if (xs[j + i * nIND] != NA){
          if (var > 0){
            xs[j + i * nIND] -= mean;
            xs[j + i * nIND] /= sqrt(var);
          }
        } else {
          xs[j + i * nIND] = 0; 
        }
      }
    } else {
      for (int j = 0; j < nIND; j++){
        xs[j + i * nIND] = 0;  
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

//' Covariance for genotype data stored in an external file
//' 
//' \code{cmpt_cov_file} computes the covariance matrix of a genotype matrix when the genotype matrix is stored in an external file.
//' 
//' @param path a character string specifying the name of the file to be processed with \code{pcadapt}.
//' @param min_maf a value between \code{0} and \code{0.45} specifying the threshold of minor allele frequencies above which p-values are computed.
//' @param ploidy an integer specifying the ploidy of the individuals.
//' 
//' @return The returned value is a Rcpp::List containing the covariance matrix, the number of individuals and the number of genetic markers present in the data.
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List cmpt_cov_file(std::string path, double min_maf, int ploidy){
  FILE *input;
  input = fopen(path.c_str(), "r");
  Rprintf("Reading file %s...\n",path.c_str());
  NumericVector file_size = get_size_file(path);
  int nSNP = file_size[0];
  int nIND = file_size[1];
  Rprintf("Number of SNPs: %i\n", nSNP);
  Rprintf("Number of individuals: %i\n", nIND);
  
  int unused_na = 0;
  double unused_maf = 0;
  int blocksize = 100;
  
  double *scratchcov = new double[nIND * nIND]();
  double *xcov_ptr = new double[nIND * nIND]();
  double *geno_ptr = new double[blocksize * nIND]();
  for (int i = 0; i < nSNP; i += blocksize){
    if (nSNP - i < blocksize){
      blocksize = nSNP - i;
    }
    unused_na = get_rows_file(geno_ptr, input, nIND, ploidy, min_maf, blocksize, &unused_maf);
    add_to_cov(xcov_ptr, scratchcov, nIND, geno_ptr, blocksize);
  }
  NumericMatrix xcov(nIND, nIND);
  for (int k = 0; k < nIND; k++){
    for (int l = 0; l < nIND; l++){
      xcov(k, l) = xcov_ptr[l + k * nIND];
    }
  }
  
  fclose(input);
  delete[] scratchcov;
  delete[] xcov_ptr;
  delete[] geno_ptr;
  return Rcpp::List::create(Rcpp::Named("xcov") = xcov,
                            Rcpp::Named("nIND") = nIND,
                            Rcpp::Named("nSNP") = nSNP);
}


//' Linear regression
//' 
//' \code{lrfunc_file} performs the multiple linear regression of the genotype matrix on the scores when the genotype matrix is stored in an external file.
//' 
//' @param filename a character string specifying the name of the file to be processed with \code{pcadapt}.
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
Rcpp::List lrfunc_file(std::string filename, NumericMatrix scores, int nIND, int nSNP, int K, int ploidy, double min_maf){
  FILE *input;
  input = fopen(filename.c_str(), "r");
  double maf_i;
  double *residuals = new double[nSNP]();
  double *Y = new double[nIND]();
  double *GenoRowScale = new double[nIND]();
  double *sum_scores_sq = new double[K]();
  int *check_na = new int[nIND]();
  NumericMatrix Z(nSNP, K);
  NumericVector missing(nSNP);
  NumericVector maf(nSNP);
  float value;
  
  for (int i = 0; i < nSNP; i++){
    missing[i] = get_rows_file(GenoRowScale, input, nIND, ploidy, min_maf, 1, &maf_i);
    maf[i] = maf_i;
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
      residuals[i] += (GenoRowScale[j] - Y[j])*(GenoRowScale[j] - Y[j]);
      Y[j] = 0;
    }
    if ((nIND - K - missing[i]) <= 0){
      residuals[i] = 0.0;
    } else {
      residuals[i] /= nIND - K - missing[i];
    }
  }
  fclose(input);
  
  /* Correcting for missing values */
  input = fopen(filename.c_str(), "r");
  
  for (int i = 0; i < nSNP; i++){
    for (int l = 0; l < K; l++){
      sum_scores_sq[l] = 0;
    }
    for (int k = 0; k < K; k++){
      for (int j = 0; j < nIND; j++){
        if (k == 0){
          if (fscanf(input, "%g", &value) != EOF){
            if (value != NA){
              sum_scores_sq[k] +=  (double) scores(j ,k) * scores(j, k);
              check_na[j] = 1;
            } else {
              check_na[j] = 0;
            }
          } 
        } else {
          if (check_na[j] == 1){
            sum_scores_sq[k] += (double) scores(j, k) * scores(j, k);
          }
        }
      }
      if (residuals[i] == 0.0){
        Z(i, k) = 0.0;
      } else {
        Z(i ,k) /= sqrt(residuals[i]);
      }
      if (sum_scores_sq[k] > 0){
        Z(i ,k) /= sqrt(sum_scores_sq[k]);
      }
    }
  }
  fclose(input);
  delete[] residuals;
  delete[] Y;
  delete[] GenoRowScale;
  delete[] sum_scores_sq;
  delete[] check_na;
  return Rcpp::List::create(Rcpp::Named("zscores") = Z,
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
  NumericVector file_size = get_size_file(input);
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













