#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

#define NA 9

using namespace Rcpp;

//' Median
//' 
//' \code{median_row_i} computes the median of a specific row.
//' 
//' @param x a genotype matrix.
//' @param i an integer.
//' 
//' @return The returned value is a real number.
//' 
//' @export
//' 
// [[Rcpp::export]]
double median_row_i(const arma::mat &x, int i){
  int nIND = x.n_cols;
  NumericVector row_i;
  for (int j = 0; j < nIND; j++){
    if (!NumericVector::is_na(x.at(i, j)) && (x.at(i, j) != NA)){
      row_i.push_back(x.at(i, j));
    }      
  }
  return(Rcpp::median(row_i));
}

//' Median
//' 
//' \code{median_per_pop} computes the median of each population for a specific row.
//' 
//' @param x a genotype matrix.
//' @param lab a vector of integers.
//' @param pop a vector of integers.
//' @param i an integer.
//' 
//' @return The returned value is a real-valued vector.
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericVector median_per_pop(const arma::mat &x, const arma::vec &lab, const arma::vec &pop, int i){
  int nIND = x.n_cols;
  int nPOP = pop.n_elem;
  NumericVector out(nPOP);
  for (int k = 0; k < nPOP; k++){
    NumericVector row_i_pop_k;
    for (int j = 0; j < nIND; j++){
      if (lab[j] == pop[k]){
        if (!NumericVector::is_na(x.at(i, j)) && (x.at(i, j) != NA)){
          row_i_pop_k.push_back(x.at(i, j));
        }
      }
    }
    out[k] = Rcpp::median(row_i_pop_k);
  }
  return(out);
}

//' Skip or discard 
//' 
//' \code{check_row} returns 0 for markers to be kept and 1 for markers to be discarded.
//' 
//' @param x a genotype matrix.
//' @param i an integer.
//' 
//' @return The returned value is an integer.
//' 
//' @export
//' 
// [[Rcpp::export]]
int check_row(const arma::mat &x, int i){
  int ncol = x.n_cols;  
  double fe = NA_REAL;
  int counter_fe = 0;
  int count = 0;
  int skip = 0;
  while ((NumericVector::is_na(fe) || (fe == NA)) && (counter_fe < ncol)){
    fe = x(i, counter_fe);
    counter_fe++;
  } 
  if (counter_fe == ncol){
    skip = 1;
  } else {
    for (int j = 0; j < ncol; j++){
      if ((x(i, j) == fe) || (NumericVector::is_na(x(i, j))) || (x(i, j) == NA)){
        count++;
      }
    }
    if (count == ncol){
      skip = 1;
    }
  }
  return(skip);
}

//' Genotype matrix imputation
//' 
//' \code{impute_geno} imputes values based on the median per row.
//' 
//' @param x a genotype matrix.
//' 
//' @return The returned value is a list containing the imputed genotype matrix and a vector indicating which markers 
//' have to be discarded.
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List impute_geno(const arma::mat &x){
  int nrow = x.n_rows;
  int ncol = x.n_cols;
  arma::vec skip(nrow);
  arma::mat y = x;
  for (int i = 0; i < nrow; i++){
    skip[i] = check_row(x, i);
    for (int j = 0; j < ncol; j++){
      if (NumericVector::is_na(x.at(i, j)) || (x.at(i, j) == NA)){
        y.at(i, j) = median_row_i(x, i);  
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("x") = y,
                            Rcpp::Named("skip") = skip);
}

//' Genotype matrix imputation
//' 
//' \code{impute_geno_pop} imputes values based on the median per row and per population.
//' 
//' @param x a genotype matrix.
//' @param lab a vector of integers.
//' @param pop a vector of integers.
//' 
//' @return The returned value is a list containing the imputed genotype matrix and a vector indicating which markers 
//' have to be discarded.
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List impute_geno_pop(const arma::mat &x, const arma::vec &lab, const arma::vec &pop){
  int nSNP = x.n_rows;
  int nIND = x.n_cols;
  int nPOP = pop.n_elem;
  arma::vec skip(nSNP);
  arma::mat y(nSNP, nIND);
  
  for (int i = 0; i < nSNP; i++){
    skip[i] = check_row(x, i);
    NumericVector mpp_i = median_per_pop(x, lab, pop, i);
    double m_i = median_row_i(x, i);
    for (int j = 0; j < nIND; j++){
      if (NumericVector::is_na(x.at(i, j)) || (x.at(i, j) == NA)){
        for (int k = 0; k < nPOP; k++){
          if (lab[j] == pop[k]){
            if (!NumericVector::is_na(mpp_i[k])){
              y.at(i, j) = mpp_i[k];    
            } else if (!NumericVector::is_na(m_i)){
              y.at(i, j) = m_i;    
            } else {
              y.at(i, j) = 0.0;
              Rprintf("SNP %d is not informative. No variation or too many missing values.\n", i);
            }
          }
        }
      } else {
        y.at(i, j) = x.at(i, j);
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("x") = y,
                            Rcpp::Named("skip") = skip);
}


