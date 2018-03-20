#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int count_na1(NumericVector x) {
  
  int n = x.size();
  int c = 0;
  
  for (int i = 0; i < n; i++) {
    if (R_IsNA(x[i])) c++;
  }
  
  return c;
}


// [[Rcpp::export]]
int count_na2(NumericVector x) {
  
  int n = x.size();
  int c = 0;
  double x_i;
  
  for (int i = 0; i < n; i++) {
    x_i = x[i];
    if (x_i == 0) {
      
    } else if (x_i == 1) {
      
    } else if (x_i == 2) {
      
    } else {
      c++;
    }
  }
  
  return c;
}

/*** R
N <- 1e5
x1 <- sample(0:2, size = N, replace = TRUE)
x2 <- sample(c(0:2, NA), size = N, replace = TRUE)
microbenchmark::microbenchmark(
  count_na1(x1), count_na2(x1),
  count_na1(x2), count_na2(x2)
)
*/
