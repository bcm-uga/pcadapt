/******************************************************************************/

#include <pcadapt/bed-acc.h>
#include <pcadapt/mat-acc.h>

/******************************************************************************/

// 3 are coding NAs
template <class C>
LogicalVector clumping(C macc,
                       const IntegerVector& ord,
                       LogicalVector& remain,
                       int size, 
                       double thr) {
  
  int n = macc.nrow();
  int p = macc.ncol();
  LogicalVector keep(p);
  
  int i, j, k, j0, N;
  double x, y;
  double xSum, xxSum, deno_x;
  double ySum, yySum, deno_y;
  double xySum, num, r2;
  
  // pre-computation
  NumericVector sumX(p), sumXX(p);
  for (j = 0; j < p; j++) {
    
    if (remain[j]) {
      xSum = xxSum = 0;
      for (i = 0; i < n; i++) {
        x = macc(i, j);
        
        if (x != 3) {
          xSum += x;
          xxSum += x*x;
        }
      }
      sumX[j] = xSum;
      sumXX[j] = xxSum;
    }
  }
  
  for (k = 0; k < p; k++) {
    
    j0 = ord[k] - 1; // C++ index
    
    if (remain[j0]) {
      
      remain[j0] = false;
      keep[j0] = true;
      int j_min = std::max(0, j0 - size);
      int j_max = std::min(p, j0 + size + 1);
      for (j = j_min; j < j_max; j++) {
        if (remain[j]) {
          
          N = n;
          xySum = 0;
          xSum = sumX[j];
          ySum = sumX[j0];
          xxSum = sumXX[j];
          yySum = sumXX[j0];
          
          for (i = 0; i < n; i++) {
            x = macc(i, j);
            y = macc(i, j0);
            
            if (y == 3) {
              N--;
              if (x == 3) { // both missing
                // nothing to do
              } else { // y is missing but not x
                xSum  -= x;
                xxSum -= x*x;
              }
            } else {
              if (x == 3) { // x is missing but not y
                ySum  -= y;
                yySum -= y*y;
                N--;
              } else { // both not missing
                xySum += x * y;
              }
            }
          }
          
          num = xySum - xSum * ySum / N;
          deno_x = xxSum - xSum * xSum / N;
          deno_y = yySum - ySum * ySum / N;
          r2 = num * num / (deno_x * deno_y);
          if (r2 > thr) remain[j] = false;
        }
      }
    }
  }
  
  return(keep);
}

/******************************************************************************/

// Dispatch function for clumping
// [[Rcpp::export]]
LogicalVector clumping(SEXP obj,
                       const NumericMatrix& lookup,
                       const IntegerMatrix& lookup_byte,
                       const IntegerVector& colInd,
                       const IntegerVector& ord,
                       LogicalVector& remain,
                       int size, 
                       double thr) {
  
  if (Rf_isMatrix(obj)) {
    matAcc macc(obj, lookup, colInd);
    return clumping(macc, ord, remain, size, thr);
  } else {
    XPtr<bed> xpMat(obj);
    bedAcc macc(xpMat, lookup, lookup_byte, colInd);
    return clumping(macc, ord, remain, size, thr);
  }
}

/******************************************************************************/