// [[Rcpp::depends(BH, pcadapt)]]
#include <pcadapt/bed-acc.h>

// [[Rcpp::export]]
IntegerVector bed_col(SEXP xptr, IntegerVector col) {
  
  XPtr<bed> xpMat(xptr);
  bedAcc macc(xpMat, col);
  
  size_t n = macc.nrow();
  IntegerVector res(n);
  
  Rcout << n << " : " << macc.ncol() << std::endl;
  
  for (size_t i = 0; i < n; i++)
    res[i] = macc(i, 0);
  
  return res;
}


/*** R
bedfile <- system.file("extdata", "geno3pops.bed", package = "pcadapt")
bed0 <- pcadapt::read.pcadapt(bedfile, type = "bed")
xptr <- pcadapt:::bedXPtr(unclass(bed0), attr(bed0, "n"), attr(bed0, "p"))
test <- bed_col(xptr, 1)
test2 <- scan(text = readLines(sub("\\.bed$$", ".pcadapt", bedfile), n = 1), what = 1L)
all.equal(test, test2)
*/
