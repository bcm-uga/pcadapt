#ifndef MAT_ACC_H
#define MAT_ACC_H

#include <Rcpp.h>
using namespace Rcpp;
using std::size_t;

class matAcc {
public:
  matAcc(const IntegerMatrix& mat, 
         const NumericMatrix& lookup_scale) {
    n = mat.nrow();
    p = mat.ncol();
    _lookup_scale = lookup_scale;
    _pMat = &(mat(0, 0));
  };
  
  size_t nrow() const { return n; }
  size_t ncol() const { return p; }
  
  inline double operator() (size_t i, size_t j) {
    int geno = _pMat[i + j * n];
    return IntegerVector::is_na(geno) ? 3 : _lookup_scale(geno, j);
  }
  
private:
  const int * _pMat;
  size_t n;
  size_t p;
  NumericMatrix _lookup_scale;
};

#endif // MAT_ACC_H