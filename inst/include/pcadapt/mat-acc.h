#ifndef MAT_ACC_H
#define MAT_ACC_H

/******************************************************************************/

#include <Rcpp.h>

using namespace Rcpp;
using std::size_t;

/******************************************************************************/

class matAcc {
public:
  matAcc(const IntegerMatrix& mat,
         const IntegerVector& col_ind) {
    
    n = mat.nrow();
    p = col_ind.size();
    _pMat = &(mat(0, 0));
    
    std::vector<size_t> col_ind2(p);
    for (size_t j = 0; j < p; j++) {
      // 'col_ind' indices comes from R, so begins at 1
      col_ind2[j] = static_cast<size_t>(col_ind[j] - 1);;
    }
    _col_ind = col_ind2;
  };
  
  size_t nrow() const { return n; }
  size_t ncol() const { return p; }
  
  inline int operator() (size_t i, size_t j) {
    int geno = _pMat[i + _col_ind[j] * n];
    return IntegerVector::is_na(geno) ? 3 : geno;
  }
  
protected:
  const int * _pMat;
  size_t n;
  size_t p;
  std::vector<size_t> _col_ind;
};

/******************************************************************************/

class matAccScaled : public matAcc {
public:
  matAccScaled(const IntegerMatrix& mat,
               const IntegerVector& col_ind,
               // af should be ALL allele frequencies
               const NumericVector& af,
               double ploidy,
               double NA_VAL) : matAcc(mat, col_ind) {
    
    _lookup_scale = NumericMatrix(4, p);
    for (size_t j = 0; j < p; j++) {
      double af_j = af[_col_ind[j]];
      for (size_t i = 0; i < 3; i++) {
        _lookup_scale(i, j) = 
          (i - ploidy * af_j) / sqrt(ploidy * af_j * (1 - af_j));
      }
      _lookup_scale(3, j) = NA_VAL;
    }
  };
  
  inline double operator() (size_t i, size_t j) {
    int geno = matAcc::operator()(i, j);
    return _lookup_scale(geno, j);
  }
  
protected:
  NumericMatrix _lookup_scale;
};

/******************************************************************************/

#endif // MAT_ACC_H
