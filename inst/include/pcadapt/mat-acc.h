#ifndef MAT_ACC_H
#define MAT_ACC_H

#include <Rcpp.h>
using namespace Rcpp;
using std::size_t;

class matAcc {
public:
  matAcc(const IntegerMatrix& mat,  // lookup must be corresponding to 'col_ind'
         const NumericMatrix& lookup_scale,
         const IntegerVector& col_ind) {
    n = mat.nrow();
    p = col_ind.size();
    _lookup_scale = lookup_scale;
    _pMat = &(mat(0, 0));
    _val_na = lookup_scale(3, 0);
    
    std::vector<size_t> col_ind2(p);
    for (size_t j = 0; j < p; j++) {
      // 'col_ind' indices comes from R, so begins at 1
      col_ind2[j] = static_cast<size_t>(col_ind[j] - 1);;
    }
    _col_ind = col_ind2;
  };
  
  size_t nrow() const { return n; }
  size_t ncol() const { return p; }
  
  inline double operator() (size_t i, size_t j) {
    int geno = _pMat[i + _col_ind[j] * n];
    return IntegerVector::is_na(geno) ? _val_na : _lookup_scale(geno, j);
  }
  
private:
  const int * _pMat;
  size_t n;
  size_t p;
  NumericMatrix _lookup_scale;
  double _val_na;
  std::vector<size_t> _col_ind;
};

#endif // MAT_ACC_H