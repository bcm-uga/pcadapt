#ifndef BED_ACC_H
#define BED_ACC_H

/******************************************************************************/

#include <mio/mmap.hpp>
#include <system_error> // for std::error_code
#include <Rcpp.h>

using namespace Rcpp;
using std::size_t;

/******************************************************************************/

class bed {
public:
  bed(const std::string path, int n, int p);
  
  IntegerMatrix get_code(int NA_VAL = 3) const {
    
    IntegerVector num = IntegerVector::create(0, NA_VAL, 1, 2);
    IntegerMatrix code(4, 256);
    
    int i, k, k2;
    int coeff = 1;
    for (i = 0; i < 4; i++) {
      for (k = 0; k < 256; k++) {
        k2 = k / coeff;
        code(i, k) = num[k2 % 4];
      }
      coeff *= 4;
    }
    
    return code;
  }
  
  const unsigned char* matrix() const { return ro_ummap.data(); }
  
  size_t nrow()  const { return n; }
  size_t ncol()  const { return p; }
  size_t nbyte() const { return n_byte; }
  
private:
  mio::ummap_source ro_ummap;
  size_t n, p, n_byte;
};

/******************************************************************************/

class bedAcc {
public:
  bedAcc(const bed * bedPtr,
         const IntegerVector& col_ind,
         int NA_VAL = 3) {
    
    n = bedPtr->nrow();
    p = col_ind.size();
    n_byte = bedPtr->nbyte(); 
    _pMat = bedPtr->matrix();
    
    _lookup_byte = bedPtr->get_code(NA_VAL);
    
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
    const unsigned char byte = _pMat[i / 4 + _col_ind[j] * n_byte];
    return _lookup_byte(i % 4, byte);
  }
  
protected:
  const unsigned char* _pMat;
  size_t n;
  size_t p;
  size_t n_byte;
  IntegerMatrix _lookup_byte;
  std::vector<size_t> _col_ind;
};

/******************************************************************************/

class bedAccScaled : public bedAcc {
public:
  bedAccScaled(const bed * bedPtr,
               const IntegerVector& col_ind,
               // af should be ALL allele frequencies
               const NumericVector& af,
               double ploidy,
               double NA_VAL) : bedAcc(bedPtr, col_ind) {
    
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
    int geno = bedAcc::operator()(i, j);
    return _lookup_scale(geno, j);
  }
  
protected:
  NumericMatrix _lookup_scale;
};

/******************************************************************************/

#endif // BED_ACC_H
