#ifndef BED_ACC_H
#define BED_ACC_H

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/noncopyable.hpp>
#include <Rcpp.h>

using namespace boost::interprocess;
using namespace Rcpp;
using std::size_t;

class bed : private boost::noncopyable {
public:
  bed(const std::string path, int n, int p);
  
  const unsigned char* matrix() const { return file_data; }
  size_t nrow() const { return n; }
  size_t ncol() const { return p; }
  size_t nbyte() const { return n_byte; }
  
private:
  boost::interprocess::file_mapping file;
  boost::interprocess::mapped_region file_region;
  const unsigned char* file_data;
  size_t n;
  size_t p;
  size_t n_byte;
};


class bedAcc {
public:
  bedAcc(const bed * bedPtr, 
         const NumericMatrix& lookup_scale,
         const IntegerMatrix& lookup_byte) {
    
    n = bedPtr->nrow();
    p = bedPtr->ncol();
    n_byte = bedPtr->nbyte(); 
    _lookup_scale = lookup_scale;
    _lookup_byte = lookup_byte;
    _pMat = bedPtr->matrix();
  };
  
  size_t nrow() const { return n; }
  size_t ncol() const { return p; }
  
  inline double operator() (size_t i, size_t j) {
    size_t i2 = i / 4;
    size_t i3 = i % 4;
    const unsigned char byte = _pMat[i2 + j * n_byte];
    int geno = _lookup_byte(i3, byte);
    return _lookup_scale(geno, j);
  }
  
private:
  const unsigned char* _pMat;
  size_t n;
  size_t p;
  size_t n_byte;
  NumericMatrix _lookup_scale;
  IntegerMatrix _lookup_byte;
};

#endif // BED_ACC_H