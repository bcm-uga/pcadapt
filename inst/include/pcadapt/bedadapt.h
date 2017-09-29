#ifndef bedadapt_H
#define bedadapt_H

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <Rcpp.h>

using namespace Rcpp;

class bedadapt {
public:
  bedadapt(std::string path, std::size_t n, std::size_t p);
  size_t nrow() const { return n; }
  size_t ncol() const { return p; }
  NumericVector AF();
  IntegerVector extractSNP(int j);
  NumericVector prodMatVec(const NumericVector &x, 
                           const NumericVector &m, 
                           const NumericVector &s);
  NumericVector prodtMatVec(const NumericVector &x, 
                            const NumericVector &m, 
                            const NumericVector &s);
  NumericMatrix linReg(const NumericMatrix &u,
                       const NumericVector &d,
                       const NumericMatrix &v,
                       const NumericVector &m);
private:
  boost::interprocess::file_mapping file;
  boost::interprocess::mapped_region file_region;
  const char* file_data;
  std::size_t n;
  std::size_t p;
  unsigned short int byte_padding; // Each new "row" starts a new byte
  static const unsigned short int length_header;
  char get_byte(std::size_t i);
};

#endif