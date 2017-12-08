/******************************************************************************/

#include <mmapcharr/charsep-acc.h>
#include <fstream>

using namespace Rcpp;
using namespace std;

/******************************************************************************/

template <class C>
void writebed(const char * filename,
              C macc,
              const RawVector& tab) {
  
  int n = macc.nrow();
  int m = macc.ncol();
  int length = ceil((double)n / 4); // DO NOT USE INTEGERS WITH CEIL
  
  char *buffer = new char[std::max(3, length)];
  ofstream myFile(filename, ios::out | ios::binary);
  
  // magic number
  buffer[0] = 108; buffer[1] = 27; buffer[2] = 1;
  myFile.write(buffer, 3);
  
  int i, j, k, ind, coef;
  
  for (j = 0; j < m; j++) {
    k = 0;
    for (i = 0; i <= n-4; i += 4) {
      ind = (macc(i, j) + 4 * macc(i+1, j)) +
        (16 * macc(i+2, j) + 64 * macc(i+3, j));
      buffer[k++] = tab[ind];
    }
    ind = 0; coef = 1;
    for (; i < n; i++) {
      ind += coef * macc(i, j);
      coef *= 4;
    }
    buffer[k] = tab[ind];
    myFile.write(buffer, length); // faster to use (char*)? Nop
  }
  
  myFile.close();
  delete[] buffer;
}

/******************************************************************************/

// Dispatch function for writebed
// [[Rcpp::export]]
void writebed(const char * filename,
              Environment e,
              const RawVector& tab, 
              bool is_pcadapt) {
  
  XPtr<charSep> xpMat = e["address"];
  
  if (is_pcadapt) {
    charSepAccTranspose<int, INTSXP> macc(xpMat, e["code"]);
    writebed(filename, macc, tab);
  } else {  // lfmm
    charSepAcc<int, INTSXP> macc(xpMat, e["code"]);
    writebed(filename, macc, tab);
  }
}

/******************************************************************************/