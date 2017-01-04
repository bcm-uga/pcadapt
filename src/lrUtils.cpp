#include <Rcpp.h>
#include <fstream>
#include <string>
#include <sstream>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix lrfile(NumericMatrix scores, int nSNP){
  int i = 0;
  int na = 0;
  int low_AF_tot = 0;
  int nIND = scores.nrow();
  int K = scores.ncol();
  NumericMatrix z(nSNP,K);
  NumericVector y_pred(nIND);
  NumericVector residuals(nSNP);
  NumericMatrix t_scores(K,nIND);
  for (int a=0;a<K;a++){
    for (int b=0;b<nIND;b++){
      t_scores(a,b) = scores(b,a);
    }
  }
  return t_scores;
}


// [[Rcpp::export]]
CharacterVector read_file_cpp2(std::string path) {
  std::ifstream in(path.c_str());
  std::string contents;
  in.seekg(0, std::ios::end);
  contents.resize(in.tellg());
  in.seekg(0, std::ios::beg);
  in.read(&contents[0], contents.size());
  in.close();
  return(contents);
}
