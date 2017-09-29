/******************************************************************************/
/* Inspired from
* https://github.com/QuantGen/BEDMatrix/blob/master/src/BEDMatrix.cpp **********
* https://github.com/privefl/bigstatsr/blob/master/src/FBM.cpp ****************/

#include <pcadapt/bedadapt.h>

using namespace Rcpp;
using namespace boost::interprocess;
using std::size_t;

/******************************************************************************/

bedadapt::bedadapt(std::string path, size_t n, size_t p) : n(n), p(p), byte_padding((n % 4 == 0) ? 0 : 4 - (n % 4)) {
  
  try {
    this->file = file_mapping(path.c_str(), read_only);
  } catch(interprocess_exception& e) {
    throw std::runtime_error("File not found.");
  }
  
  this->file_region = mapped_region(this->file, read_only);
  this->file_data = static_cast<const char*>(this->file_region.get_address());
  
  if (!(this->file_data[0] == '\x6C' && this->file_data[1] == '\x1B')) {
    throw std::runtime_error("File is not a binary PED file.");
  }
  
  // Check mode: 00000001 indicates the default variant-major mode (i.e.
  // list all samples for first variant, all samples for second variant,
  // etc), 00000000 indicates the unsupported sample-major mode (i.e. list
  // all variants for the first sample, list all variants for the second
  // sample, etc)
  if (this->file_data[2] != '\x01') {
    throw std::runtime_error("Sample-major mode is not supported.");
  }
  
  // Get number of bytes
  const size_t num_bytes = this->file_region.get_size();
  
  // Check if given dimensions match the file
  if ((this->n * this->p) + (this->byte_padding * this->p) != (num_bytes - this->length_header) * 4) {
    throw std::runtime_error("n or p does not match the dimensions of the file.");
  }
}

const unsigned short int bedadapt::length_header = 3;

/******************************************************************************/

// [[Rcpp::export]]
SEXP bedadaptXPtr(std::string path, int n, int p) {
  
  // http://gallery.rcpp.org/articles/intro-to-exceptions/
  try {
    // Create a pointer to a bedadapt object and wrap it as an external pointer
    // http://dirk.eddelbuettel.com/code/rcpp/Rcpp-modules.pdf
    XPtr<bedadapt> ptr(new bedadapt(path, n, p), true);
    // Return the external pointer to the R side
    return ptr;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
    return 0;
  }
}

/******************************************************************************/

char bedadapt::get_byte(std::size_t i) {
  char genotypes = this->file_data[i + this->length_header];
  return genotypes;
}

std::size_t byte_position(int npbp, int j) {
  // Every byte encodes 4 genotypes, find the one of interest
  std::size_t res = std::floor(npbp * j / 4);
  return(res);
}

IntegerVector bedadapt::extractSNP(int j) {
  IntegerVector snp(this->n);
  size_t byte_start = byte_position(this->n + this->byte_padding, j);
  size_t byte_end = byte_position(this->n + this->byte_padding, j + 1);  
  char raw_element;
  char genotype;
  int i = 0;
  for (size_t byte = byte_start; byte < byte_end; byte++) {
    for (size_t byte = byte_start; byte < byte_end; byte++) {
      raw_element = get_byte(byte);
      for (int idx = 0; idx < 4; idx++) {
        genotype = raw_element >> (idx * 2) & 3;
        if (i < this->n) {
          if (genotype == 0) {
            snp[i] = 2; // homozygous AA
          } else if (genotype == 2) {
            snp[i] = 1; // heterozygous AB
          } else if (genotype == 3) {
            snp[i] = 0; // homozygous BB
          } else {
            snp[i] = NA_INTEGER;
          }
          i++;
        }
      }
    }  
  }
  return(snp);
}

/******************************************************************************/

NumericVector bedadapt::AF() {
  Rcpp::NumericVector maf(this->p);
  for (int j = 0; j < this->p; j++) {
    int n_available = 0; // Counts the number of available values for SNP j
    size_t byte_start = byte_position(this->n + this->byte_padding, j);
    size_t byte_end = byte_position(this->n + this->byte_padding, j + 1);
    int acc = 0;
    char raw_element;
    char genotype;
    for (size_t byte = byte_start; byte < byte_end; byte++) {
      raw_element = get_byte(byte);
      for (int idx = 0; idx < 4; idx++) {
        genotype = raw_element >> (idx * 2) & 3;
        if (genotype == 0 && acc < this->n) {
          maf[j] += 2.0; // homozygous AA
          n_available++;
        } else if (genotype == 2 && acc < this->n) {
          maf[j] += 1.0; // heterozygous AB
          n_available++;
        } else if (genotype == 3 && acc < this->n) {
          n_available++;
        }
        acc++;
      }
    }
    maf[j] = (double) maf[j] / (2 * n_available);
  }
  return maf;
}

// [[Rcpp::export]]
RObject cmpt_af(RObject xp_) {
  XPtr<bedadapt> ptr(xp_);
  NumericVector res = ptr->AF();
  return res;
}

/******************************************************************************/

// [[Rcpp::export]]
RObject prodMatVec(RObject xp_, 
                   const NumericVector &x, 
                   const NumericVector &m, 
                   const NumericVector &s) {
  // Convert inputs to appropriate C++ types
  XPtr<bedadapt> ptr(xp_);
  NumericVector res = ptr->prodMatVec(x, m, s);
  return res;
}

NumericVector bedadapt::prodMatVec(const NumericVector &x, 
                                   const NumericVector &m, 
                                   const NumericVector &s) {
  // Input vector of length p
  // Output vector of length n
  NumericVector y(this->n);
  for (int j = 0; j < this->p; j++) {
    size_t byte_start = byte_position(this->n + this->byte_padding, j);
    size_t byte_end = byte_position(this->n + this->byte_padding, j + 1);
    int acc = 0;
    char raw_element;
    char genotype;
    for (size_t byte = byte_start; byte < byte_end; byte++) {
      raw_element = get_byte(byte);
      for (int idx = 0; idx < 4; idx++) {
        genotype = raw_element >> (idx * 2) & 3;
        if (genotype == 0 && acc < this->n) {
          y[acc] += ((2.0 - 2 * m[j]) * x[j]) / s[j]; // homozygous AA
        } else if (genotype == 2 && acc < this->n) {
          y[acc] += ((1.0 - 2 * m[j]) * x[j]) / s[j]; // heterozygous AB
        } else if (genotype == 3 && acc < this->n) {
          y[acc] -= 2 * m[j] * x[j] / s[j]; // heterozygous BB
        } 
        acc ++;
      }
    }
  }
  return y;
}

// [[Rcpp::export]]
RObject prodtMatVec(RObject xp_, 
                    const NumericVector &x, 
                    const NumericVector &m, 
                    const NumericVector &s) {
  // Convert inputs to appropriate C++ types
  XPtr<bedadapt> ptr(xp_);
  NumericVector res = ptr->prodtMatVec(x, m, s);
  return res;
}

NumericVector bedadapt::prodtMatVec(const NumericVector &x, 
                                    const NumericVector &m, 
                                    const NumericVector &s) {
  // Input vector of length n
  // Output vector of length p
  Rcpp::NumericVector y(this->p);
  for (int j = 0; j < this->p; j++) {
    size_t byte_start = byte_position(this->n + this->byte_padding, j);
    size_t byte_end = byte_position(this->n + this->byte_padding, j + 1);
    int acc = 0;
    char raw_element;
    char genotype;
    for (size_t byte = byte_start; byte < byte_end; byte++) {
      raw_element = get_byte(byte);
      for (int idx = 0; idx < 4; idx++) {
        genotype = raw_element >> (idx * 2) & 3;
        if (genotype == 0 && acc < this->n) {
          y[j] += ((2.0 - 2 * m[j]) * x[acc]) / s[j]; // homozygous AA
        } else if (genotype == 2 && acc < this->n) {
          y[j] += ((1.0 - 2 * m[j]) * x[acc]) / s[j]; // heterozygous AB
        } else if (genotype == 3 && acc < this->n) {
          y[j] -= 2 * m[j] * x[acc] / s[j]; // heterozygous BB
        }
        acc ++;
      }
    }
  }
  return y;
}

/******************************************************************************/

NumericMatrix bedadapt::linReg(const NumericMatrix &u,
                               const NumericVector &d,
                               const NumericMatrix &v,
                               const NumericVector &m) {
  int K = u.ncol();
  NumericMatrix Z(this->p, K); // z-scores
  IntegerVector G(this->n);

  for (int j = 0; j < this->p; j++) {
    /* Compute the residuals */
    G = this->extractSNP(j);
    double residual = 0;
    int n_available = 0;
    NumericVector sum_squared_u(K);   
    for (int i = 0; i < this->n; i++) {
      double y = 0;
      for (int k = 0; k < K; k++) {
        y += u(i, k) * d[k] * v(j, k) ; // Y = UDV
      }
      if (!IntegerVector::is_na(G[i])) {
        double tmp = y - (G[i] - 2 * m[j]) / sqrt(2 * m[j] * (1 - m[j]));
        residual += tmp * tmp;
        n_available++;
        for (int k = 0; k < K; k++) {
          // sum_squared_u = 1 if SNP j has no missing value
          // sum_squared_u < 1 if SNP j has at least one missing value
          sum_squared_u[k] += u(i, k) * u(i, k); 
        }
      } 
    }
    
    /* t-score */
    for (int k = 0; k < K; k++) {
      if (residual > 0 && n_available > K) {
        Z(j, k) = v(j, k) * d[k] / sqrt(residual / (n_available - K));
      }
      if (sum_squared_u[k] > 0) {
        // this should never happen
        Z(j, k) /= sqrt(sum_squared_u[k]);
      }
    }
  }
  return(Z);
}

// [[Rcpp::export]]
RObject linReg(RObject xp_, 
               const NumericMatrix &u, 
               const NumericVector &d, 
               const NumericMatrix &v, 
               const NumericVector &m) {
  // Convert inputs to appropriate C++ types
  XPtr<bedadapt> ptr(xp_);
  NumericVector res = ptr->linReg(u, d, v, m);
  return res;
}

