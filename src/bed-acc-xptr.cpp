/******************************************************************************/

#include <pcadapt/bed-acc.h>

/******************************************************************************/

inline size_t ceil4(size_t n) {
  return (n + 3) / 4;
}

/******************************************************************************/

bed::bed(std::string path, int n, int p) : 
  n(n), p(p), n_byte(ceil4(n)) {
  
  // Rcout << this->n << " ; " << this->p << " ; " << this->n_byte << std::endl;
  
  try {
    this->file = file_mapping(path.c_str(), read_only);
  } catch(interprocess_exception& e) {
    throw std::runtime_error("File not found.");
  }
  
  this->file_region = mapped_region(this->file, read_only);
  this->file_data = 
    static_cast<const unsigned char*>(this->file_region.get_address());
  
  if (!(this->file_data[0] == '\x6C' && this->file_data[1] == '\x1B')) {
    throw std::runtime_error("File is not a binary PED file.");
  }
  
  /* Check mode: 00000001 indicates the default variant-major mode (i.e.
   list all samples for first variant, all samples for second variant,
   etc), 00000000 indicates the unsupported sample-major mode (i.e. list
   all variants for the first sample, list all variants for the second
   sample, etc */
  if (this->file_data[2] != '\x01') {
    throw std::runtime_error("Sample-major mode is not supported.");
  }
  
  // Point after this magic number
  this->file_data += 3;
  
  // Check if given dimensions match the file
  if ((3 + this->n_byte * this->p) != this->file_region.get_size()) {
    throw std::runtime_error("n or p does not match the dimensions of the file.");
  }
}

/******************************************************************************/

// [[Rcpp::export]]
SEXP bedXPtr(std::string path, int n, int p) {
  
  // http://gallery.rcpp.org/articles/intro-to-exceptions/
  try {
    /* Create a pointer to a bedAcc object and wrap it as an external pointer
     http://dirk.eddelbuettel.com/code/rcpp/Rcpp-modules.pdf */
    XPtr<bed> ptr(new bed(path, n, p), true);
    // Return the external pointer to the R side
    return ptr;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
    return 0;
  }
}

/******************************************************************************/