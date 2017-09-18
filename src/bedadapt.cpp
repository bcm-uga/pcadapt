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
 
 /******************************************************************************/
 
 NumericVector bedadapt::minor_AF() {
   Rcpp::NumericVector maf(this->p);
   for (int j = 0; j < this->p; j++) {
     int n_available = 0; // Counts the number of available values for SNP j
     // No need to convert from 1-index to 0-index
     size_t SNP_start = (this->n + this->byte_padding) * j;
     size_t SNP_end = (this->n + this->byte_padding) * (j + 1);
     
     // Every byte encodes 4 genotypes, find the one of interest
     size_t byte_start = std::floor(SNP_start / 4);
     size_t byte_end = std::floor(SNP_end / 4);
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
 
 /******************************************************************************/
 
 // [[Rcpp::export]]
 RObject cmpt_minor_af_BED(RObject xp_) {
   XPtr<bedadapt> ptr(xp_);
   NumericVector res = ptr->minor_AF();
   return res;
 }
 
 // [[Rcpp::export]]
 RObject prodMatVec_export(RObject xp_, NumericVector x) {
   // Convert inputs to appropriate C++ types
   XPtr<bedadapt> ptr(xp_);
   NumericVector res = ptr->prodMatVec(x);
   return res;
 }
 
 NumericVector bedadapt::prodMatVec(NumericVector x) {
   NumericVector y(this->n);
   for (int j = 0; j < this->p; j++) {
     // No need to convert from 1-index to 0-index
     size_t SNP_start = (this->n + this->byte_padding) * j;
     size_t SNP_end = (this->n + this->byte_padding) * (j + 1);
     
     // Every byte encodes 4 genotypes, find the one of interest
     size_t byte_start = std::floor(SNP_start / 4);
     size_t byte_end = std::floor(SNP_end / 4);
     int acc = 0;
     char raw_element;
     char genotype;
     for (size_t byte = byte_start; byte < byte_end; byte++) {
       raw_element = get_byte(byte);
       for (int idx = 0; idx < 4; idx++) {
         genotype = raw_element >> (idx * 2) & 3;
         if (genotype == 0 && acc < this->n) {
           y[acc] += (2.0 * x[j]); // homozygous AA
         } else if (genotype == 2 && acc < this->n) {
           y[acc] += (1.0 * x[j]); // heterozygous AB
         }
         acc ++;
       }
     }
   }
   return y;
 }
 
 // [[Rcpp::export]]
 RObject prodtMatVec_export(RObject xp_, NumericVector x) {
   // Convert inputs to appropriate C++ types
   XPtr<bedadapt> ptr(xp_);
   NumericVector res = ptr->prodtMatVec(x);
   return res;
 }
 
 NumericVector bedadapt::prodtMatVec(NumericVector x) {
   Rcpp::NumericVector y(this->p);
   for (int j = 0; j < this->p; j++) {
     // No need to convert from 1-index to 0-index
     std::size_t SNP_start = (this->n + this->byte_padding) * j;
     std::size_t SNP_end = (this->n + this->byte_padding) * (j + 1);
     
     // Every byte encodes 4 genotypes, find the one of interest
     std::size_t byte_start = std::floor(SNP_start / 4);
     std::size_t byte_end = std::floor(SNP_end / 4);
     int acc = 0;
     char raw_element;
     char genotype;
     for (size_t byte = byte_start; byte < byte_end; byte++) {
       raw_element = get_byte(byte);
       for (int idx = 0; idx < 4; idx++) {
         genotype = raw_element >> (idx * 2) & 3;
         if (genotype == 0 && acc < this->n) {
           y[j] += 2.0 * x[acc]; // homozygous AA
         } else if (genotype == 2 && acc < this->n) {
           y[j] += 1.0 * x[acc]; // heterozygous AB
         }
         acc ++;
       }
     }
   }
   return y;
 }