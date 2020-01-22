#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;

CharacterVector get_geno_char(CharacterVector allele_sep){
  int n_delim = allele_sep.size();  
  CharacterVector geno_char(n_delim * 4 + 2);
  for (int k = 0; k < n_delim; k++){
    geno_char[k * 4] = "0" + allele_sep[k] + "0";
    geno_char[k * 4 + 1] = "0" + allele_sep[k] + "1";
    geno_char[k * 4 + 2] = "1" + allele_sep[k] + "0";
    geno_char[k * 4 + 3] = "1" + allele_sep[k] + "1";
  }
  geno_char[n_delim * 4] = "0";
  geno_char[n_delim * 4 + 1] = "1";
  return(geno_char);
}

IntegerVector get_geno_int(CharacterVector allele_sep){
  int n_delim = allele_sep.size();  
  IntegerVector geno_int(n_delim * 4 + 2);
  for (int k = 0; k < n_delim; k++){
    geno_int[k * 4] = 0;
    geno_int[k * 4 + 1] = 1;
    geno_int[k * 4 + 2] = 1;
    geno_int[k * 4 + 3] = 2;
  }
  geno_int[n_delim * 4] = 0;
  geno_int[n_delim * 4 + 1] = 1;
  return(geno_int);
}

int check_line_na(CharacterVector string_geno_row, CharacterVector geno_char){
  int nIND = string_geno_row.size();
  int M = geno_char.size();
  int count = 0;
  int na_tot = 0;
  int skip_bool = 0;
  
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < M; k++){
      if (string_geno_row[j] == geno_char[k]){
        count += 1;
      }
    }
    if (count == 0){
      na_tot += 1;
    }
  }
  if (na_tot >= nIND){
    skip_bool = 1;
  }
  return(skip_bool);
}

//' Convert vcfR genotype matrices
//'
//' \code{vcf_convert} converts outputs of \code{extract.gt} to the format 
//' \code{pcadapt}.
//'
//' @param string_geno a genotype matrix extracted from a VCF file with `vcfR`. 
//' @param output a character string indicating the name of the output file.
//' @param allele.sep a vector of characters indicating what delimiters are used 
//' to separate alleles.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
IntegerVector vcf_convert(CharacterMatrix string_geno, 
                          std::string output, 
                          CharacterVector allele_sep){
  int nIND = string_geno.ncol();
  int nSNP = string_geno.nrow();
  int n_delim = allele_sep.size();
  CharacterVector geno_char(n_delim * 4 + 2);
  IntegerVector geno_int(n_delim * 4 + 2);
  geno_char = get_geno_char(allele_sep);
  geno_int = get_geno_int(allele_sep);
  CharacterVector geno_row(nIND);
  IntegerVector res(nSNP);
  
  int count_ij;
  int skip_bool;
  int na_tot = 0;
  
  FILE *file = fopen(output.c_str(), "w");
  int i = 0;
  while (i < nSNP){
    for (int j = 0; j < nIND; j++){
      geno_row[j] = string_geno(i, j);  
    }
    skip_bool = check_line_na(geno_row, geno_char);
    if (skip_bool == 0){
      for (int j = 0; j < nIND; j++){
        count_ij = 0;
        for (int k = 0; k < (n_delim * 4 + 2); k++){
          if (string_geno(i, j) == geno_char[k]){
            if (j < (nIND - 1)){
              fprintf(file, "%d ", geno_int[k]);
            } else if (j == (nIND - 1)){
              fprintf(file, "%d", geno_int[k]);
            }
            count_ij += 1;
          }
        }
        if (count_ij == 0){
          if (j < (nIND - 1)){
            fprintf(file, "%d ", 9);
          } else if (j == (nIND - 1)){
            fprintf(file, "%d", 9);
          }  
        }
      }
      fprintf(file, "\n");
      i++;
    } else {
      res[i] = 1;
      i++;
      na_tot++;
    }
  }
  fclose(file);
  return(res);
}

