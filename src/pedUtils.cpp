/******************************************************************************/

#include <Rcpp.h>

using namespace Rcpp;

#define SEP " \t\n"

/******************************************************************************/

void print_error_global(const char *msg, char *file, int n){
  Rprintf("\n");
  if (!strcmp(msg, "open")){
    Rprintf("Error: unable to open file %s. Please check that the name of the file is correct.\n", file);
  } else if (!strcmp(msg, "read")){
    Rprintf("Error: unable to read file %s. Please check that the format is correct.\n", file);
  } else if (!strcmp(msg, "interne")){
    Rprintf("Error: internal error. Please run the program again.\n");
  } else if (!strcmp(msg, "constant")){
    Rprintf("Error: %d SNPs are invariant. Please remove these SNPs before running the program.\n", n);
  } else if (!strcmp(msg, "nan")){
    Rprintf("Error: internal error. Please run the program again.\n");
  } else {
    Rprintf("Error: internal error.\n");
  }
  Rprintf("\n");
  stop("File conversion aborted.");
}

FILE* fopen_read(char *file_data){
  FILE *m_File = fopen(file_data, "r");
  const char *msg = "open";
  if (!m_File){
    print_error_global(msg, file_data, 0);
  }
  return m_File;
}

FILE* fopen_write(char *file_data){
  FILE *m_File = fopen(file_data, "w");
  const char *msg = "open";
  if (!m_File){
    print_error_global(msg, file_data, 0);
  }
  return m_File;
}

void write_geno(char *output_file, int N, int M, int *data){
  FILE *file = NULL;
  file = fopen_write(output_file);
  for (int j = 0; j < M; j++){
    for (int i = 0; i < N; i++){
      if (i < N - 1){
        fprintf(file, "%d ", data[i * M + j]);
      }
      if (i == (N - 1)){
        fprintf(file, "%d", data[i * M + j]);
      }
    }
    fprintf(file, "\n");
  }
  fclose(file);
}

int nb_cols_lfmm(char *file){
  FILE *fp = fopen_read(file);
  int cols = 0;
  int c;
  char *token;
  c = fgetc(fp);
  while ((c != EOF) && (c != 10)){
    c = fgetc(fp);
    cols++;
  }
  fclose(fp);
  fp = fopen_read(file);
  char *szbuff = new char[2 * cols];
  token = fgets(szbuff, 2 * cols, fp);
  cols = 0;
  token = strtok(szbuff, SEP);
  while (token){
    cols++;
    token = strtok(NULL, SEP);
  }
  fclose(fp);
  delete[] szbuff;
  return cols;
}

int nb_lines(char *file, int M){
  FILE *fp = fopen_read(file);
  int lines = 0;
  int max_char_per_line = 20 * M + 10;
  char *szbuff = new char[max_char_per_line];
  while (fgets(szbuff, max_char_per_line, fp)){
    lines++;
  }  
  fclose(fp);
  delete[] szbuff;
  return lines;
}

void test_column(char *file, FILE *m_File, int i, int j, int N, char *token){
  if (i != N){
    Rprintf("Error: unable to read file %s. Inconsistent number of columns.\n", file);
    fclose(m_File);
    stop("File conversion aborted.");
  }
  if (token && *token != '\n' && *token != EOF){
    Rprintf("Error: unable to read file %s. Inconsistent number of columns.\n", file);
    fclose(m_File);
    stop("File conversion aborted.");
  }
}

void test_line(char *file, FILE *m_File, int i, int N){
  if (i != N){
    Rprintf("Error: unable to read file %s. Inconsistent number of lines.\n", file);
    fclose(m_File);
    stop("File conversion aborted.");
  }
  if (!feof(m_File)){
    Rprintf("Error: unable to read file %s. Inconsistent number of lines.\n", file);
    fclose(m_File);
    stop("File conversion aborted.");
  }
}

char* next_token(char* input_file, int i, int j){
  char* token = strtok(NULL, SEP);
  if (!token){
    if (j == 0){
      Rprintf("Error while reading individual information at line %d.\n", i);
    }
    else {
      Rprintf("Error while reading file %s.\n", input_file);
      stop("File conversion aborted.");
    }
  }
  return token;
}

/******************************************************************************/

void test_token_ped(char token, int j, int i, char* input_file){
  if (!(token == '0' || token == '1' || token == '2' || token == 'A' || token == 'C' || token == 'T' || token == 'G')){
    Rprintf("Error: in file %s, line %d, one allele of SNP %d is '%c' and not 0, 1, 2, A, C, T, or G.\n", input_file, i, j, token);
    stop("File conversion aborted.");
  }
}

/******************************************************************************/

void fill_line_ped(int *data, char* szbuff, int M, int i, char* input_file, FILE *File, char *ref){
  char *token1, *token2;
  int tmp;
  token1 = strtok(szbuff, SEP);
  if (!token1){
    Rprintf("Error while reading individual information at line %d.\n", i + 1);
    stop("File conversion aborted.");
  }
  for(int j = 0; j < 5; j++){
    token1 = next_token(input_file, i + 1, 0);
  }
  int j = 0;
  token1 = strtok(NULL, SEP);
  token2 = strtok(NULL, SEP);
  while (token1 && token2 && token1[0] != EOF && token2[0] != EOF && token1[0] != 10 && token2[0] != 10 && j < M){
    test_token_ped(token1[0], j + 1, i + 1, input_file);
    test_token_ped(token2[0], j + 1, i + 1, input_file);
    tmp = 0;
    if (ref[j] == '0') {
      if (token1[0] == '0' && token2[0] == '0'){
        tmp = 9;
      } else if (token2[0] != '0' && token1[0] == '0'){
        ref[j] = token2[0];
        tmp = 9;
      } else if (token1[0] != '0' && token2[0] == '0'){
        ref[j] = token1[0];
        tmp = 9;
      } else {
        ref[j] = token2[0];
        tmp = 1;
        if (token1[0] == ref[j]){
          tmp ++;
        }
      }
    } else {
      if(token1[0] == '0' || token2[0] == '0'){
        tmp = 9;		
      }
      else {
        if (token1[0] == ref[j]){
          tmp ++;
        }
        if (token2[0] == ref[j]){
          tmp ++;
        }
      }
    }
    data[i * M + j] = tmp;
    token1 = strtok(NULL, SEP);
    token2 = strtok(NULL, SEP);
    j ++;
  }
  test_column(input_file, File, j, i + 1, M, token1);
}

/******************************************************************************/

void read_ped(char* input_file, int N, int M, int* data){
  FILE *File = NULL;
  int max_char_per_line = 5 * (M + 6) + 20;
  char *szbuff = new char[max_char_per_line];
  char *ref = new char[M];
  for (int i = 0; i < M; i++){
    ref[i] = '0';
  }
  File = fopen_read(input_file);
  int i = 0;
  while (fgets(szbuff, max_char_per_line, File) && i < N){
    fill_line_ped(data, szbuff, M, i, input_file, File, ref);	
    i++;
  }
  test_line(input_file, File, i, N);
  fclose(File);
  delete[] szbuff;
  delete[] ref;
}

/******************************************************************************/

void ped2geno(char *input_file, char* output_file, int *N, int *M){
  int nb;
  nb = nb_cols_lfmm(input_file);
  *M = (nb - 6) / 2;
  *N = nb_lines(input_file, nb);
  int *data = new int[(*N)*(*M)];
  read_ped(input_file, *N, *M, data);
  write_geno(output_file, *N, *M, data);
  delete[] data;
}

/******************************************************************************/

//' Summary
//'
//' \code{print_convert} prints out a summary of the file conversion.
//'
//' @param input a genotype matrix or a character string specifying the name of the file to be converted.
//' @param output a character string specifying the name of the output file.
//' @param M an integer specifying the number of genetic markers present in the data.
//' @param N an integer specifying the number of individuals present in the data.
//' @param pool an integer specifying the type of data. `0` for genotype data, `1` for pooled data.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
void print_convert(std::string input, std::string output, int M, int N, int pool){
  char *writable_in = new char[input.size() + 1];
  std::copy(input.begin(), input.end(), writable_in);
  writable_in[input.size()] = '\0';
  char *writable_out = new char[output.size() + 1];
  std::copy(output.begin(), output.end(), writable_out);
  writable_out[output.size()] = '\0';
  if (pool == 0){
    Rprintf("Summary:\n\n");
    Rprintf("\t- input file:\t\t\t\t%s\n", writable_in);
    Rprintf("\t- output file:\t\t\t\t%s\n\n", writable_out);
    Rprintf("\t- number of individuals detected:\t%d\n", N);
    Rprintf("\t- number of loci detected:\t\t%d\n\n", M);
  } else if (pool == 1){
    Rprintf("Summary:\n\n");
    Rprintf("\t- input file:\t\t\t\t%s\n", writable_in);
    Rprintf("\t- output file:\t\t\t\t%s\n\n", writable_out);
    Rprintf("\t- number of pools detected:\t%d\n", N);
    Rprintf("\t- number of loci detected:\t\t%d\n\n", M);  
  }
}

/******************************************************************************/

//' Convert ped files
//'
//' \code{ped2pcadapt} converts \code{ped} files to the format \code{pcadapt}.
//'
//' @param input a character string specifying the name of the file to be converted.
//' @param output a character string specifying the name of the output file.
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
int ped2pcadapt(std::string input, std::string output){
  int M;
  int N;
  char *writable_in = new char[input.size() + 1];
  std::copy(input.begin(), input.end(), writable_in);
  writable_in[input.size()] = '\0';
  char *writable_out = new char[output.size() + 1];
  std::copy(output.begin(), output.end(), writable_out);
  writable_out[output.size()] = '\0';
  ped2geno(writable_in, writable_out, &N, &M);
  print_convert(input, output, M, N, 0);
  delete[] writable_in;
  delete[] writable_out;
  return(0);
}

/******************************************************************************/
