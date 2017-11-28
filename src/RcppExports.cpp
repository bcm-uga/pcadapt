// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// af
NumericVector af(SEXP obj, const NumericMatrix& lookup_scale, const IntegerMatrix& lookup_byte);
RcppExport SEXP _pcadapt_af(SEXP objSEXP, SEXP lookup_scaleSEXP, SEXP lookup_byteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type lookup_scale(lookup_scaleSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type lookup_byte(lookup_byteSEXP);
    rcpp_result_gen = Rcpp::wrap(af(obj, lookup_scale, lookup_byte));
    return rcpp_result_gen;
END_RCPP
}
// bedXPtr
SEXP bedXPtr(std::string path, int n, int p);
RcppExport SEXP _pcadapt_bedXPtr(SEXP pathSEXP, SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(bedXPtr(path, n, p));
    return rcpp_result_gen;
END_RCPP
}
// bedadaptXPtr
SEXP bedadaptXPtr(std::string path, int n, int p);
RcppExport SEXP _pcadapt_bedadaptXPtr(SEXP pathSEXP, SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(bedadaptXPtr(path, n, p));
    return rcpp_result_gen;
END_RCPP
}
// cmpt_af
RObject cmpt_af(RObject xp_);
RcppExport SEXP _pcadapt_cmpt_af(SEXP xp_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type xp_(xp_SEXP);
    rcpp_result_gen = Rcpp::wrap(cmpt_af(xp_));
    return rcpp_result_gen;
END_RCPP
}
// prodMatVec
RObject prodMatVec(RObject xp_, const NumericVector& x, const NumericVector& m, const NumericVector& s, const LogicalVector& pass);
RcppExport SEXP _pcadapt_prodMatVec(SEXP xp_SEXP, SEXP xSEXP, SEXP mSEXP, SEXP sSEXP, SEXP passSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type xp_(xp_SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type pass(passSEXP);
    rcpp_result_gen = Rcpp::wrap(prodMatVec(xp_, x, m, s, pass));
    return rcpp_result_gen;
END_RCPP
}
// prodtMatVec
RObject prodtMatVec(RObject xp_, const NumericVector& x, const NumericVector& m, const NumericVector& s, const LogicalVector& pass);
RcppExport SEXP _pcadapt_prodtMatVec(SEXP xp_SEXP, SEXP xSEXP, SEXP mSEXP, SEXP sSEXP, SEXP passSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type xp_(xp_SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type pass(passSEXP);
    rcpp_result_gen = Rcpp::wrap(prodtMatVec(xp_, x, m, s, pass));
    return rcpp_result_gen;
END_RCPP
}
// linReg
RObject linReg(RObject xp_, const NumericMatrix& u, const NumericVector& d, const NumericMatrix& v, const NumericVector& m);
RcppExport SEXP _pcadapt_linReg(SEXP xp_SEXP, SEXP uSEXP, SEXP dSEXP, SEXP vSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type xp_(xp_SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(linReg(xp_, u, d, v, m));
    return rcpp_result_gen;
END_RCPP
}
// get_sumX
NumericVector get_sumX(SEXP obj, const NumericMatrix& lookup_scale, const IntegerMatrix& lookup_byte);
RcppExport SEXP _pcadapt_get_sumX(SEXP objSEXP, SEXP lookup_scaleSEXP, SEXP lookup_byteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type lookup_scale(lookup_scaleSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type lookup_byte(lookup_byteSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sumX(obj, lookup_scale, lookup_byte));
    return rcpp_result_gen;
END_RCPP
}
// get_denoX
NumericVector get_denoX(SEXP obj, const NumericMatrix& lookup_scale, const IntegerMatrix& lookup_byte, const NumericVector& m);
RcppExport SEXP _pcadapt_get_denoX(SEXP objSEXP, SEXP lookup_scaleSEXP, SEXP lookup_byteSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type lookup_scale(lookup_scaleSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type lookup_byte(lookup_byteSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(get_denoX(obj, lookup_scale, lookup_byte, m));
    return rcpp_result_gen;
END_RCPP
}
// clumping
LogicalVector clumping(SEXP obj, const NumericMatrix& lookup, const IntegerMatrix& lookup_byte, const IntegerVector& ord, LogicalVector& remain, const NumericVector& sumX, const NumericVector& denoX, int size, double thr);
RcppExport SEXP _pcadapt_clumping(SEXP objSEXP, SEXP lookupSEXP, SEXP lookup_byteSEXP, SEXP ordSEXP, SEXP remainSEXP, SEXP sumXSEXP, SEXP denoXSEXP, SEXP sizeSEXP, SEXP thrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type lookup(lookupSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type lookup_byte(lookup_byteSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ord(ordSEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type remain(remainSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sumX(sumXSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type denoX(denoXSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type thr(thrSEXP);
    rcpp_result_gen = Rcpp::wrap(clumping(obj, lookup, lookup_byte, ord, remain, sumX, denoX, size, thr));
    return rcpp_result_gen;
END_RCPP
}
// cmpt_cov_cpp
Rcpp::List cmpt_cov_cpp(std::string filename, arma::mat& xmatrix, double min_maf, int ploidy, int type, int blocksize);
RcppExport SEXP _pcadapt_cmpt_cov_cpp(SEXP filenameSEXP, SEXP xmatrixSEXP, SEXP min_mafSEXP, SEXP ploidySEXP, SEXP typeSEXP, SEXP blocksizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xmatrix(xmatrixSEXP);
    Rcpp::traits::input_parameter< double >::type min_maf(min_mafSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type blocksize(blocksizeSEXP);
    rcpp_result_gen = Rcpp::wrap(cmpt_cov_cpp(filename, xmatrix, min_maf, ploidy, type, blocksize));
    return rcpp_result_gen;
END_RCPP
}
// cart2bary_cpp
arma::mat cart2bary_cpp(arma::mat& X, arma::mat& P);
RcppExport SEXP _pcadapt_cart2bary_cpp(SEXP XSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(cart2bary_cpp(X, P));
    return rcpp_result_gen;
END_RCPP
}
// impute_geno
Rcpp::List impute_geno(const arma::mat& x);
RcppExport SEXP _pcadapt_impute_geno(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(impute_geno(x));
    return rcpp_result_gen;
END_RCPP
}
// impute_geno_pop
Rcpp::List impute_geno_pop(const arma::mat& x, const arma::vec& lab, const arma::vec& pop);
RcppExport SEXP _pcadapt_impute_geno_pop(SEXP xSEXP, SEXP labSEXP, SEXP popSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lab(labSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pop(popSEXP);
    rcpp_result_gen = Rcpp::wrap(impute_geno_pop(x, lab, pop));
    return rcpp_result_gen;
END_RCPP
}
// get_window
IntegerVector get_window(int i, const arma::vec& map, const double window_size);
RcppExport SEXP _pcadapt_get_window(SEXP iSEXP, SEXP mapSEXP, SEXP window_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const double >::type window_size(window_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(get_window(i, map, window_size));
    return rcpp_result_gen;
END_RCPP
}
// cmpt_global_pca
arma::mat cmpt_global_pca(const arma::mat& geno, const arma::mat& V, const arma::vec& sigma);
RcppExport SEXP _pcadapt_cmpt_global_pca(SEXP genoSEXP, SEXP VSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(cmpt_global_pca(geno, V, sigma));
    return rcpp_result_gen;
END_RCPP
}
// cmpt_local_pca
arma::mat cmpt_local_pca(const arma::mat& geno, const arma::mat& V, const arma::vec& sigma, const int beg, const int end);
RcppExport SEXP _pcadapt_cmpt_local_pca(SEXP genoSEXP, SEXP VSEXP, SEXP sigmaSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< const int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(cmpt_local_pca(geno, V, sigma, beg, end));
    return rcpp_result_gen;
END_RCPP
}
// updt_local_scores
void updt_local_scores(arma::mat& u, const arma::mat& geno, const arma::mat& V, const arma::vec& sigma, const int beg_old, const int end_old, const int beg_new, const int end_new);
RcppExport SEXP _pcadapt_updt_local_scores(SEXP uSEXP, SEXP genoSEXP, SEXP VSEXP, SEXP sigmaSEXP, SEXP beg_oldSEXP, SEXP end_oldSEXP, SEXP beg_newSEXP, SEXP end_newSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const int >::type beg_old(beg_oldSEXP);
    Rcpp::traits::input_parameter< const int >::type end_old(end_oldSEXP);
    Rcpp::traits::input_parameter< const int >::type beg_new(beg_newSEXP);
    Rcpp::traits::input_parameter< const int >::type end_new(end_newSEXP);
    updt_local_scores(u, geno, V, sigma, beg_old, end_old, beg_new, end_new);
    return R_NilValue;
END_RCPP
}
// multLinReg
NumericMatrix multLinReg(SEXP obj, const NumericMatrix& lookup_scale, const IntegerMatrix& lookup_byte, const NumericMatrix& u, const NumericVector& d, const NumericMatrix& v);
RcppExport SEXP _pcadapt_multLinReg(SEXP objSEXP, SEXP lookup_scaleSEXP, SEXP lookup_byteSEXP, SEXP uSEXP, SEXP dSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type lookup_scale(lookup_scaleSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type lookup_byte(lookup_byteSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(multLinReg(obj, lookup_scale, lookup_byte, u, d, v));
    return rcpp_result_gen;
END_RCPP
}
// cmpt_af_matrix
NumericVector cmpt_af_matrix(const NumericMatrix& G);
RcppExport SEXP _pcadapt_cmpt_af_matrix(SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(cmpt_af_matrix(G));
    return rcpp_result_gen;
END_RCPP
}
// prodGx_matrix
NumericVector prodGx_matrix(const NumericMatrix& G, const NumericVector& x, const NumericVector& m, const NumericVector& s, const LogicalVector& pass);
RcppExport SEXP _pcadapt_prodGx_matrix(SEXP GSEXP, SEXP xSEXP, SEXP mSEXP, SEXP sSEXP, SEXP passSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type pass(passSEXP);
    rcpp_result_gen = Rcpp::wrap(prodGx_matrix(G, x, m, s, pass));
    return rcpp_result_gen;
END_RCPP
}
// prodtGx_matrix
NumericVector prodtGx_matrix(const NumericMatrix& G, const NumericVector& x, const NumericVector& m, const NumericVector& s, const LogicalVector& pass);
RcppExport SEXP _pcadapt_prodtGx_matrix(SEXP GSEXP, SEXP xSEXP, SEXP mSEXP, SEXP sSEXP, SEXP passSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type pass(passSEXP);
    rcpp_result_gen = Rcpp::wrap(prodtGx_matrix(G, x, m, s, pass));
    return rcpp_result_gen;
END_RCPP
}
// nb_nona
ListOf<NumericVector> nb_nona(SEXP obj, const NumericMatrix& lookup_scale, const IntegerMatrix& lookup_byte);
RcppExport SEXP _pcadapt_nb_nona(SEXP objSEXP, SEXP lookup_scaleSEXP, SEXP lookup_byteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type lookup_scale(lookup_scaleSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type lookup_byte(lookup_byteSEXP);
    rcpp_result_gen = Rcpp::wrap(nb_nona(obj, lookup_scale, lookup_byte));
    return rcpp_result_gen;
END_RCPP
}
// covRob_cpp
Rcpp::List covRob_cpp(arma::mat& x);
RcppExport SEXP _pcadapt_covRob_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(covRob_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// print_convert
void print_convert(std::string input, std::string output, int M, int N, int pool);
RcppExport SEXP _pcadapt_print_convert(SEXP inputSEXP, SEXP outputSEXP, SEXP MSEXP, SEXP NSEXP, SEXP poolSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type input(inputSEXP);
    Rcpp::traits::input_parameter< std::string >::type output(outputSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type pool(poolSEXP);
    print_convert(input, output, M, N, pool);
    return R_NilValue;
END_RCPP
}
// ped2pcadapt
int ped2pcadapt(std::string input, std::string output);
RcppExport SEXP _pcadapt_ped2pcadapt(SEXP inputSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type input(inputSEXP);
    Rcpp::traits::input_parameter< std::string >::type output(outputSEXP);
    rcpp_result_gen = Rcpp::wrap(ped2pcadapt(input, output));
    return rcpp_result_gen;
END_RCPP
}
// lfmm2pcadapt
int lfmm2pcadapt(std::string input, std::string output);
RcppExport SEXP _pcadapt_lfmm2pcadapt(SEXP inputSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type input(inputSEXP);
    Rcpp::traits::input_parameter< std::string >::type output(outputSEXP);
    rcpp_result_gen = Rcpp::wrap(lfmm2pcadapt(input, output));
    return rcpp_result_gen;
END_RCPP
}
// sample_geno_matrix
NumericMatrix sample_geno_matrix(NumericMatrix freq, double ploidy, IntegerVector sample_size);
RcppExport SEXP _pcadapt_sample_geno_matrix(SEXP freqSEXP, SEXP ploidySEXP, SEXP sample_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< double >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sample_size(sample_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_geno_matrix(freq, ploidy, sample_size));
    return rcpp_result_gen;
END_RCPP
}
// pMatVec4
NumericVector pMatVec4(SEXP obj, const NumericVector& x, const NumericMatrix& lookup_scale, const IntegerMatrix& lookup_byte, const IntegerVector& pass);
RcppExport SEXP _pcadapt_pMatVec4(SEXP objSEXP, SEXP xSEXP, SEXP lookup_scaleSEXP, SEXP lookup_byteSEXP, SEXP passSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type lookup_scale(lookup_scaleSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type lookup_byte(lookup_byteSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type pass(passSEXP);
    rcpp_result_gen = Rcpp::wrap(pMatVec4(obj, x, lookup_scale, lookup_byte, pass));
    return rcpp_result_gen;
END_RCPP
}
// cpMatVec4
NumericVector cpMatVec4(SEXP obj, const NumericVector& x, const NumericMatrix& lookup_scale, const IntegerMatrix& lookup_byte, const IntegerVector& pass);
RcppExport SEXP _pcadapt_cpMatVec4(SEXP objSEXP, SEXP xSEXP, SEXP lookup_scaleSEXP, SEXP lookup_byteSEXP, SEXP passSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type obj(objSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type lookup_scale(lookup_scaleSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type lookup_byte(lookup_byteSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type pass(passSEXP);
    rcpp_result_gen = Rcpp::wrap(cpMatVec4(obj, x, lookup_scale, lookup_byte, pass));
    return rcpp_result_gen;
END_RCPP
}
// get_pop_size
IntegerVector get_pop_size(const StringVector& pop, const StringVector& popUnique);
RcppExport SEXP _pcadapt_get_pop_size(SEXP popSEXP, SEXP popUniqueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< const StringVector& >::type popUnique(popUniqueSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pop_size(pop, popUnique));
    return rcpp_result_gen;
END_RCPP
}
// slidingWindows_fast
arma::mat slidingWindows_fast(const arma::mat& sgeno, const arma::vec& d, const arma::mat& v, const StringVector& pop, const StringVector& popUnique, const CharacterVector& admixed, const int window_size, const arma::vec map, const int with_map);
RcppExport SEXP _pcadapt_slidingWindows_fast(SEXP sgenoSEXP, SEXP dSEXP, SEXP vSEXP, SEXP popSEXP, SEXP popUniqueSEXP, SEXP admixedSEXP, SEXP window_sizeSEXP, SEXP mapSEXP, SEXP with_mapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type sgeno(sgenoSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const StringVector& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< const StringVector& >::type popUnique(popUniqueSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type admixed(admixedSEXP);
    Rcpp::traits::input_parameter< const int >::type window_size(window_sizeSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const int >::type with_map(with_mapSEXP);
    rcpp_result_gen = Rcpp::wrap(slidingWindows_fast(sgeno, d, v, pop, popUnique, admixed, window_size, map, with_map));
    return rcpp_result_gen;
END_RCPP
}
// slidingWindows_new
arma::mat slidingWindows_new(const arma::mat& sgeno, const arma::vec& d, const arma::mat& v, const StringVector& pop, const StringVector& popUnique, const CharacterVector& admixed, const int window_size, const arma::vec map, const int with_map);
RcppExport SEXP _pcadapt_slidingWindows_new(SEXP sgenoSEXP, SEXP dSEXP, SEXP vSEXP, SEXP popSEXP, SEXP popUniqueSEXP, SEXP admixedSEXP, SEXP window_sizeSEXP, SEXP mapSEXP, SEXP with_mapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type sgeno(sgenoSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const StringVector& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< const StringVector& >::type popUnique(popUniqueSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type admixed(admixedSEXP);
    Rcpp::traits::input_parameter< const int >::type window_size(window_sizeSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const int >::type with_map(with_mapSEXP);
    rcpp_result_gen = Rcpp::wrap(slidingWindows_new(sgeno, d, v, pop, popUnique, admixed, window_size, map, with_map));
    return rcpp_result_gen;
END_RCPP
}
// get_fitted_matrix
arma::mat get_fitted_matrix(arma::mat& Y, arma::mat& U);
RcppExport SEXP _pcadapt_get_fitted_matrix(SEXP YSEXP, SEXP USEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type U(USEXP);
    rcpp_result_gen = Rcpp::wrap(get_fitted_matrix(Y, U));
    return rcpp_result_gen;
END_RCPP
}
// clumping_cpp
LogicalVector clumping_cpp(NumericMatrix& G, const IntegerVector& ord, LogicalVector& remain, const NumericVector& sumX, const NumericVector& denoX, int size, double thr);
RcppExport SEXP _pcadapt_clumping_cpp(SEXP GSEXP, SEXP ordSEXP, SEXP remainSEXP, SEXP sumXSEXP, SEXP denoXSEXP, SEXP sizeSEXP, SEXP thrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type ord(ordSEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type remain(remainSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sumX(sumXSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type denoX(denoXSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type thr(thrSEXP);
    rcpp_result_gen = Rcpp::wrap(clumping_cpp(G, ord, remain, sumX, denoX, size, thr));
    return rcpp_result_gen;
END_RCPP
}
// pca_rotation
arma::mat pca_rotation(arma::mat& a, arma::mat& b);
RcppExport SEXP _pcadapt_pca_rotation(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(pca_rotation(a, b));
    return rcpp_result_gen;
END_RCPP
}
// get_size_cpp
NumericVector get_size_cpp(std::string filename);
RcppExport SEXP _pcadapt_get_size_cpp(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(get_size_cpp(filename));
    return rcpp_result_gen;
END_RCPP
}
// get_nb_ind
int get_nb_ind(const StringVector& pop, const CharacterVector& name);
RcppExport SEXP _pcadapt_get_nb_ind(SEXP popSEXP, SEXP nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const StringVector& >::type pop(popSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type name(nameSEXP);
    rcpp_result_gen = Rcpp::wrap(get_nb_ind(pop, name));
    return rcpp_result_gen;
END_RCPP
}
// cmpt_minor_af
NumericVector cmpt_minor_af(arma::mat& xmatrix, int ploidy);
RcppExport SEXP _pcadapt_cmpt_minor_af(SEXP xmatrixSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type xmatrix(xmatrixSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    rcpp_result_gen = Rcpp::wrap(cmpt_minor_af(xmatrix, ploidy));
    return rcpp_result_gen;
END_RCPP
}
// scale_geno
arma::mat scale_geno(arma::mat& xmatrix, int ploidy, arma::vec maf, arma::vec keep_or_not);
RcppExport SEXP _pcadapt_scale_geno(SEXP xmatrixSEXP, SEXP ploidySEXP, SEXP mafSEXP, SEXP keep_or_notSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type xmatrix(xmatrixSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type keep_or_not(keep_or_notSEXP);
    rcpp_result_gen = Rcpp::wrap(scale_geno(xmatrix, ploidy, maf, keep_or_not));
    return rcpp_result_gen;
END_RCPP
}
// cmpt_loadings
arma::mat cmpt_loadings(std::string filename, arma::mat& xmatrix, arma::mat& scores, int nIND, int nSNP, int K, int ploidy, double min_maf, arma::vec& sigma, int type);
RcppExport SEXP _pcadapt_cmpt_loadings(SEXP filenameSEXP, SEXP xmatrixSEXP, SEXP scoresSEXP, SEXP nINDSEXP, SEXP nSNPSEXP, SEXP KSEXP, SEXP ploidySEXP, SEXP min_mafSEXP, SEXP sigmaSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xmatrix(xmatrixSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type scores(scoresSEXP);
    Rcpp::traits::input_parameter< int >::type nIND(nINDSEXP);
    Rcpp::traits::input_parameter< int >::type nSNP(nSNPSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< double >::type min_maf(min_mafSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(cmpt_loadings(filename, xmatrix, scores, nIND, nSNP, K, ploidy, min_maf, sigma, type));
    return rcpp_result_gen;
END_RCPP
}
// lrfunc_cpp
Rcpp::List lrfunc_cpp(std::string filename, arma::mat& xmatrix, arma::mat& scores, int nIND, int nSNP, int K, int ploidy, double min_maf, arma::vec& sigma, int type);
RcppExport SEXP _pcadapt_lrfunc_cpp(SEXP filenameSEXP, SEXP xmatrixSEXP, SEXP scoresSEXP, SEXP nINDSEXP, SEXP nSNPSEXP, SEXP KSEXP, SEXP ploidySEXP, SEXP min_mafSEXP, SEXP sigmaSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xmatrix(xmatrixSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type scores(scoresSEXP);
    Rcpp::traits::input_parameter< int >::type nIND(nINDSEXP);
    Rcpp::traits::input_parameter< int >::type nSNP(nSNPSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< double >::type min_maf(min_mafSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(lrfunc_cpp(filename, xmatrix, scores, nIND, nSNP, K, ploidy, min_maf, sigma, type));
    return rcpp_result_gen;
END_RCPP
}
// sample_geno_file
NumericVector sample_geno_file(std::string input, std::string output, double ploidy, IntegerVector sample_size);
RcppExport SEXP _pcadapt_sample_geno_file(SEXP inputSEXP, SEXP outputSEXP, SEXP ploidySEXP, SEXP sample_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type input(inputSEXP);
    Rcpp::traits::input_parameter< std::string >::type output(outputSEXP);
    Rcpp::traits::input_parameter< double >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sample_size(sample_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_geno_file(input, output, ploidy, sample_size));
    return rcpp_result_gen;
END_RCPP
}
// vcf_convert
IntegerVector vcf_convert(CharacterMatrix string_geno, std::string output, CharacterVector allele_sep);
RcppExport SEXP _pcadapt_vcf_convert(SEXP string_genoSEXP, SEXP outputSEXP, SEXP allele_sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterMatrix >::type string_geno(string_genoSEXP);
    Rcpp::traits::input_parameter< std::string >::type output(outputSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type allele_sep(allele_sepSEXP);
    rcpp_result_gen = Rcpp::wrap(vcf_convert(string_geno, output, allele_sep));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pcadapt_af", (DL_FUNC) &_pcadapt_af, 3},
    {"_pcadapt_bedXPtr", (DL_FUNC) &_pcadapt_bedXPtr, 3},
    {"_pcadapt_bedadaptXPtr", (DL_FUNC) &_pcadapt_bedadaptXPtr, 3},
    {"_pcadapt_cmpt_af", (DL_FUNC) &_pcadapt_cmpt_af, 1},
    {"_pcadapt_prodMatVec", (DL_FUNC) &_pcadapt_prodMatVec, 5},
    {"_pcadapt_prodtMatVec", (DL_FUNC) &_pcadapt_prodtMatVec, 5},
    {"_pcadapt_linReg", (DL_FUNC) &_pcadapt_linReg, 5},
    {"_pcadapt_get_sumX", (DL_FUNC) &_pcadapt_get_sumX, 3},
    {"_pcadapt_get_denoX", (DL_FUNC) &_pcadapt_get_denoX, 4},
    {"_pcadapt_clumping", (DL_FUNC) &_pcadapt_clumping, 9},
    {"_pcadapt_cmpt_cov_cpp", (DL_FUNC) &_pcadapt_cmpt_cov_cpp, 6},
    {"_pcadapt_cart2bary_cpp", (DL_FUNC) &_pcadapt_cart2bary_cpp, 2},
    {"_pcadapt_impute_geno", (DL_FUNC) &_pcadapt_impute_geno, 1},
    {"_pcadapt_impute_geno_pop", (DL_FUNC) &_pcadapt_impute_geno_pop, 3},
    {"_pcadapt_get_window", (DL_FUNC) &_pcadapt_get_window, 3},
    {"_pcadapt_cmpt_global_pca", (DL_FUNC) &_pcadapt_cmpt_global_pca, 3},
    {"_pcadapt_cmpt_local_pca", (DL_FUNC) &_pcadapt_cmpt_local_pca, 5},
    {"_pcadapt_updt_local_scores", (DL_FUNC) &_pcadapt_updt_local_scores, 8},
    {"_pcadapt_multLinReg", (DL_FUNC) &_pcadapt_multLinReg, 6},
    {"_pcadapt_cmpt_af_matrix", (DL_FUNC) &_pcadapt_cmpt_af_matrix, 1},
    {"_pcadapt_prodGx_matrix", (DL_FUNC) &_pcadapt_prodGx_matrix, 5},
    {"_pcadapt_prodtGx_matrix", (DL_FUNC) &_pcadapt_prodtGx_matrix, 5},
    {"_pcadapt_nb_nona", (DL_FUNC) &_pcadapt_nb_nona, 3},
    {"_pcadapt_covRob_cpp", (DL_FUNC) &_pcadapt_covRob_cpp, 1},
    {"_pcadapt_print_convert", (DL_FUNC) &_pcadapt_print_convert, 5},
    {"_pcadapt_ped2pcadapt", (DL_FUNC) &_pcadapt_ped2pcadapt, 2},
    {"_pcadapt_lfmm2pcadapt", (DL_FUNC) &_pcadapt_lfmm2pcadapt, 2},
    {"_pcadapt_sample_geno_matrix", (DL_FUNC) &_pcadapt_sample_geno_matrix, 3},
    {"_pcadapt_pMatVec4", (DL_FUNC) &_pcadapt_pMatVec4, 5},
    {"_pcadapt_cpMatVec4", (DL_FUNC) &_pcadapt_cpMatVec4, 5},
    {"_pcadapt_get_pop_size", (DL_FUNC) &_pcadapt_get_pop_size, 2},
    {"_pcadapt_slidingWindows_fast", (DL_FUNC) &_pcadapt_slidingWindows_fast, 9},
    {"_pcadapt_slidingWindows_new", (DL_FUNC) &_pcadapt_slidingWindows_new, 9},
    {"_pcadapt_get_fitted_matrix", (DL_FUNC) &_pcadapt_get_fitted_matrix, 2},
    {"_pcadapt_clumping_cpp", (DL_FUNC) &_pcadapt_clumping_cpp, 7},
    {"_pcadapt_pca_rotation", (DL_FUNC) &_pcadapt_pca_rotation, 2},
    {"_pcadapt_get_size_cpp", (DL_FUNC) &_pcadapt_get_size_cpp, 1},
    {"_pcadapt_get_nb_ind", (DL_FUNC) &_pcadapt_get_nb_ind, 2},
    {"_pcadapt_cmpt_minor_af", (DL_FUNC) &_pcadapt_cmpt_minor_af, 2},
    {"_pcadapt_scale_geno", (DL_FUNC) &_pcadapt_scale_geno, 4},
    {"_pcadapt_cmpt_loadings", (DL_FUNC) &_pcadapt_cmpt_loadings, 10},
    {"_pcadapt_lrfunc_cpp", (DL_FUNC) &_pcadapt_lrfunc_cpp, 10},
    {"_pcadapt_sample_geno_file", (DL_FUNC) &_pcadapt_sample_geno_file, 4},
    {"_pcadapt_vcf_convert", (DL_FUNC) &_pcadapt_vcf_convert, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_pcadapt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
