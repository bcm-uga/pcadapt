#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

#define NA 9

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat fJ_cpp(int n){
  arma::mat zz(n, n);
  zz.ones();
  zz /= n;
  arma::mat H(n, n);
  H.eye();
  H -= zz;
  return(H);
}

// [[Rcpp::export]]
arma::mat fcnt_cpp(arma::mat &a){
  int nrow = a.n_rows;
  arma::mat aa = fJ_cpp(nrow);
  return(aa * a);
}

// [[Rcpp::export]]
arma::mat pca_rotation(arma::mat &a, arma::mat &b){
  arma::mat fcnt_a = fcnt_cpp(a);
  arma::mat fcnt_b = fcnt_cpp(b);
  arma::mat x = fcnt_a.t() * fcnt_b;
  arma::mat u;
  arma::vec s;
  arma::mat v;
  svd(u, s, v, x);
  arma::mat R;
  R = v * u.t();
  arma::cx_vec eigval_v;
  arma::cx_mat eigvec_v;
  arma::cx_vec eigval_u;
  arma::cx_mat eigvec_u;
  eig_gen(eigval_v, eigvec_v, v);
  eig_gen(eigval_u, eigvec_u, u);
  arma::cx_double tmp_v = prod(eigval_v);
  arma::cx_double tmp_u = prod(eigval_u);
  double chk1 = real(tmp_v);
  double chk2 = real(tmp_u);
  if ((chk1 < 0) && (chk2 > 0)) {
    for (int i = 0; i < v.n_rows; i++){
      v(i, v.n_cols - 1) *= (-1);
    }
    R = v * u.t();
  }
  if ((chk2 < 0) && (chk1 > 0)) {
    for (int i = 0; i < u.n_rows; i++){
      u(i, u.n_cols - 1) *= (-1);
    }
    R = v * u.t();
  }
  return(R);
}


//' Number of individuals in a specific population
//' 
//' \code{get_nb_ind} returns the number of individuals in a specific population.
//' 
//' @param lab a vector of integers.
//' @param anc an integer.
//' 
//' @return The returned value is an integer.
//' 
//' @export
//' 
// [[Rcpp::export]]
int get_nb_ind(const arma::vec &lab, const int anc){
  int n = lab.n_elem;
  int c = 0;
  for (int i = 0; i < n; i++){
    if (lab[i] == anc){
      c++;  
    }
  }
  return(c);
}


//' @export
//' 
// [[Rcpp::export]]
Rcpp::List cmpt_centroids(arma::mat u, arma::vec lab, int anc1, int anc2){
  int nIND = u.n_rows;
  int K = u.n_cols;
  arma::vec m1(K);
  m1.zeros();
  arma::vec m2(K);
  m2.zeros();
  int c1 = get_nb_ind(lab, anc1);
  int c2 = get_nb_ind(lab, anc2);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      if (lab[j] == anc1){
        m1[k] += u(j, k) / c1;
      } else if (lab[j] == anc2){
        m2[k] += u(j, k) / c2;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("m1") = m1,
                            Rcpp::Named("m2") = m2);
}

//' @export
//' 
// [[Rcpp::export]]
void cmpt_transformation(arma::mat &uloc, 
                         arma::mat &uglob, 
                         const arma::vec &lab, 
                         const int ancstrl1, 
                         const int ancstrl2, 
                         arma::vec &s, 
                         arma::vec &dloc, 
                         arma::vec &dglob, 
                         arma::mat &R){
  Rcpp::List mglob = cmpt_centroids(uglob, lab, ancstrl1, ancstrl2);
  Rcpp::List mloc = cmpt_centroids(uloc, lab, ancstrl1, ancstrl2);
  arma::vec mglob1 = mglob[0];
  arma::vec mglob2 = mglob[1];
  arma::vec mloc1 = mloc[0];
  arma::vec mloc2 = mloc[1];
  dglob = (mglob1 + mglob2);
  dglob /= 2.0;
  dloc =  (mloc1 + mloc2);
  dloc /= 2.0;
  int K = s.n_elem;
  for (int k = 0; k < K; k++){
    s[k] = fabs(mglob1[k] - mglob2[k]) / fabs(mloc1[k] - mloc2[k]);
  }
}

//' @export
//' 
// [[Rcpp::export]]
arma::mat rescale_local_pca(arma::mat &u, arma::vec &s, arma::vec &dep_loc, arma::vec &dep_glob, arma::mat &R){
  int nIND = u.n_rows;
  int K = u.n_cols;
  arma::mat usc(nIND, K);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      usc(j, k) = u(j, k) * s[k];
      usc(j, k) = usc(j, k) + dep_glob[k] - dep_loc[k];
    }
  }
  usc *= R;
  return(usc);
}


//' Global Principal Component Analysis
//' 
//' \code{cmpt_global_pca} computes the scores using all genetic markers.
//' 
//' @param geno a genotype matrix.
//' @param V a loading matrix.
//' @param sigma a vector of singular values.
//' 
//' @return The returned value is a matrix of scores.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::mat cmpt_global_pca(const arma::mat &geno, const arma::mat &V, const arma::vec &sigma){
  int nIND = geno.n_rows;
  int nSNP = V.n_rows;
  int K = V.n_cols;
  arma::mat u(nIND, K, arma::fill::zeros);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      for (int i = 0; i < nSNP; i++){
        u(j, k) += geno(j, i) * V(i, k) / sigma[k];
      }
    }
  }
  return(u);
}

//' Local Principal Component Analysis
//' 
//' \code{cmpt_local_pca} computes the scores using a subset of genetic markers.
//' 
//' @param geno a genotype matrix.
//' @param V a loading matrix.
//' @param sigma a vector of singular values.
//' @param beg an integer specifying the first marker to be included.
//' @param end an integer specifying the first marker to be excluded.
//' 
//' @return The returned value is a matrix of scores.
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::mat cmpt_local_pca(const arma::mat &geno, const arma::mat &V, const arma::vec &sigma, const int beg, const int end){
  // [beg, end) 
  int nIND = geno.n_rows;
  int nSNP = V.n_rows;
  int K = V.n_cols;
  arma::mat u(nIND, K, arma::fill::zeros);
  double cst = (double) nSNP / (end - beg);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      for (int i = beg; i < end; i++){
        u(j, k) += geno(j, i) * V(i, k) * cst / sigma[k];
      }
    }
  }
  return(u);
}

//' Update local Principal Component Analysis
//' 
//' \code{updt_local_scores} computes the scores using a subset of genetic markers.
//' 
//' @param u a score matrix.
//' @param geno a genotype matrix.
//' @param V a loading matrix.
//' @param sigma a vector of singular values.
//' @param beg an integer specifying the first marker to be included.
//' @param end an integer specifying the first marker to be excluded.
//' 
//' @return The returned value is a score matrix.
//' 
//' @export
//' 
// [[Rcpp::export]]
void updt_local_scores(arma::mat &u, const arma::mat &geno, const arma::mat &V, const arma::vec &sigma, const int beg, const int end){
  int nIND = geno.n_rows; 
  int nSNP = geno.n_cols;
  int K = u.n_cols;
  double cst = (double) nSNP / (end - beg);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      u(j, k) -= (geno.at(j, beg - 1) * V(beg - 1, k)) * cst / sigma[k];
      u(j, k) += (geno.at(j, end) * V(end, k)) * cst / sigma[k];
    }
  }
}


//' @export
//' 
// [[Rcpp::export]]
double cmpt_window_stat(arma::mat &uloc,
                        arma::mat &uglob, 
                        const int direction, 
                        const arma::vec &lab, 
                        const int adm, 
                        const int axis){
  int nIND = uglob.n_rows; 
  double stat = 0;
  if (direction == 1){
    for (int j = 0; j < nIND; j++){
      if ((lab[j] == adm) && (uloc(j, axis) - uglob(j, axis)) > 0){
        stat += (uloc(j, axis) - uglob(j, axis)) * (uloc(j, axis) - uglob(j, axis));
      }
    }
  } else if (direction == (-1)){
    for (int j = 0; j < nIND; j++){
      if ((lab[j] == adm) && (uloc(j, axis) - uglob(j, axis)) < 0){
        stat += (uloc(j, axis) - uglob(j, axis)) * (uloc(j, axis) - uglob(j, axis));
      }
    }
  } else if (direction == 0){
    for (int j = 0; j < nIND; j++){
      if (lab[j] == adm){
        stat += (uloc(j, axis) - uglob(j, axis));
      }
    }
  }
  return(stat);
}

//' @export
//' 
// [[Rcpp::export]]
arma::vec get_rank(const arma::vec &v_temp){
  int n = v_temp.n_elem;
  arma::vec v_sort(n);
  for (int i = 0; i < n; i++){
    v_sort[i] = v_temp[i];
  }
  arma::uvec idx = sort_index(v_sort);
  arma::vec rank(n);
  rank.zeros();
  for (int i = 0; i < n; i++){
    int tmp = idx[i];
    rank[tmp] = i + 1;
  }
  return(rank);
}

//' @export
//' 
// [[Rcpp::export]]
arma::vec get_axis(arma::mat &uglob, const arma::vec &lab, const int anc1, const int anc2){
  Rcpp::List res = cmpt_centroids(uglob, lab, anc1, anc2);
  arma::vec m1 = res[0];
  arma::vec m2 = res[1];
  return(m2 - m1);
}

//' @export
//' 
// [[Rcpp::export]]
double cmpt_directional_stat(arma::mat &usc,
                             arma::mat &uglob, 
                             const arma::vec &lab, 
                             const int adm, 
                             arma::vec &ax){
  
  int nIND = uglob.n_rows; 
  double stat = 0;
  for (int j = 0; j < nIND; j++){
    if (lab[j] == adm){
      stat += arma::dot(usc.row(j) - uglob.row(j), ax);
    }
  }
  return(stat);
}


//' @export
//' 
// [[Rcpp::export]]
double cmpt_window_wilcoxon(arma::mat &uloc,
                            arma::mat &uglob, 
                            int direction, 
                            arma::vec &lab, 
                            int adm, 
                            int axis){
  int nIND = uglob.n_rows;
  int nAND = get_nb_ind(lab, adm);
  arma::vec tmp(nIND);
  arma::vec Z(nIND);
  Z.zeros();
  double diff;
  double m = 0;
  for (int j = 0; j < nIND; j++){
    if (lab[j] == adm){
      m += uglob(j, axis) / nAND;
    }
  }
  
  for (int j = 0; j < nIND; j++){
    diff = uloc(j, axis) - m;
    tmp[j] = fabs(diff);
    if ((direction == 1) && (diff > 0)){
      Z[j] = 1;
    } else if ((direction == (-1)) && (diff < 0)){
      Z[j] = 1;
    }
  }
  arma::vec tmp_sort = get_rank(tmp);
  double W = 0;
  for (int j = 0; j < nIND; j++){
    if (lab[j] == adm){
      W += (double) Z[j] * tmp_sort[j];
    }
  }
  return(W);
}

//' @export
//' 
// [[Rcpp::export]]
arma::vec cmpt_all_stat(const arma::mat &geno, 
                        const arma::mat &V, 
                        const arma::vec &sigma, 
                        const int window_size,  
                        const int direction, 
                        const arma::vec lab, 
                        const int ancstrl1,
                        const int ancstrl2,
                        const int adm, 
                        const arma::vec axis){
  int nSNP = geno.n_cols;
  int nIND = geno.n_rows;
  int K = V.n_cols;
  arma::vec s(K, arma::fill::zeros);
  arma::vec dglob(K, arma::fill::zeros);
  arma::vec dloc(K, arma::fill::zeros);
  arma::mat usc(nIND, K, arma::fill::zeros);
  arma::vec stat(nSNP, arma::fill::zeros);
  arma::vec ax(K, arma::fill::zeros);
  
  arma::mat R(K, K); // ROTATION CORRECTION
  R.eye();
  
  arma::mat uglob = cmpt_global_pca(geno, V, sigma);
  arma::mat uloc = cmpt_local_pca(geno, V, sigma, 0, window_size);
  ax = get_axis(uglob, lab, ancstrl1, ancstrl2);
  
  for (int k = 0; k < K; k++){
    ax[k] *= axis[k];   
  }
  
  cmpt_transformation(uloc, uglob, lab, ancstrl1, ancstrl2, s, dloc, dglob, R);
  usc = rescale_local_pca(uloc, s, dglob, dloc, R);
  for (int i = 1; i < (nSNP - window_size); i++){
    updt_local_scores(uloc, geno, V, sigma, i, i + window_size);
    cmpt_transformation(uloc, uglob, lab, ancstrl1, ancstrl2, s, dloc, dglob, R);
    usc = rescale_local_pca(uloc, s, dloc, dglob, R);
    stat[i] = cmpt_directional_stat(usc, uglob, lab, adm, ax);
    //stat[i] = cmpt_window_stat(usc, uglob, direction, lab, adm, axis);
    //stat[i] = cmpt_window_wilcoxon(usc, uglob, direction, lab, adm, axis);
  }
  stat[0] = stat[1];
  for (int i = (nSNP - window_size); i < nSNP; i++){
    stat[i] = stat[nSNP - window_size - 1];
  }
  return(stat);
}
