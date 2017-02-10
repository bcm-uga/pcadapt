#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

#define NA 9

using namespace Rcpp;

//' @export
//' 
// [[Rcpp::export]]
arma::mat cmpt_global_pca(arma::mat &geno, arma::mat &V, arma::vec &sigma){
  int nIND = geno.n_rows;
  int nSNP = V.n_rows;
  int K = V.n_cols;
  arma::mat u(nIND, K);
  u.zeros();
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      for (int i = 0; i < nSNP; i++){
        u(j, k) += geno(j, i) * V(i, k) / sigma[k];
      }
    }
  }
  return(u);
}

//' @export
//' 
// [[Rcpp::export]]
arma::mat cmpt_local_pca(arma::mat &geno, arma::mat &V, arma::vec &sigma, int beg, int end){
  int nIND = geno.n_rows;
  int nSNP = V.n_rows;
  int K = V.n_cols;
  arma::mat uloc(nIND, K);
  uloc.zeros();
  double cst = (double) nSNP;
  cst /=  (double) end - beg;
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      for (int i = beg; i < end; i++){
        uloc(j, k) += geno(j, i) * V(i, k) * cst / sigma[k];
      }
    }
  }
  return(uloc);
}

//' @export
//' 
// [[Rcpp::export]]
void updt_local_scores(arma::mat &u, arma::mat &geno, arma::mat &V, arma::vec &sigma, int window_size, int i){
  int nIND = geno.n_rows; 
  int nSNP = geno.n_cols;
  int K = u.n_cols;
  int beg = i - 1;
  int end = i + window_size; 
  double cst = (double) nSNP;
  cst /=  (double) end - beg;
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      u(j, k) -= (geno.at(j, beg) * V(beg, k)) * cst / sigma[k];
      u(j, k) += (geno.at(j, end) * V(end, k)) * cst / sigma[k];
    }
  }
}

//' @export
//' 
// [[Rcpp::export]]
int get_nb_ind(arma::vec lab, int anc){
  int n = lab.n_elem;
  int c = 0;
  for (int i = 0; i < n; i++){
    if (lab[i] == anc){
      c += 1;  
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
void cmpt_transformation(arma::mat &uloc, arma::mat &uglob, arma::vec &lab, int ancstrl1, int ancstrl2, arma::vec &s, arma::vec &dloc, arma::vec &dglob){
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
arma::mat rescale_local_pca(arma::mat &u, arma::vec &s, arma::vec &dep_loc, arma::vec &dep_glob){
  int nIND = u.n_rows;
  int K = u.n_cols;
  arma::mat usc(nIND, K);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      usc(j, k) = u(j, k) * s[k];
      usc(j, k) = usc(j, k) + dep_glob[k] - dep_loc[k];
    }
  }
  return(usc);
}

//' @export
//' 
// [[Rcpp::export]]
double cmpt_window_stat(arma::mat &uloc,
                        arma::mat &uglob, 
                        int direction, 
                        arma::vec &lab, 
                        int adm, 
                        int axis){
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
arma::vec cmpt_all_stat(arma::mat &geno, 
                        arma::mat &V, 
                        arma::vec &sigma, 
                        int window_size,  
                        int direction, 
                        arma::vec lab, 
                        int ancstrl1,
                        int ancstrl2,
                        int adm, 
                        int axis){
  int nSNP = geno.n_cols;
  int nIND = geno.n_rows;
  int K = V.n_cols;
  arma::vec s(K);
  s.zeros();
  arma::vec dglob(K);
  dglob.zeros();
  arma::vec dloc(K);
  dloc.zeros();
  
  arma::vec stat(nSNP);
  
  arma::mat uglob = cmpt_global_pca(geno, V, sigma);
  arma::mat uloc = cmpt_local_pca(geno, V, sigma, 0, window_size);
  cmpt_transformation(uloc, uglob, lab, ancstrl1, ancstrl2, s, dloc, dglob);
  
  arma::mat usc(nIND, K);
  usc = rescale_local_pca(uloc, s, dglob, dloc);
  stat[0] = cmpt_window_stat(usc, uglob, direction, lab, adm, axis);
  //stat[0] = cmpt_window_wilcoxon(usc, uglob, direction, lab, adm, axis);
  
  for (int i = 1; i < (nSNP - window_size); i++){
    updt_local_scores(uloc, geno, V, sigma, window_size, i);
    cmpt_transformation(uloc, uglob, lab, ancstrl1, ancstrl2, s, dloc, dglob);
    usc = rescale_local_pca(uloc, s, dloc, dglob);
    stat[i] = cmpt_window_stat(usc, uglob, direction, lab, adm, axis);
    //stat[i] = cmpt_window_wilcoxon(usc, uglob, direction, lab, adm, axis);
  }
  for (int i = (nSNP - window_size); i < nSNP; i++){
    stat[i] = stat[nSNP - window_size - 1];
  }
  return(stat);
}
