#ifndef RAND_H
#define RAND_H

//#define PI (3.141592653589793)

/*
 * Library of function to generate random variables or vectors, from different distributions.
 */

void init_random();

double drand();

int rand_int(int size);

double rand_double(double min, double max);

void rand_matrix_double(double* A, int M, int N);

double rand_normal( double mean, double var);

double rand_normal_r();

void mvn_rand(double* mu, double* L, int D, double* y);

int rand_vector(double *Pi, int size);

double rand_gamma(int alpha, double beta);

double rand_invgamma(int alpha, double beta);

double rand_chisq(int n);

double rand_beta(double alpha, double beta);

void rand_dirichlet(double *x, double *alpha, int K);

#endif // RAND_H 
