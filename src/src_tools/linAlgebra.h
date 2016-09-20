#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include "../src_Lapack/lpk.h"

/* invMatrix
 * Calculate the inverse of Matrix, in Matrix.
 */
void invMatrix(double *Matrix, int nrow);

/* SVD
 * Calculate Singular value decomposition of G, and writes U in Lambda, and Sigma * V in Factors.
 */
double SVD(double *G, double *Factors, double *Lambda, int K, int nSNP, int nIND);

/* tAA
 * Optimized way to calculate A tranpose * A in tAA.
 */
void tAA(double *A, double *tAA, int nrow, int ncol);

/* diagonalize
 * Diagoanlize matrix cov and write squarred root eigenvalues in val, and eigenvectors in vect.
 */
void diagonalize(double *cov, int N, int K, double *val, double *vect);

