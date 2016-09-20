#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include "../src_Lapack/lpk.h"

/* displayMatrix
 * display Matrix Matrix.
 */
void displayMatrix(double *Matrix, int nrow, int ncol);

/* rowMeans
 * Calculate the vector of row means of matrix in Means.
 */
void rowMeans(double *Matrix, int nrow, int ncol, double *Means);

/* colMeans
 * Calculate the vector of column means of matrix in Means.
 */
void colMeans(double *Matrix, int nrow, int ncol, double *Means);

/* sumMatrix
 * Sum of matrices A and B.
 */
void sumMatrix(double *A, double *B, double *res, int nrow, int ncol);

/* difMatrix
 * Difference between matrices A and B.
 */
void difMatrix(double *A, double *B, double *res, int nrow, int ncol);

/* prodMatrix
 * Calculate A * B in res.
 */
void prodMatrix(double *A, double *B, double *res, int nrowA, int ncolA, int nrowB, int ncolB);

/* tr
 * Calculate the transpose of Matrix, in Matrix.
 */
void tr(double *Matrix, int nrow, int ncol);

/* findMax
 * find Maximum of matrice M and return its indice.
 */
int findMax(double *M, int nrow, int ncol);

/* findMin
 * find Minimum of matrice M and return its indice.
 */
int findMin(double *M, int nrow, int ncol);

/* scale
 * Centers and scales (if std) a matrix by lines 
 */
void scale(double *Matrix, double *rowSd, int nrow, int ncol, int center, int std);

