#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>

#define PRECISION_FLOAT 16
#define NA 9
#define MAXFILE 64

int initializeVariables__f(double **U, double **Sigma, double **V, double **SNPSd, double **Cov, double **miss, double **mAF, int K, int nSNP, int nIND);

void writeMatrix__f(double *U, double *Sigma, double *V, int nSNP, int K, int nIND, int nF, char *OutputFileName);

