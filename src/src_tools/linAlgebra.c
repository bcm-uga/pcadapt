/*
   PCAdapt linAlgebra.c
   Copyright (C) 2014 Nicolas Duforet Frebourg

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include "../src_Lapack/lpk.h"
#include "matrix.h"
#include "linAlgebra.h"

#define NA 9

void invMatrix(double *Matrix, int nrow){

	if (nrow == 1){ 
		Matrix[0] = 1.0/Matrix[0];
	} else {
        long int m= (long int) nrow, n= (long int) nrow;
        long int *pivot = calloc(((long int) nrow), sizeof(long int));
        long int lda= (long int) nrow;
        double *work = calloc((nrow)*(nrow), sizeof(double));
        long int info;
        long int length = nrow*nrow;

        dgetrf_(&m, &n, Matrix, &lda, pivot, &info);
        if (info > 0){
            Rprintf("WARNING singular matrix inversion...%d\n", info);
        }
        dgetri_(&m, Matrix, &lda, pivot, work, &length, &info);
        free(pivot);
        free(work);
	}
}


double SVD(double *G, double *Factors, double *Lambda, int K, int nSNP, int nIND){

    double *U = calloc(nSNP*nSNP, sizeof(double));
    double *Vt = calloc(nIND*nIND, sizeof(double));
    double *A = calloc(nIND*nSNP, sizeof(double));
    double *FLambda = calloc(nIND*nSNP, sizeof(double));
    long int i, j, m = nSNP, n = nIND;
    long int min_mn = (m < n) ? m : n;
    long int lda= m, ldu=m, ldv=n, lwork=m*n;
    long int *iwork = calloc(8*min_mn, sizeof(long int));
    double *singValues = calloc(min_mn, sizeof(double));
    double *work = calloc(lwork, sizeof(double));
    double sqerr = 0;
    double *err = calloc(n, sizeof(double));
    long int info;

    for (i=0; i<nSNP; i++){
        for (j=0; j<nIND; j++){
            A[j*nSNP + i] = G[i*nIND + j];
            }
    }
//    int res = dgesvd_("a", "a", (integer *) &m, (integer *) &n, (doublereal *) A, (integer *) &lda, (doublereal *) singValues, (doublereal *) U, (integer *) &ldu, (doublereal *) Vt, (integer *) &ldv, (doublereal *) work, (integer *) &lwork, (integer *) &info);
    dgesvd_("a", "a", (integer *) &m, (integer *) &n, (doublereal *) A, (integer *) &lda, (doublereal *) singValues, (doublereal *) U, (integer *) &ldu, (doublereal *) Vt, (integer *) &ldv, (doublereal *) work, (integer *) &lwork, (integer *) &info);
    for (i=0; i<nSNP; i++){
        for (j=0; j<K; j++){
            Factors[i*K + j] = U[j*ldu + i]*singValues[j];
        }
    }

    for (i=0; i<K; i++){
        for (j=0; j<nIND; j++){
            Lambda[i*nIND + j] = Vt[j*ldv + i];
        }
    }

    prodMatrix(Factors, Lambda, FLambda, nSNP, K, K, nIND);
    for (i=0; i<nSNP*nIND; i++){
        sqerr += (FLambda[i] - G[i])*(FLambda[i] - G[i])/(nSNP*nIND);
    }
    free(U);
    free(err);
    free(Vt);
    free(work);
    free(iwork);
    free(singValues);
    free(FLambda);
    free(A);
    return sqerr;
}

void tAA(double *A, double *tAA, int nrow, int ncol){

    /* matrix given A is column lead */
    int i;
    double *tA = malloc(sizeof(double)*nrow*ncol);
    for (i=0; i<nrow*ncol; i++){
        tA[i] = A[i];
    }
    long int nrtA = ncol, nctB = ncol, nrtB = nrow;
    long int ldc = nrtA, lda = nrtA, ldb = nrtB;
    double alpha = 1, beta = 0;
    tr(A, nrow, ncol);
//    int r = dgemm_("N", "N", &nrtA, &nctB, &nrtB, &alpha, tA, &lda, A, &ldb, &beta, tAA, &ldc);
    dgemm_("N", "N", &nrtA, &nctB, &nrtB, &alpha, tA, &lda, A, &ldb, &beta, tAA, &ldc);
    tr(A, ncol, nrow);
    free(tA);
}

/* Matrix must be symmetric */
void diagonalize(double *cov, int N, int K, double *val, double *vect){

    long int n = (long int)N;
    long int M = (long int)K;
    double abstol = 1e-10;
    long int *supp = (long int *) malloc (2 * N * sizeof(long int));
    long int lwork = 26*N;
    double *work = (double *)calloc(lwork, sizeof(double));
    long int liwork = 10*N;
    long int *iwork = (long int *)calloc(liwork, sizeof(double));
    long int info;
    double vl = 0.0, vu = 0.0;
    char jobz = 'V', range = 'I', uplo = 'U';
    long int il =  (long int)N-K+1;
    long int ul = (long int)N;
    double *valp = (double *) calloc(N, sizeof(double));
    double *vectp = (double *) calloc(N * N, sizeof(double));
    int i, k;
	double trCov = 0;

    for (i=0; i<N; i++){
        trCov += cov[i*N + i];
    }

    dsyevr_((char *) (&jobz), (char *) (&range) , (char *) (&uplo), (integer *) (&n), (doublereal *) cov, (integer *) (&n), (doublereal *) (&vl), (doublereal *) (&vu), (integer *) (&il) , (integer *) (&ul), (doublereal *) (&abstol), (integer *) (&M), (doublereal *) valp, (doublereal *) vectp, (integer *) (&n), (integer *)supp, (doublereal *)work, (integer *) (&lwork), (integer *)iwork, (integer *) (&liwork), (integer *) (&info));

    for (k = 0; k < K; k++){
        //Rprintf("Percentage of variance in PC%i: %g\n", k + 1, valp[K - (k + 1)]/trCov);
        val[k] = sqrt(valp[K-(k+1)]);
	}

    for (k = 0; k < K; k++){
        for (i = 0; i < N; i++){
            vect[i * K + k] = vectp[(K - (k + 1)) * N + i];
        }
    }
    free(valp);
    free(vectp);
    free(supp);
    free(work);
    free(iwork);
}

