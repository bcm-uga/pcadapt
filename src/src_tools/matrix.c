/*
   PCAdapt matrix.c
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

#define NA 9

void displayMatrix(double *Matrix, int nrow, int ncol){

    int i=0, j=0;
    for(i=0; i<nrow; i++){
        for(j=0; j<ncol; j++){
            Rprintf("%f ", Matrix[i*ncol + j]);
            }
        Rprintf("\n");
    }
}

void rowMeans(double *Matrix, int nrow, int ncol, double *Means){

	int row, col, na;
	for (row=0; row<nrow; row++){
		Means[row] = 0;
		na = 0;
		for (col=0; col<ncol; col++){
			if(Matrix[row*ncol+ col] != NA){
				Means[row]= Means[row] + Matrix[row*ncol+ col];
			} else {
				na++;
			}
		}
		Means[row] = Means[row]/(ncol - na);
	}
}

void colMeans(double *Matrix, int nrow, int ncol, double *Means){

    int row, col, na;
    for (col=0; col<ncol; col++){
        Means[col] = 0;
		na = 0;
        for (row=0; row<nrow; row++){
			if (Matrix[row*ncol + col] != NA){
				Means[col]+=Matrix[row*ncol + col];
			} else {
				na++;
			}
		}
		Means[col] = Means[col]/(nrow - na);
    }
}

void sumMatrix(double *A, double *B, double *res, int nrow, int ncol){

    int row, col;
    for (row=0; row<nrow; row++){
        for (col=0; col<ncol; col++){
			res[row*ncol + col] = A[row*ncol + col] + B[row*ncol + col];
            if (A[row*ncol + col] == NA || B[row*ncol + col] == NA){
                res[row*ncol + col] = NA;
            }
        }
    }
}

void difMatrix(double *A, double *B, double *res, int nrow, int ncol){
        
    int row, col;
    for (row=0; row<nrow; row++){
        for (col=0; col<ncol; col++){
			res[row*ncol + col] = A[row*ncol + col] - B[row*ncol + col];
            if ((A[row*ncol + col] == NA) || (B[row*ncol + col] == NA)) {
                res[row*ncol + col] = NA;
            }
		}
    }
}

void prodMatrix(double *A, double *B, double *res, int nrowA, int ncolA, int nrowB, int ncolB){

    long int nrtA = nrowA, nctB = ncolB, nctA = ncolA;
    long int ldc = nctB, lda = nctA, ldb = nctB;
    double alpha = 1, beta = 0;
    dgemmlig_("T", "T", &nrtA, &nctB, &nctA, &alpha, A, &lda, B, &ldb, &beta, res, &ldc);
}

void tr(double *Matrix, int nrow, int ncol){

	int i, j;
	double *tmp = calloc(nrow*ncol, sizeof(double));
	for (i=0; i<ncol*nrow; i++) tmp[i] = Matrix[i];
	for (i=0; i<nrow; i++){
		for (j=0; j<ncol; j++){
			Matrix[j*nrow+ i] = tmp[i*ncol + j];
		}
	}
	free(tmp);
}

int findMax(double *M, int nrow, int ncol){

	int i, imax = 0;
    for (i=0; i<ncol*nrow; i++){
		imax = ((fabs(M[i]) > fabs(M[imax])) ? i : imax);
    }
	return imax;
}

int findMin(double *M, int nrow, int ncol){

    int i, imin = 0;
    for (i=0; i<ncol*nrow; i++){
        if (fabs(M[i]) < fabs(M[imin])){
            imin = i;
        }
    }
    return imin;
}

void scale(double *Matrix, double *rowSd, int nrow, int ncol, int center, int std){

    double *rMeans = calloc(nrow, sizeof(double));
    double *rSdev = calloc(nrow, sizeof(double));
    int i, j;
    if (center){
        rowMeans(Matrix, nrow, ncol, rMeans);
        for (i=0; i<nrow; i++){
            for(j=0; j<ncol; j++){
                Matrix[i*ncol + j] -= rMeans[i];
            }
        }
        for (i=0; i<nrow; i++){
            for (j=0; j<ncol; j++){
                rSdev[i] += Matrix[i*ncol + j]*Matrix[i*ncol + j]/(ncol - 1);
            }
            rowSd[i] = sqrt(rSdev[i]);
            if (std){
                for (j=0; j<ncol; j++){
                    if(rowSd[i] != 0){
                        Matrix[i*ncol + j] /= rowSd[i];
                    }
                }
                rowSd[i] = 1;
            }
        }
    }
    free(rMeans);
    free(rSdev);
}

