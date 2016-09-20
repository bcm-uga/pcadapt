/*
   FastPCAdapt Data__f.c
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <R.h>
#include "Data__f.h"
#include "random.h"

#define PRECISION_FLOAT 16
#define NA 9

int initializeVariables__f(double **U, double **Sigma, double **V, double **SNPSd, double **Cov, double **miss, double **mAF, int K, int nSNP, int nIND){

	int i, j;
    init_random();

	*Sigma = calloc(K, sizeof(double));
    *U = malloc(sizeof(double)*nSNP*K);
    for (i=0; i<K; i++){
    	for (j=0; j<nSNP; j++){
    		*(*(U) + i*nSNP + j) = rand_normal(0, 1);
        }
    }
	*SNPSd = malloc(sizeof(double)*nSNP);
    *V = calloc(K*nIND, sizeof(double));
	*Cov = calloc(nIND*nIND, sizeof(double));
    for (i=0; i<K*nIND; i++){
        *(*(V) + i) = rand_normal_r();
    }
	*miss = calloc(nSNP, sizeof(double));
	*mAF = calloc(nSNP, sizeof(double));
	return 0;
}

void writeMatrix__f(double *U, double *Sigma, double *V, int nSNP, int K, int nIND, int nF, char *OutputFileName){

    FILE *resultFile;
	char *fileU = malloc(sizeof(char)*256);
    char *fileS = malloc(sizeof(char)*256);
    char *fileV = malloc(sizeof(char)*256);
    int i, j;
    strcpy(fileU, OutputFileName);
    strcat(fileU, ".loadings");

    if((resultFile = fopen(fileU, "w")) == NULL){
        Rprintf("Error, unable to open %s.\n", fileU);
    }
    for (i=0; i<nSNP; i++){
    	for (j=0; j<nF; j++){
    		fprintf(resultFile, "%f ", U[i*nF + j]);
        }
        fprintf(resultFile, "\n");
    }
	fclose(resultFile);

    strcpy(fileS, OutputFileName);
    strcat(fileS, ".sigma");

    if((resultFile = fopen(fileS, "w")) == NULL){
        Rprintf("Error, unable to open %s.\n", fileS);
    }
    for (i=0; i<nF; i++){
    	fprintf(resultFile, "%f ", Sigma[i]/sqrt(nIND - 1));
    }
	fprintf(resultFile, "\n");
    fclose(resultFile);

    strcpy(fileV, OutputFileName);
    strcat(fileV, ".scores");

    if((resultFile = fopen(fileV, "w")) == NULL){
        Rprintf("Error, unable to open %s.\n", fileV);
    }
    for (i=0; i<nF; i++){
    	for (j=0; j<nIND; j++){
    		fprintf(resultFile, "%f ", V[i*nIND + j]);
        }
        fprintf(resultFile, "\n");
    }
    fclose(resultFile);
	free(fileU);
    free(fileS);
    free(fileV);

}

