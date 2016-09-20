/*
   FastPCAdapt FastPCAdapt.c
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
#include <unistd.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include "src_tools/Cov_line.h"
#include "src_tools/Data__f.h"
#include "src_tools/linAlgebra.h"
#include "src_tools/matrix.h"
#include "pcadapt.h"
#define MAXFILE 64

void wrapper_pcadapt(char **inputfilename,int *npc, double *minmaf, int *ploidy, char **outputfilename){
	int number_of_pc;
	char *name_of_the_input;
	double minor_allele_freq_threshold;
	int haploid_logical;
    char *name_of_the_output;
	number_of_pc = *npc;
	name_of_the_input = *inputfilename;
	minor_allele_freq_threshold = *minmaf;
	haploid_logical = *ploidy;
    name_of_the_output = *outputfilename;
	FastPCAdapt(name_of_the_input,number_of_pc,minor_allele_freq_threshold,haploid_logical,name_of_the_output);
}

void FastPCAdapt(char *filename, int K, double min_AF, int haploid, char *OUTFILE){

    FILE *GenoFile;
	char *OutputFileName = malloc(sizeof(char)*256);
	int nSNP, nIND, nF = K;
	int sc = 1, i, snp;
	double *U, *V, *Sigma, *SNPSd, *Cov, *mAF, *miss;
	Rprintf("Reading file %s...\n",filename);
	FILE *INPUTFile;
	int nbbn = 0;
	int nbsp = 0;
	int res;
	if((INPUTFile = fopen(filename, "r")) == NULL){
		Rprintf("Error, invalid input file.\n");
	}
	int currentchar;
	int prevchar;
	currentchar = fgetc(INPUTFile);
	while(currentchar != EOF){
		if (currentchar == 10){
			nbbn++;
	        if (prevchar != 32 && prevchar != '\t'){
	        	nbsp++;
	        }
	    }
	    if ((currentchar == 32 || prevchar == '\t') && (prevchar != 32 || prevchar != '\t')) nbsp++;
	    prevchar = currentchar;
	    currentchar = fgetc(INPUTFile);
	}
	fclose(INPUTFile);
	snp = (int) nbbn;
	res = nbsp/nbbn;
	nIND = res;
	nSNP = snp;
	Rprintf("Number of SNPs: %i\n",snp);
	Rprintf("Number of individuals: %i\n",nIND);

	/* Allocate memory */
    initializeVariables__f(&U, &Sigma, &V, &SNPSd, &Cov, &miss, &mAF, nF, nSNP, nIND);
    OutputFileName = OUTFILE;

	/* Compute nxn covariance matrix */
	double mean;
	int i1, na1, na_tot = 0, low_AF_tot = 0;
	int blocksize = 120;
	double *Geno = calloc(nIND*blocksize, sizeof(double));
	double *scratchCov = calloc(nIND*nIND, sizeof(double));
	blocksize = 120;
	if((GenoFile = fopen(filename, "r")) == NULL){
		Rprintf("Error, invalid input file.\n");
	}
	GenoFile = fopen(filename,"r");
	for (i1=0; i1<nSNP; i1 += blocksize){
		if (nSNP - i1 < blocksize) blocksize = nSNP - i1;
		na1 = get_row(Geno, GenoFile, nIND, &mean, SNPSd + i1, sc, blocksize, haploid, min_AF, &low_AF_tot);
		add_to_cov(Cov, scratchCov, nIND, Geno, blocksize);
		na_tot += na1;
	}
	// CHECK BLOCSIZE AND GENO
	Rprintf("Number of SNPs with minor allele frequency lower than %g ignored: %i\n", min_AF,low_AF_tot);
	if (na_tot) Rprintf("%i out of %i missing data ignored.\n", na_tot, nSNP*nIND);
    /* END of Cov_line*/

	/* Covariance matrix diagonalization */
	/* Compute singular values */
	/* Compute scores nxK matrix V */
	diagonalize(Cov, nIND, K, Sigma, V);
	/* END */

	/* Compute loadings Kxp matrix U */
	int i2, j, na2, low_AF_tot2 = 0;
	double var2, mean2;
	for (i=0; i<nIND; i++){
		for (j=0; j<K; j++){
			V[i*K + j] /= Sigma[j];
        }
	}
    if((GenoFile = fopen(filename, "r")) == NULL){
    	Rprintf("Error, invalid input file.\n");
    }
    //double *Geno = calloc(nIND, sizeof(double));
	for (i2=0; i2<nSNP; i2++){
		na2 = get_row(Geno, GenoFile, nIND, &mean2, &var2, sc, 1, haploid, min_AF, &low_AF_tot2);
		miss[i2] = na2;
		mAF[i2] = mean2;
		prodMatrix(Geno, V, (U + K*(i2)), 1, nIND, nIND, K);
	}
	fclose(GenoFile);
	int i3, j3;
	for (i3=0; i3<nIND; i3++){
		for (j3=0; j3<K; j3++){
			V[i3*K + j3] *= Sigma[j3];
        }
    }
	/*END*/

	/* Linear regression */
	int i4=0, na4=0, low_AF_tot4 = 0, i6;
	double *Z,*Ypred, *tV, *residuals;
	Z = malloc(sizeof(double)*nSNP*K);
	Ypred = malloc(sizeof(double)*nIND);
	tV = calloc(K*nIND, sizeof(double));
    for (i6 = 0;i6<(K*nIND);i6++){
    	tV[i6] = V[i6];
    }
    tr(tV,nIND,K);
	residuals = malloc(sizeof(double)*nSNP);
	double var4, mean4;
    if((GenoFile = fopen(filename, "r")) == NULL){
    	Rprintf("Error, invalid input file.\n");
    }
    double *Genotype = calloc(nIND, sizeof(double));
    int ax;
	for (i4=0; i4<nSNP; i4++){
		residuals[i4] = 0;
		na4 = get_row(Genotype, GenoFile, nIND, &mean4, &var4, 0, 1, haploid, min_AF, &low_AF_tot4);
		miss[i4] = na4;
		mAF[i4] = mean4;
		prodMatrix(Genotype, V, (Z + K*(i4)), 1, nIND, nIND, K);
		prodMatrix((Z+K*(i4)),tV,Ypred,1,K,K,nIND);
	    for (ax=0;ax<nIND;ax++){
	        residuals[i4] += (Genotype[ax]-Ypred[ax])*(Genotype[ax]-Ypred[ax]);
	    }
	    if (nIND-K-na4<=0){
	        residuals[i4] = 0.0;
	    } else {
	        residuals[i4] /= (double) nIND-K-na4;
	    }
	}
	fclose(GenoFile);
    
    
    /* Correcting for missing values */
    if((GenoFile = fopen(filename, "r")) == NULL){
        Rprintf("Error, invalid input file.\n");
    }
    int ii=0,jj=0,ll=0,lb=0;
    double *xx = calloc(nF, sizeof(double));
    float value;
    double *info_na = calloc(nIND, sizeof(double));
    
    for (ii=0;ii<nSNP;ii++){
        for (lb=0;lb<nF;lb++){
            xx[lb]=0;
        }
        for (ll=0;ll<nF;ll++){
            for (jj=0;jj<nIND;jj++){
                if (ll==0){
                    if(fscanf(GenoFile, "%g", &value) != EOF){
                        if(value != NA){
                            xx[ll] += (double) V[nF*jj+ll]*V[nF*jj+ll];
                            info_na[jj] = 1;
                        } else {
                            info_na[jj] = 0;
                        }
                    }
                } else {
                    if (info_na[jj] == 1){
                        xx[ll] += (double) V[nF*jj+ll]*V[nF*jj+ll];
                    }
                }
            }
            if (xx[ll] > 0){
                Z[ii*nF+ll] /= sqrt(xx[ll]);
            }
        }
        //Rprintf("aux1 = %f, aux2 = %f, aux3 = %f\n",xx[0],xx[1],xx[2]);
    }
    fclose(GenoFile);
    
    
	FILE *zFile;
    FILE *mafFile;
    char *fileMAF = malloc(sizeof(char)*256);
	char *fileZ = malloc(sizeof(char)*256);
    strcpy(fileZ, OutputFileName);
    strcat(fileZ, ".zscores");
    strcpy(fileMAF, OutputFileName);
    strcat(fileMAF, ".maf");
    if((zFile = fopen(fileZ, "w")) == NULL) Rprintf("Error, unable to open %s.\n", fileZ);
    if((mafFile = fopen(fileMAF, "w")) == NULL) Rprintf("Error, unable to open %s.\n", fileMAF);
    int i5,j5;
    for (i5=0; i5<nSNP; i5++){
        fprintf(mafFile,"%f",mAF[i5]);
    	for (j5=0; j5<nF; j5++){
            if (residuals[i5]==0.0){
                fprintf(zFile, "NA ");
            } else {
                fprintf(zFile, "%f ", Z[i5*nF + j5]/sqrt(residuals[i5]));
            }
        }
        fprintf(mafFile,"\n");
        fprintf(zFile, "\n");
    }
    fclose(zFile);
    fclose(mafFile);
	/* now we refer to V as tV, the K * nIND matrix */
	tr(V, nIND, K);

    writeMatrix__f(U, Sigma, V, nSNP, K, nIND, nF, OutputFileName);
    
    free(Geno);
	free(miss);
	free(mAF);
	free(Cov);
}

