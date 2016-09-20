/*
 *    FastPCAdapt Cov_line.c
 *    Copyright (C) 2014 Nicolas Duforet Frebourg
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    this program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <R.h>
#include "math.h"
#include "Cov_line.h"
#include "Data__f.h"
#include "linAlgebra.h"
#include "matrix.h"

#define PRECISION_FLOAT 16
#define NA 9
#define TOL .0001

//get_row
//Read a block of lines (snp), scale it and store it.

int get_row(double *Geno, FILE *GenoFile, int nIND, double *mean, double *SNPSd, int sc, int blocksize, int haploid, double min_AF, int *low_AF_tot){

	float value;
	int snp, ind, na_tot = 0, na, low_AF = 0;
	double var, maf;
	for(snp=0; snp<blocksize; snp++){
		ind = 0;
		var = 0;
		na = 0;
		*mean = 0;
		maf = 0;
		while(ind < nIND){
            if(fscanf(GenoFile, "%g", &value) != EOF){
                Geno[snp*nIND + ind] = (double) value;
                if(value != NA){
                    *mean += (double) value;
                } else {
                    na++;
                }
            }
            ind++;
		}
		if (nIND <= na)	{
			*mean = NA;
			maf = NA;
		} else {
			*mean /= (nIND - na);
			maf = *mean;
			if (*mean > 1) maf = 2 - *mean;
			maf = (double) maf/2.0;
			if (haploid){
				maf = *mean;
				if(*mean > .5) maf = 1.0 - *mean;
			}
		}
		if (maf >= min_AF){
			for (ind=0; ind<nIND; ind++) {
                if(Geno[snp*nIND + ind] != NA){
                    Geno[snp*nIND + ind] -= *mean;
                    var += Geno[snp*nIND + ind]*Geno[snp*nIND + ind];
                } else {
                    Geno[snp*nIND + ind] = 0;
                }
            }
			if (sc){
				if (var > TOL){
                    //Use parametric variance
					var = 2*maf*(1.0 - maf);
                    if (haploid){
                        var = maf*(1 - maf);
                    }
					for (ind=0; ind<nIND; ind++){
                        Geno[snp*nIND + ind] /= sqrt(var);
                    }
					SNPSd[snp] = sqrt(var);
				} else {
					SNPSd[snp] = 1;
				}
			} else {
				SNPSd[snp] = sqrt(var/(nIND - na));
			}
			na_tot += na;
		} else {
			SNPSd[snp] = 1;
			for (ind=0; ind<nIND; ind++) Geno[snp*nIND + ind] = 0;
			low_AF++;
		}
	}
	*mean = maf;
	*low_AF_tot += low_AF;
	return na_tot;
}

//add_to_cov
//Given a block of snps, computes the correlation matrix of the block and add it to the existing matrix.

void add_to_cov(double *Cov, double *scratchCov, int nIND, double *Geno, int blocksize){

	tAA(Geno, scratchCov, blocksize, nIND);
	int i;
    for (i=0; i<nIND*nIND; i++){
        Cov[i] += scratchCov[i];
    }
}

//Cov_line
//Learn covariance (or cor) matrix, while reading a data file, by storing blocks of lines.

int Cov_line(double *Cov, double *SNPSd, int nSNP, int *nSNP_file, int nIND, int sc, char **GenoFileName, int nfile, int haploid, double min_AF){
	FILE *GenoFile;
	double mean;
	int i, na, na_tot = 0, file, low_AF_tot = 0;
	int blocksize = 120, snp_count = 0;
	double *Geno = calloc(nIND*blocksize, sizeof(double));
	double *scratchCov = calloc(nIND*nIND, sizeof(double));

	for (file=0; file<nfile; file++){
		blocksize = 120;
        if((GenoFile = fopen(GenoFileName[file], "r")) == NULL){
            Rprintf("Error, invalid input file\n");
            return 1;
        }
		for (i=0; i<nSNP_file[file] ; i += blocksize){
			if (nSNP_file[file] - i < blocksize) blocksize = nSNP_file[file] - i;
			na = get_row(Geno, GenoFile, nIND, &mean, SNPSd + i + snp_count, sc, blocksize, haploid, min_AF, &low_AF_tot);
			add_to_cov(Cov, scratchCov, nIND, Geno, blocksize);
			na_tot += na;
		}
		snp_count += nSNP_file[file];
		fclose(GenoFile);
	}

	Rprintf("%i SNPs with maf lower than %g ignored (set to 0)\n", low_AF_tot, min_AF);
	if (na_tot) Rprintf("%i out of %i missing data ignored\n", na_tot, nSNP*nIND);
	free(Geno);
	return 0;
}
