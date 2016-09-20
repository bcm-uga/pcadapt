#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "Data__f.h"
#include "linAlgebra.h"
#include "matrix.h"

#define NA 9
#define TOL .0001

//get_row
//read a block of lines (snp), scale it and store it

int get_row(double *Geno, FILE *GenoFile, int nIND, double *mean, double *SNPSd, int sc, int blocksize, int haploid, double min_AF, int *low_AF_tot);

//add_to_cov
//Given a block of snps, computes correlation matrix of the block and add it to the existing matrix.

void add_to_cov(double *Cov, double *scratchCov, int nIND, double *Geno, int blocksize);

//Cov_line
//function to learn covariance (or cor) matrix, while reading a data file, by storing blocks of lines.
//TODO: parallelize

int Cov_line(double *Cov, double *SNPSd, int nSNP, int *nSNP_file, int nIND, int sc, char **GenoFileName, int nfile, int haploid, double min_AF);
