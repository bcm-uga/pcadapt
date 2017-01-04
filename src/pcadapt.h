#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "src_tools/Data__f.h"
#include "src_tools/linAlgebra.h"
#include "src_tools/matrix.h"

#define MAXFILE 64

void wrapper_pcadapt(char **inputfilename,int *npc, double *minmaf, int *ploidy, char **outputfilename);

void FastPCAdapt(char *filename, int K, double min_AF, int haploid, char *OUTFILE);

void compute_covariance(char **inputfilename,int *npc, double *minmaf, int *ploidy, char **outputfilename, double *result);

void get_size(char **inputfilename, int *size);

int lrfunc_(double *scores, char **inputfilename, int *num_ind, int *num_snp, int *num_pc, int *ploidy, double *minmaf);




