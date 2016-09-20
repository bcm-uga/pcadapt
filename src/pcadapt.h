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




