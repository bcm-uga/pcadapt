#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "src_convert/ped.h"
#include "src_convert/register_convert.h"
#include "src_convert/read.h"
#include "src_convert/vcf2geno.h"

#ifndef MAIN_PED2GENO_H_
#define MAIN_PED2GENO_H_

void wrapper_converter(char **filename, int *filetype);

void converter_all(char *filename, int t);

int ped_convert(int argc, char *argv[]);

int vcf_convert(int argc, char *argv[]);

int lfmm_convert(int argc, char *argv[]);

#endif /* MAIN_PED2GENO_H_ */
