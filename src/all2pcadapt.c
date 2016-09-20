/*
 *     ped2geno, file: main_ped2geno.c
 *     Copyright (C) 2013 Eric Frichot
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <R.h>
#include "src_convert/io_tools.h"
#include "src_convert/io_error.h"
#include "src_convert/ped.h"
#include "src_convert/register_convert.h"
#include "src_convert/read.h"
#include "src_convert/vcf2geno.h"
#include "all2pcadapt.h"

void wrapper_converter(char **filename,int *filetype){
    char *input_name;
    int input_type;
    input_name = *filename;
    input_type = *filetype;
    converter_all(input_name,input_type);
    Rprintf("File has been sucessfully converted.\n");
}

void converter_all(char *filename, int t){
	int argc = 2;
	char *aux[2];
	aux[1] = filename;
    if (t == 0){
        ped_convert(argc,aux);
    } else if (t == 1){
        vcf_convert(argc,aux);
    } else if (t==2){
        lfmm_convert(argc,aux);
    }
}

int ped_convert(int argc, char *argv[]) {

	int M;                      // number of SNPs
	int N;                      // number of individuals
	char input_file[512];       // input file
	char output_file[512];      // output genotype file

	analyse_param_convert(argc, argv, input_file, output_file, "pcadapt");

	ped2geno(input_file, output_file, &N, &M);
    
	print_convert(N, M);

	return 0;
}

int vcf_convert(int argc, char *argv[]) {
    
    int M;                      // number of SNPs
    int N;                      // number of individuals
    int removed;                // number of removed CNVs
    char input_file[512];		// input file
    char output_file[512];		// output genotype file
    char snp_bp_file[512];		// snp file
    char removed_bp_file[512];	// removed snp file
    char *tmp;
    
    analyse_param_convert(argc, argv, input_file, output_file, "pcadapt");
    
    tmp = remove_ext(output_file,'.','/');
    strcpy(snp_bp_file, tmp);
    strcat(snp_bp_file, ".vcfsnp");
    
    strcpy(removed_bp_file, tmp);
    strcat(removed_bp_file, ".removed");
    
    vcf2geno(input_file, output_file, &N, &M, snp_bp_file, removed_bp_file, &removed);
    
    print_convert(N, M);
    
    Rprintf("For SNP info, please check %s.\n\n",snp_bp_file);
    Rprintf("%d line(s) have been removed because these are not SNPs.\n",removed);
    Rprintf("Please, check %s file, for more information.\n\n",removed_bp_file);
    
    free(tmp);
    
    return 0;
    
}

int lfmm_convert(int argc, char *argv[]) {
    
    int M;                  // number of SNPs
    int N;                  // number of individuals
    char input_file[512];	// input file
    char output_file[512];	// output genotype file

    analyse_param_convert(argc, argv, input_file, output_file, "pcadapt");

    lfmm2geno(input_file, output_file, &N, &M);
    
    print_convert(N, M);
    
    return 0;
}


