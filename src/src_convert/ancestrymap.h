/*
 * @file ancestrymap.h
 *
 * @brief functions to read/write files with ancestrymap format.
 */

#ifndef ANCESTRYMAP_H
#define ANCESTRYMAP_H

#include <stdint.h>

/**
 * count the number of individuals in ancestrymap files
 *
 * @param input_file	input file in ancestrymap format 
 *
 * @return the number of individuals
 */
int nb_ind_ancestrymap(char *input_file);

/**
 * read a file of N individuals and M SNPs in ancestrymap format
 *
 * @param input_file	input file in ancestrymap format 
 * @param N		number of individuals 
 * @param M		number of SNPs
 * @param data		output data set (of size NxM)
 *
 */
void read_ancestrymap (char* input_file, int N, int M, int* data);

/**
 * read a line from an ancestrymap file
 *
 * @param szbuff	line to read
 * @param allele	variable to read allele
 * @param name		name of the SNP
 * @param i		number of the individual
 * @param j		number of the SNP
 * @param input		input file
 * @param warning       Boolean: true if a warning has already been diplayed.
 *                               false otherwise
 */
void read_line_ancestrymap(char *szbuff, int *allele, char *name, int i, int j, 
	char *input, int *warning);

/**
 * read a file of N individuals and M SNPs in ancestrymap format
 *
 * @param input_file	input file in ancestrymap format 
 * @param output_file	output file in geno format 
 * @param N		output number of individuals 
 * @param M		output number of SNPs
 */
void ancestrymap2geno (char *input_file, char* output_file, int *N, int *M);

/**
 * read a file of N individuals and M SNPs in ancestrymap format
 *
 * @param input_file	input file in ancestrymap format 
 * @param output_file	output file in lfmm format 
 * @param N		output number of individuals 
 * @param M		output number of SNPs
 */
void ancestrymap2lfmm (char *input_file, char* output_file, int *N, int *M);

#endif // ANCESTRYMAP_H

