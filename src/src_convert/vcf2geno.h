/**
 * @file vcf2geno.h
 *
 * @brief functions to read a file in vcf format and write in geno format.
 */

#ifndef VCF2GENO_H
#define VCF2GENO_H

#include <stdint.h>

/**
 * count the number of individual in a vcf file
 *
 * @param file	file name (vcf format)
 *
 * @return number of individuals
 */
int nb_cols_vcf(char *file);

/**
 * read cnv informations for the file and and save them
 *
 * @param token		to read infos
 * @param infos		to save infos
 * @param szbuff	informations to read
 * @param j		number of the line (from 1)
 */
void read_cnv_info(char *token, char **infos, char* szbuff, int j);

/**
 * write snp info in the file
 *
 * @param output_file	output file
 * @param infos		informations to write 
 * @param removed	boolean if snp removed, write REMOVED at the end of the line
 */
void write_snp_info (FILE* output_File, char **infos, int removed);

/**
 * convert a file from vcf to geno format and store additional informatio
 *
 * @param input_file	input_file in vcf format
 * @param output_file	output file in geno format
 * @param *N		output number of individuals
 * @param *M		output number of snps
 * @param snp_bp_file	output informations about SNPs
 * @param removed_bp_file	output informations about removed SNVs
 */ 
void vcf2geno (char *input_file, char* output_file, int *N, int *M, char* snp_file, 
	char* removed_file, int *removed);

/**
 *
 */
void fill_line_vcf(char *token, int *allele, int j, int N, char* input_file, FILE* input_File);

#endif // VCF2GENO_H
