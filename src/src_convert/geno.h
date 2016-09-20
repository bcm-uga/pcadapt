/**
 * @file geno.h
 *
 * @brief functions to read/write files with geno format.
 */

#ifndef GENO_H
#define GENO_H

/**
 * read a line of the file in geno format and fill data
 *
 * @param data  	the data set (of size NxM)
 * @param M     	the number of lines
 * @param N     	the number of columns
 * @param j    		the number of the current line (from 0)
 * @param file_data     data file name
 * @param m_file	data file
 * @param szbuff     	line to read
 * @param warning	Boolean: true if a warning has already been diplayed.
 *				 false otherwise
 */
void fill_line_geno (int* data, int M, int N, int j, char *file_data, FILE* m_File, 
	char *szbuff, int *warning);

/**
 * read the file in geno format of size (NxM) and fill data
 *
 * @param file_data     data file name
 * @param data  	the data set (of size NxM)
 * @param N     	the number of columns
 * @param M     	the number of lines
 */
void read_geno(char *input_file, int *data, int N, int M);

/**
 * write the file in geno format from data (of size NxM)
 *
 * @param output_file   output file name
 * @param N     	the number of columns
 * @param M     	the number of lines
 * @param data  	the data set (of size NxM)
 */
void write_geno(char* output_file, int N, int M, int *data);

/**
 *
 * write a line of a geno data into geno file
 *
 * @param File		geno file opened
 * @param line		geno line data
 * @param N		size of line
 */
void write_geno_line(FILE* File, int* line, int N);

#endif // GENO_H
