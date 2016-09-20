/*
 *     convert, file ped.c
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
#include "ped.h"
#include "geno.h"
#include "io_tools.h"
#include "io_error.h"
#include "io_data_int.h"
#include "read.h"

// next_token

char* next_token(char* input_file, int i, int j)
{
    char* token = strtok(NULL, SEP);
    if (!token) {
        if (j == 0) {
            Rprintf("Error while reading individual information of %s file at line %d.\n\n", input_file, i);
        }
        else {
			Rprintf("Error while reading %s file at line %d, SNP %d.\n\n", input_file, i, j);
            error("File conversion aborted.");
        }
    }
	return token;
}

// pep2geno

void ped2geno (char *input_file, char* output_file, int *N, int *M)
{
    int *data;
	int nb;

    // number of lines and columns
	nb = nb_cols_lfmm(input_file);
	*M = (nb - 6) / 2;
    *N = nb_lines(input_file, nb);

	// allocate_memory
    data = (int *) malloc((*N)*(*M) * sizeof(int));

    // read in ancestrymap format
    read_ped(input_file, *N, *M, data);

    // write in geno format
    write_geno(output_file, *N, *M, data);

	// free memory
    free(data);
}


// ped2lfmm

void ped2lfmm (char *input_file, char* output_file, int *N, int *M)
{
    int *data;
	int nb;

    // number of lines and columns
	nb = nb_cols_lfmm(input_file);
	*M = (nb - 6) / 2;
    *N = nb_lines(input_file, nb);

	// allocate memory
    data = (int *) malloc((*N)*(*M) * sizeof(int));

    // read in ancestrymap format
    read_ped(input_file, *N, *M, data);

    // write in lfmm format
    write_data_int(output_file, *N, *M, data);

	// free memory
    free(data);
}


// read_ped

void read_ped (char* input_file, int N, int M, int* data)
{
    FILE *File = NULL;
	int i;
    int max_char_per_line = 5 * (M+6) + 20;
    char *szbuff;
	char *ref;

	// allocate memory
	szbuff = (char *) malloc(max_char_per_line * sizeof(char));
	ref = (char *) malloc(M * sizeof(char));

	// init ref
    for (i = 0; i < M; i++){
		ref[i] = '0';
    }
    
    // open input file
    File = fopen_read(input_file);

	i = 0;
	while (fgets(szbuff, max_char_per_line, File) && i < N) {
		// fill current line
		fill_line_ped (data, szbuff, M, i, input_file, File, ref);	
		i++;
	}
	// test the number of lines
	test_line(input_file, File, i, N);

    fclose(File);
	// free memory
	free(szbuff);
	free(ref);
}

// read_line

void fill_line_ped (int *data, char* szbuff, int M, int i, char* input_file, FILE *File, char *ref)
{
	char* token1, *token2;
	int tmp;
	int j;

	// read individual informations
    token1 = strtok(szbuff, SEP);
    if (!token1) {
    Rprintf("Error while reading individual informations of %s file at line %d.\n\n", input_file, i+1);
    error("File conversion aborted.");
	}
    for(j = 0; j < 5; j++){
		token1 = next_token(input_file, i+1, 0);
    }
    
	j = 0;
	// read first SNP
    token1 = strtok(NULL, SEP);
    token2 = strtok(NULL, SEP);
	while (token1 && token2 && token1[0] != EOF && token2[0] != EOF && 
		token1[0] != 10 && token2[0] != 10 && j < M) {

		// test if a token is 0,1,2,3,4, A, C, T, G
		test_token_ped(token1[0], j+1, i+1, input_file);
		test_token_ped(token2[0], j+1, i+1, input_file);

		tmp = 0;
		// process genotype
		if (ref[j] == '0') {
			if (token1[0] == '0' && token2[0] == '0') {
				tmp = 9;
			} else if (token2[0] != '0' && token1[0] == '0') {
				ref[j] = token2[0];
				tmp = 9;
			} else if (token1[0] != '0' && token2[0] == '0') {
				ref[j] = token1[0];
				tmp = 9;
			} else {
				ref[j] = token2[0];
				tmp = 1;
                if (token1[0] == ref[j]){
					tmp ++;
                }
			}
		} else {
			if(token1[0] == '0' || token2[0] == '0')
				tmp = 9;		
			else {
                if (token1[0] == ref[j]){
					tmp ++;
                }
                if (token2[0] == ref[j]){
					tmp ++;
                }
			}
		}

		data[i * M + j] = tmp;

		// read next SNP
        token1 = strtok(NULL, SEP);
        token2 = strtok(NULL, SEP);
		j ++;
	}
	// test the number of columns 
    test_column(input_file, File, j, i+1, M, token1);
}

// test_token_ped

void test_token_ped(char token, int j, int i, char* input_file)
{
	if (!(token == '0' || token == '1' || token == '2' || token == '3' || 
		token == '4' || token == 'A' || token == 'C' || token == 'T' ||
		token == 'G')) {
		Rprintf("Error: in file %s, line %d, one allele of SNP %d is '%c'"
		" and not 0, 1, 2, 3, 4, A, C, T, or G.\n\n",input_file, i, j,token);
		error("File conversion aborted.");
	}
}

void test_token_pcad(char token, int j, int i, char* input_file)
{
    if (!(token == '0' || token == '1' || token == '2' || token == '3' ||
          token == '4' || token == '9')) {
        Rprintf("Error: in file %s, line %d, one allele of SNP %d is '%c'"
                " and not 0, 1, 2, 3, 4, or 9.\n\n",input_file, i, j,token);
        error("File conversion aborted.");
    }
}

void lfmm2geno(char *input_file, char* output_file, int *N, int *M)
{
    int *data;
    int nb;
    
    // number of lines and columns
    nb = nb_cols_lfmm(input_file);
    *M = nb;
    *N = nb_lines(input_file, nb);
    
    // allocate_memory
    data = (int *) malloc((*N)*(*M) * sizeof(int));
    
    // read in ancestrymap format
    read_pcad(input_file, *N, *M, data);
    
    // write in geno format
    write_geno(output_file, *N, *M, data);
    
    // free memory
    free(data);
}

void read_pcad(char* input_file, int N, int M, int* data)
{
    FILE *File = NULL;
    int i;
    int max_char_per_line = 5 * (M) + 20;
    char *szbuff;
    char *ref;
    
    // allocate memory
    szbuff = (char *) malloc(max_char_per_line * sizeof(char));
    ref = (char *) malloc(M * sizeof(char));
    
    // init ref
    for (i = 0; i < M; i++){
        ref[i] = '0';
    }
    
    // open input file
    File = fopen(input_file,"r");
    
    i = 0;
    while (fgets(szbuff, max_char_per_line, File) && i < N) {
        // fill current line
        fill_line_pcad(data, szbuff, M, i, input_file, File, ref);
        i++;
    }
    // test the number of lines
    test_line(input_file, File, i, N);
    
    fclose(File);
    // free memory
    free(szbuff);
    free(ref);
}


void fill_line_pcad(int *data, char* szbuff, int M, int i, char* input_file, FILE *File, char *ref)
{
    char* token1;
    int tmp;
    int j;
    
    token1 = strtok(szbuff, SEP);
    j = 0;
    
    while (token1 && token1[0] != EOF && token1[0] != 10 && j < M){
        
        // test if a token is 0,1,2,3,4,9
        test_token_pcad(token1[0], j+1, i+1, input_file);
        
        tmp = 0;
        // process genotype
        
        if (token1[0] == '0'){
            tmp = 0;
        } else if (token1[0] == '1'){
            tmp = 1;
        } else if (token1[0] == '2'){
            tmp = 2;
        } else {
            tmp = 9;
        }
        
        data[i * M + j] = tmp;
        
        // read next SNP
        token1 = strtok(NULL, SEP);
        j ++;
    }
    // test the number of columns
    test_column(input_file, File, j, i+1, M, token1);
}
