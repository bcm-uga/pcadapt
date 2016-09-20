/*
   PCAdapt random.c
   Copyright (C) 2014 Nicolas Duforet Frebourg

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <time.h>
#include "random.h"
#include "matrix.h"

// init_random

void init_random(long long *seed)
{
    unsigned long s = (unsigned long)(*seed);
    
    // dans le cas d'un code R, l'initialisation de la seed est gérée par
    // la partie R du code (avant l'appelle à la fonction en C .C(...))
    // Dans le cas du code R, il est nécessaire que seed > 0
#ifndef USING_R
    if (*seed < 0)
        s = (unsigned long) mix(clock(), time(NULL), getpid());
    
    srand(s);
#endif
    
    *seed = s;
}

// drand


//double drand()   /* uniform distribution, (0..1] */
//{
//  return (rand()+1.0)/(RAND_MAX+1.0);
//}

double drand()
{				/* uniform distribution, (0..1] */
#ifdef USING_R
    double r;
    
    r = unif_rand();
    
    return r;
#else
    return (rand() + 1.0) / (RAND_MAX + 1.0);
#endif
}

// rand_ind

int rand_int(int size) 
{
        int i;
        double r = drand();
        double sum = 0;

        for(i=0; i<size; i++) {
                sum += 1.0;
                if ( r <= sum) {
                        return i;
                }
        }
        return -1;
}

// rand_double

double rand_double(double min, double max) 
{
        return drand() * (max - min) + min;
}

// rand_matrix_double

void rand_matrix_double(double* A, int M, int N) 
{
	int i;
	for(i=0;i<N*M;i++) {
		A[i] = drand();
	}
}

// rand_normal

double rand_normal( double mean, double var) 
{
	double x;
        x = sqrt(var) * sqrt(-2 * log (drand())) * cos(2 * PI * drand()) + mean;
	return x;
}

// rand_normal_r

double rand_normal_r() 
{
        double x = sqrt(-2 * log (drand())) * cos(2 * PI * drand());
        return x;
}

// mvn_rand

void mvn_rand(double* mu, double* L, int D, double* y) {

	int i,j;
	double *x = (double *)calloc(D, sizeof(double));

	for(i=0;i<D;i++)
		x[i] = rand_normal_r();

	for(i=0;i<D;i++) {
		y[i] = mu[i];
		for(j=0;j<D;j++) {
			y[i] += L[j*D+i]*x[j]; 
		}
	}
	free(x);
}

// rand_vector

int rand_vector(double *Pi, int size) 
{
        int i;
        double r = drand();
        double sum = 0;

        for(i=0; i<size; i++) {
                sum += Pi[i];
                if ( r <= sum) {
                        return i;
                }
        }
//	printf("Warning invalid probilities:\n");
//        displayMatrix(Pi, 1, size);
        return -1;
}

// rand_gamma

double rand_gamma(int alpha, double beta) 
{
        int i = 0;
        double y = 0;
	
        for(i=0; i<alpha; i++) {
        	y += log(drand());
	}
        y *= -(1.0/beta);

	return y;
}

double rand_invgamma(int alpha, double beta){

	return 1.0/rand_gamma(alpha, 1/beta);
	
}

double rand_chisq(int n){

	int i;
	double x = 0;
	for (i=0; i<n; i++) x += pow(rand_normal_r(), 2);
	return x;

}

double rand_beta(double alpha, double beta){

	double X = rand_gamma((int) alpha, 1);
	double Y = rand_gamma((int) beta, 1);
	return X/(X + Y);

}

void rand_dirichlet(double *x, double *alpha, int K){

	int i;
	double s = 0;
	for (i=0; i<K; i++){
		x[i] = rand_gamma(alpha[i], 1);
		s += x[i];
	}
	for (i=0; i<K; i++) x[i] /= s;
}

