#ifndef __FOURN_H
#define __FOURN_H

/******************************************************************************
 * FOURN.H
 * C support functions for FFTs
 *
 * R. Kwan
 * August 31, 1995
 *****************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

void fourn(double data[], unsigned long nn[], int ndim, int isign);
void four1(double data[], int nn, int isign);

void fournf(float data[], unsigned long nn[], int ndim, int isign);
void four1f(float data[], int nn, int isign);

void dftgl(double data[], int nn, double dataf[], int L, int points[]);

void chirpdft(double data[], int N, double dataf[], int L, 
              double theta, double phi);

void *memswap(void *s1, void *s2, size_t n);
int  _power_of_two(unsigned int n);
void  _fftshift(void *mat, unsigned int nrows, unsigned int ncols, size_t el_size);
void _fftshift_1d(void *mat, unsigned int N, size_t el_size);
int  _next_power_of_two(unsigned int n);

#ifdef __cplusplus
}
#endif

#endif
