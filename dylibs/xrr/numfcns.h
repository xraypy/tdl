/******************************************************************************
* Header for numfcns.c
******************************************************************************/
#ifndef _MATH_UTILS_H
#define _MATH_UTILS_H

// include gsl
// for complex math etc..
#include "complex.h"

// function prototypes
double square(double a);
int    d2i(double a);
double line(double x, double offset, double slope );
double gauss(double x, double xcen, double fwhm, double mag);
double lor(double x, double xcen, double fwhm, double mag);
double trapz_int(double *x, double *y, int n);
int    convolve(double *x, double *y, int num, double wconv);
int    convolve_pad(double *x, double *y, int num, double wconv);
int    norm_array(double *x, double *y, int num, double xnorm, double norm_val);

//void print_matrix(int m, int n, double *y );

// complex fnc prototypes
typedef struct FCOMPLEX {float r,i;} fcomplex;

fcomplex Cadd(fcomplex a, fcomplex b);
fcomplex Csub(fcomplex a, fcomplex b);
fcomplex Cmul(fcomplex a, fcomplex b);
fcomplex Complex(float re, float im);
fcomplex Conjg(fcomplex z);
fcomplex Cdiv(fcomplex a, fcomplex b);
float Cabs(fcomplex z);
fcomplex Csqrt(fcomplex z);
fcomplex RCmul(float x, fcomplex a);


#endif
/*****************************************************************************/
