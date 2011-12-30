/******************************************************************************
* lmdif.h
* Public header for levenberg-marquardt code, based on 
* the minpack routines of  Garbow, Hillstrom, and More.
******************************************************************************/
#ifndef _LMDIF_H
#define _LMDIF_H

int lmfit(int (*)(int, int, double *, double *, void *), 
      int, int, double *, double *, void *, 
      double , int *);

char *lm_info(int info);

int lmdif(int (*)(int, int, double *, double *, void *),  
      int, int, double *, double *, void *,
      double, double, double, int, double, 
      double *, int, double, int, int *, double *, int,
      int *, double *,  double *, double *, double *,double *);

int fiterr(int (*)(int, int, double *, double *, void *), 
       int, int, double *, double *, void *,
       double *, double *, int *);

double sum_squares(double *, int);

#endif
/*****************************************************************************/

