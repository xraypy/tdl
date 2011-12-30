/******************************************************************************
* Header for reflectivity and reflection XSW functions
******************************************************************************/
#ifndef _WRAP_H
#define _WRAP_H

// try to use compiler macros to figure
// out what system we are on
// _MSC_VER should be defined by MSVS
// __GNUC__ should be defined by gcc
#ifdef _MSC_VER 
  #define WINDOWS
#else
  #ifdef __GNUC__
    #define LINUX
  #endif
#endif

// extern statements for wrap functions
#ifdef WINDOWS
  #define EXP extern __declspec(dllexport)
#else
  #define EXP extern
#endif

// prototypes
EXP int wrhello( double **x, double *y, int nr, int nc);

EXP int wrxref( int nlayer, int nelem, int nthet, double *calc_params, 
                 double *d, double *rho, double *sigma, double **comp,
                 double *elem_z, double *fp, double *fpp, double *amu, double *mu_at,
                 double *theta_arr, double *R, double *Y, 
                 double *del, double *bet, double *amu_t, double *mu_t,
                 double *Re_X, double *Im_X, double *Re_Ai, double *Im_Ai, 
                 double *Re_Ar, double *Im_Ar, double *Re_g, double *Im_g);


#endif