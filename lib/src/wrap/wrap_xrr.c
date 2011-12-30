/******************************************************************************
* Wrapper for x-ray reflectivity and reflection XSW functions
*
* Authors/modifications
* ---------------------
* Tom Trainor
*
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <xrr/xrr.h>
#include "wrap.h"

/******************************************************************************
* _xref()
* Wrapper for x-ray reflectivity
*
* Parameters
* ---------
*
* Returns
* -------
*
*
* Notes
* -----
*
******************************************************************************/
int  wrxref( int nlayer, int nelem, int nthet, double *calc_params, 
                 double *d, double *rho, double *sigma, double **comp,
                 double *elem_z, double *fp, double *fpp, double *amu, double *mu_at,
                 double *theta_arr, double *R, double *Y, 
                 double *del, double *bet, double *amu_t, double *mu_t,
                 double *Re_X, double *Im_X, double *Re_Ai, double *Im_Ai, 
                 double *Re_Ar, double *Im_Ar, double *Re_g, double *Im_g)
{
    int ret;
    ret = xref(nlayer, nelem, nthet, calc_params, d, rho, sigma, comp,
               elem_z, fp, fpp, amu, mu_at,theta_arr, R, Y, 
               del, bet, amu_t, mu_t, Re_X, Im_X, Re_Ai, Im_Ai, 
               Re_Ar, Im_Ar, Re_g, Im_g);
    return (ret);

}
/*****************************************************************************/
