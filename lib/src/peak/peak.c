/******************************************************************************
* Peak fitting functions 
*
* Authors/modifications
* ---------------------
* Tom Trainor, 4-02-03
*
* Notes
* -----
* Give an example of use...
* 
* 
* Todo
* ----
* Add more background options
*
******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <utils/utils.h>
#include <utils/numfcns.h>
#include <fit/fit.h>
#include <peak/peak.h>

/******************************************************************************
* init_peak()
* Create a new peak instance
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
*
******************************************************************************/
FitData_t *init_peak(int n_peaks) 
{
    FitData_t *fit;
    Peak_t    *farg;
    int        n_par, j;

    // how many params
    if (n_peaks < 1){
        n_par = 2;
        n_peaks = 0;
    } else {
        n_par = 2 + n_peaks*4;
    }

    // make the farg structure
    farg = new_array(1,Peak_t);

    // n_peaks
    farg->n_peaks = n_peaks;

    // bgr
    farg->bgr_params = new_array(2, double);
    farg->bgr_params[0] = 0.0;
    farg->bgr_params[1] = 0.0;

    // pk_include
    farg->pk_include = new_array(n_peaks + 1, int);
    for (j=0;j<n_peaks+1;j++){
        farg->pk_include[j] = TRUE;
    }

    //  peak params
    if (n_par > 2){
        farg->pk_params = new_array(n_peaks, double *);
        for (j=0;j<n_peaks;j++){
            farg->pk_params[j] = new_array(4, double);
            farg->pk_params[j][0] = 0.0;
            farg->pk_params[j][1] = 0.0;
            farg->pk_params[j][2] = 0.0;
            farg->pk_params[j][3] = 0.0;
        }
    }

    // init the fit and add func data and farg
    fit = init_fit(calc_peak, 0, farg, "Peak");

    return(fit);
}

/******************************************************************************
* clear_peak()
* Free allocated memory
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
*
******************************************************************************/
int clear_peak(FitData_t *fit) 
{
    Peak_t *farg;
    int     j=0;

    if (string_compare(fit->fname,"Peak") != 0){
        buff_print("Error, this is not a Peak function\n");
        return(FAILURE);
    }

    // get rid of farg pointer in fit
    farg = fit->farg;
    fit->farg = NULL;

    // free the peak structure data
    free(farg->bgr_params);
    free(farg->pk_include);
    for (j=0;j<farg->n_peaks;j++){
        free(farg->pk_params[j]);
    }
    free(farg->pk_params);
    free(farg);

    // clear all the fit params too
    clear_fit(fit);
    free(fit);
    fit = NULL;

    return(SUCCESS);
}

/******************************************************************************
* set_peak_bgr()
* Set bgr params for a peak
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
*
******************************************************************************/
int set_peak_bgr(FitData_t *fit, double offset, double slope) 
{
    Peak_t *farg;

    if (string_compare(fit->fname,"Peak") != 0){
        buff_print("Error, this is not a Peak function\n");
        return(FAILURE);
    }
    farg = (Peak_t *) fit->farg;
    farg->bgr_params[0] = offset;
    farg->bgr_params[1] = slope; 

    return(SUCCESS);
}

/******************************************************************************
* set_peak()
*
* Parameters
* ---------
* - pk_idx is 1->n_peaks
*
* Returns
* -------
*
* Notes
* -----
*
******************************************************************************/
int set_peak(FitData_t *fit, int pk_idx, double cen, double fwhm, 
             double mag, double flor) 
{
    Peak_t *farg;

    if (string_compare(fit->fname,"Peak") != 0){
        buff_print("Error, this is not a Peak function\n");
        return(FAILURE);
    }

    farg = (Peak_t *) fit->farg;

    if (pk_idx < 1 || pk_idx > farg->n_peaks ){
        buff_print("Error: Peak %d is out of range\n", pk_idx);
        return(FAILURE);
    }

    farg->pk_params[pk_idx - 1][0] = cen;
    farg->pk_params[pk_idx - 1][1] = fwhm;
    farg->pk_params[pk_idx - 1][2] = mag;
    farg->pk_params[pk_idx - 1][3] = flor;

    return(SUCCESS);
}

/******************************************************************************
* set_peak_include() 
*
* Parameters
* ---------
* - pk_idx = 0 -> n_peaks
*   0 is bgr
*
* Returns
* -------
*
* Notes
* -----
*
* 
******************************************************************************/
int set_peak_include(FitData_t *fit, int pk_idx, int pk_include) 
{
    Peak_t *farg;

    if (string_compare(fit->fname,"Peak") != 0){
        buff_print("Error, this is not a Peak function\n");
        return(FAILURE);
    }

    farg = (Peak_t *) fit->farg;

    if (pk_idx < 0 || pk_idx > farg->n_peaks ){
        buff_print("Peak is out of range\n");
        return(FAILURE);
    }

    farg->pk_include[pk_idx] = pk_include;

    return(SUCCESS);
}

/******************************************************************************
* set_peak_fit_param()
*
* Parameters
* ---------
* - pk_idx = 0 is bgr
*   pk_idx = 1 is first peak etc..
* - pk_param is [slope,offset,cen,fwhm,mag,flor]
* - param_idx is which param.  < 0 just add to end of list
*
* Returns
* -------
*
* Notes
* -----
*
*
******************************************************************************/
int set_peak_fit_param(FitData_t *fit, int pk_idx, char *pk_param, 
                       int fparam_idx, double mult_fac, int reset) 
{
    int     ret = SUCCESS;
    Peak_t *farg;
    double *farg_ptr = NULL;

    if (string_compare(fit->fname,"Peak") != 0){
        buff_print("Error, this is not a Peak function\n");
        return(FAILURE);
    }

    farg = (Peak_t *) fit->farg;

    if (pk_idx < 0 || pk_idx > farg->n_peaks){
        buff_print("Peak parameter is not valid\n");
        return(FAILURE);
    }

    if (pk_idx == 0){
        if ( string_compare(pk_param,"offset") == 0){
            farg_ptr = &(farg->bgr_params[0]);
        
        } else if ( string_compare(pk_param,"slope") == 0) {
            farg_ptr = &(farg->bgr_params[1]);
        
        } else {
            buff_print("Error: %s is not valid background parameter\n", pk_param);
            return(FAILURE);
        }
    } else {
        if ( string_compare(pk_param,"cen") == 0){
            farg_ptr = &(farg->pk_params[pk_idx-1][0]);

        } else if ( string_compare(pk_param,"fwhm") == 0) {
            farg_ptr = &(farg->pk_params[pk_idx-1][1]);
        
        } else if ( string_compare(pk_param,"mag") == 0) {
            farg_ptr = &(farg->pk_params[pk_idx-1][2]);
        
        } else if ( string_compare(pk_param,"flor") == 0) {
            farg_ptr = &(farg->pk_params[pk_idx-1][3]);
        
        } else {
            buff_print("Error: %s is not a valid peak parameter\n");
            return(FAILURE);
        }
    }
    // add to params
    //ret = add_fit_param(fit, fparam_idx, farg_ptr, mult_fac, reset);
    fparam_idx = new_fit_param(fit, fparam_idx );
    ret = link_fit_var(fit, fparam_idx, farg_ptr, mult_fac, reset);
    return(ret);
}

/******************************************************************************
* fit_peak()
* Fit/calc function (pass data and do calc, fit or integrate...)
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
*
******************************************************************************/
int fit_peak(FitData_t *fit, int ndat, double *xdat, double *ydat, 
             double *yerr, int integrate_flag, int fit_flag) 
{
    int ret;

    if (string_compare(fit->fname,"Peak") != 0){
        buff_print("Error, this is not a Peak function\n");
        return(FAILURE);
    }

    // check int flag
    if (integrate_flag == TRUE){
        // switch to integrate function
        fit->func = integrate_peak;
        fit->n_ret = -1;
    } else {
        // make sure its set to Peak function
        fit->func = calc_peak;
        fit->n_ret = 0;
    }
    // NOTE update_fit_data doesnt free ycalc, so maybe
    // memory leak if call repeated when switching between 
    // integrate and non.  need to change this in the 
    // update_fit_data function.  ?
    ret = update_fit_data(fit, ndat, xdat, ydat, yerr);
    ret = do_fit( fit, fit_flag);
    return(ret);
}

/******************************************************************************
* calc_peak()
* Function for arb num of voight functions
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
*
*
******************************************************************************/
int calc_peak(double *x, int num_x, double *ycalc, int n_ycalc, Peak_t *farg)
{
    int     j, k, npeaks;
    double  b, g, l, p;
    int     *pk_include;
    double  *bgr_params;
    double **pk_params;

    // check dims
    if (num_x != n_ycalc) {
        buff_print("Array mismatch in calc_peak\n");
        return(CA_FAILURE);
    }
    if (num_x < 1 ){
        buff_print("calc_peak requires more than one point\n");
        return(CA_FAILURE);
    }

    //get pointers/data from farg
    npeaks     = farg->n_peaks;
    pk_include = farg->pk_include;
    bgr_params = farg->bgr_params;
    pk_params  = farg->pk_params; 

    // do calc
    for (j=0;j<num_x;j++){
        p = 0;
        b = 0;
        
        if ( pk_include[0] == TRUE ) {
          b =   line( x[j], bgr_params[0], bgr_params[1] );
        }
        // calc peak part
        for ( k=0; k < npeaks; k++){
            if ( pk_include[k+1] == TRUE ) {

                g = gauss( x[j], pk_params[k][0], pk_params[k][1], pk_params[k][2] );
                
                l = lor( x[j], pk_params[k][0], pk_params[k][1], pk_params[k][2] );
                
                p =  p +  pk_params[k][3] * l  +  ( 1 - pk_params[k][3] ) * g   ;
            }
        }
        // put em together
        ycalc[j] =  b + p ; 
    }
    return(CA_SUCCESS);
}

/******************************************************************************
* integrate_peak()
* Function for integrating peak function
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
* This only gives back one value => integral
*
******************************************************************************/
int integrate_peak(double *x, int num_x, double *ycalc, 
                    int n_ycalc, Peak_t *farg)
{
    int     ret;
    double  *ytmp;

    if (n_ycalc != 1) { 
        buff_print("integrate_peak only rets one value\n"); 
        return(CA_FAILURE);
    }
    ycalc[0] = 0;

    if (num_x > 0) {
        ytmp = new_array(num_x,double);
    } else {
        buff_print("integrate_peak requires more than one value\n"); 
        return(CA_FAILURE);
    }

    ret = calc_peak(x, num_x, ytmp, num_x, farg);
    if (ret == CA_FAILURE) return(CA_FAILURE);
    
    ycalc[0] = trapz_int(x, ytmp, num_x);
    
    free(ytmp);
    return(CA_SUCCESS);
}

/*****************************************************************************/
