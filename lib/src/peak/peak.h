/******************************************************************************
* Header for peak functions
******************************************************************************/
#ifndef _PEAK_H
#define _PEAK_H

struct Peak {
    // number of peaks
    int        n_peaks;

    // array of 0/1.  0 turns off the peak
    // 1 includes it in the calculation
    int       *pk_include;

    // background params
    // bgr_params[0], lin bgr offset
    // bgr_params[1], lin bgr slope
    double    *bgr_params;

    // For each peak
    // pk_params[0][0], center
    // pk_params[0][1], fwhm
    // pk_params[0][2], mag
    // pk_params[0][3], fraction lor
    double   **pk_params;
};
typedef struct Peak Peak_t;

FitData_t *init_peak(int n_peaks);

int clear_peak(FitData_t *fit);

int set_peak_bgr(FitData_t *fit, double offset, double slope) ;

int set_peak(FitData_t *fit, int pk_idx, double cen, double fwhm, 
             double mag, double flor);

int set_peak_include(FitData_t *fit, int pk_idx, int pk_include);

int set_peak_fit_param(FitData_t *fit, int pk_idx, char *pk_param, 
                       int fparam_idx, double mult_fac, int reset); 

int fit_peak(FitData_t *fit, int ndat, double *xdat, double *ydat, 
             double *yerr, int integrate_flag, int fit_flag) ;

int calc_peak(double *x, int num_x, double *ycalc, int n_ycalc, Peak_t *farg);

int integrate_peak(double *x, int num_x, double *ycalc, int n_ycalc, Peak_t *farg);

#endif
/*****************************************************************************/
