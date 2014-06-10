/******************************************************************************
* Header for reflectivity and reflection XSW functions
******************************************************************************/
#ifndef _XRR_H
#define _XRR_H

// data structures
typedef struct {
    // settings
    // calc_params[0] = energy (incident energy eV)
    // calc_params[1] = width of convolution function (in degrees)
    // calc_params[2] = length of the sample (mm)
    // calc_params[3] = vert beam size (mm)
    // calc_params[4] = horz beam size (mm)
    // calc_params[5] = flag for performing area calc (1 yes, 0 no)
    // calc_params[6] = fy_idx (idx of fy element, <0 for reflectivity only)
    // calc_params[7] = fyenergy (FY energy eV)
    // calc_params[8] = detector angle (ie take off angle btwn substrate and det) (deg)
    // calc_params[9] = theta norm (theta value to use for yield normalization) (deg)
    // calc_params[10] = roughness flag (see calc_A())
    // calc_params[11] = delta z for fy integrations (ang)
    // calc_params[12] = penetration depth multiplier for fy calcs in base layer
    //                   if 0.0 then use d[0]
    // calc_params[13] = reflectivity scale factor
    double *calc_params;
    double k;

    // layer data
    // arrays are dim nlayer
    // d are the layer thicknesses (angstrom)
    // rho are the layer densities (g/cm^3)
    // comp is the mole fractions of elements
    // and is indexed comp[element_idx][layer_idx]
    int    nlayer;
    double *d;
    double *rho;
    double **comp;

    // sigmas are interface roughness
    // there are nlayer - 1 roughness 
    // values.  sigma[0] is the 0/1 
    // interface roughness etc..
    double *sigma;

    // composition stuff
    // arrays are dim nelem
    int    nelem; 
    double *elem_z;    
    double *fp; 
    double *fpp; 
    double *amu;
    double *mu_at;
} ref_model;

// calculated composition stuff
// arrays are dim nlayer
// values independant of theta
typedef struct {
    int     nlayer;
    double *del;
    double *bet;
    double *mu_t;
    double *amu_t;
} layer_calc;

// place holders for calc data
// the arrays are all dim nlayer
// values depend on theta
typedef struct {
    int     nlayer;
    double *Re_X;
    double *Im_X;
    double *Re_Ai;
    double *Im_Ai;
    double *Re_Ar;
    double *Im_Ar;
    double *Re_g;
    double *Im_g;
} angle_calc;

// function prototypes
int xref( int nlayer, int nelem, int nthet, double *calc_params, 
          double *d, double *rho, double *sigma, double **comp,
          double *elem_z, double *fp, double *fpp, double *amu, double *mu_at,
          double *theta_arr, double *R, double *Y, 
          double *del, double *bet, double *amu_t, double *mu_t,
          double *Re_X, double *Im_X, double *Re_Ai, double *Im_Ai, 
          double *Re_Ar, double *Im_Ar, double *Re_g, double *Im_g);

#endif
/*****************************************************************************/
