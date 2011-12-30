/******************************************************************************
* Header for fit.c
*
******************************************************************************/
#ifndef _FIT_UTILS_H
#define _FIT_UTILS_H

///////////////////////////////////////////////////////////////////////////////
// definitions and structures
///////////////////////////////////////////////////////////////////////////////

#define CA_FAILURE    -1   /* calc failed     */
#define CA_SUCCESS     0   /* calc succeded   */

// ParamMap 
struct ParamMap{
    int      n_ptr;
    double **ptr;
    double  *ptr_val;
    double  *scale;
};
typedef struct ParamMap ParamMap_t;

// Function typedef
typedef int Function();

// FitData
struct FitData{
    // fit function 
    Function     *func;
    char         *fname;
    int           n_ret;  
    void         *farg;  

    // data to be fit 
    int           n_dat;         
    double       *xdat;        // len of n_dat
    double       *ydat;        // len of n_dat
    double       *yerr;        // len of n_dat

    // fit parameters
    int           n_params;
    double       *del;          //len of n_params
    double       *del_err;      //len of n_params
    double       *dmin;         //len of n_params
    double       *dmax;         //len of n_params
    int          *float_idx;    //len of n_params
    ParamMap_t   *param_map;    //len of n_params
    
    // parameters set/returned by do_fit
    int           n_evals;             
    int           n_free;
    double        chi_sqr;
    double        red_chi_sqr;
    
    // size and results of calc/fit 
    int           n_ycalc;      // calculated from n_ret and n_dat
    double       *ycalc;        // len of n_ycalc       
    double       *resid;        // len of n_ycalc      

    // flags etc...
    int           check_plims;
    int           calc_resid;
};
typedef struct FitData FitData_t;

///////////////////////////////////////////////////////////////////////////////
// function prototypes
///////////////////////////////////////////////////////////////////////////////

// Setup
FitData_t *init_fit(Function *func, int n_ret, void  *farg, char *fname);
int clear_fit(FitData_t *fit);
int clear_fit_params(FitData_t *fit);

int new_fit_param(FitData_t *fit, int param_idx );
int set_fit_param(FitData_t *fit, int p_idx, double del, double dmin, double dmax, int flt);
int link_fit_var(FitData_t *fit, int param_idx, double *var_ptr, double scale, int reset);
int accept_fit_params(FitData_t *fit);

int update_fit_data(FitData_t *fit, int n_dat, double  *xdat, double  *ydat, double  *yerr);

// do fit/calc
int do_fit( FitData_t *fit, int fit_flag);

// output
double *get_fit_stats(FitData_t *fit);
double *get_fit_ycalc(FitData_t *fit, int mk_copy);
double *get_fit_resid(FitData_t *fit, int mk_copy);
double *get_fit_params(FitData_t *fit, int mk_copy);
double *get_fit_param_errs(FitData_t *fit, int mk_copy);
void show_fit(FitData_t *fit);

#endif
/*****************************************************************************/
