/******************************************************************************
* X-ray reflectivity and reflection XSW functions
*
* Authors/modifications
* ---------------------
* Tom Trainor
* 2-06-04, T2: created
* 7-11-07, T2: major update
*
*
* Notes
* -----
* This model calculates the electric field intensity at any point
* within a multilayer structure as a function of incidence angle, and
* the reflectivity and fluorescent yield profile.
*
* We assume that the multilayer has 'nlayers' and each layer
* within the model has a homogeneous composition and density.
* Hence composition and density variations are modeled by
* discrete layers with different properties.
*
* The functions/data structures assume j = 'nlayers - 1' is the top layer,
* which should generally be modelled as air/vaccum, and j = 0 is the substrate.
*
* For e.g. a 4 layer model (nlayers=4) would have:
*
*
*               Air
*          j = nlayers -1 = 3
* -------------------------------- interface 2
*          j = nlayers -2 = 2
* -------------------------------- interface 1
*          j = nlayers -3 = 1
* -------------------------------- interface 0
*          j = nlayers -4 = 0
*             Substrate
*                 |
*                 |
*                \ /
*               -inf
* ---------------------------------
* /////////////////////////////////
* /////////////////////////////////
*
* We use the following conventions:
*
* The layer thickness are in angstroms and the
* array is indexed as:
*   d[3] = Air layer thickness (Layer 3)
*   d[2] = Layer 2 thickness
*   d[1] = Layer 1 thickness
*   d[0] = Substrate thickness (Layer 0)
*
* Note that the thickness of the Air and Substrate layers are arbitraty
* for the calculation of the field intensity and reflectivity.
* The specified layer thicknesses for these layers currently determine how
* deep in the layer (ie distance away from the interface) to calculate the
* field intensity for fluorescent yield calculations.  Note if the calc
* parameter pdepth is passed, d[0] is also ignored in the FY calcs,
* rather the depth is calculated as a multiple of the penetration depth...
*
* The densities are in g/cm^3:
*   rho[3] = Air layer density (Layer 3)
*   rho[2] = Layer 2 density
*   rho[1] = Layer 1 density
*   rho[0] = Substrate density (Layer 0)
*
* The composistions are given as the elemental stiochiometric
* coefficients, or mole fractions.  The are given as:
*   comp[k][j] = fraction of element k in layer j
*
* The element Z's are passed in the array (len = nelem)
*   elem_z[k] = Z of element k
*
* We also require that fp and fpp are passed in
* (which we assume are calculated for the incident energy)
* These arrays are of len = nelem:
*   fp[k]  = f prime (electrons/atom) of element k
*   fpp[k] = f double prime (electrons/atom) of element k
*
* The amu array (len=nelem) are the elements atomic masses:
*   amu[k] = mass of element k (g/mole)
*
* The mu_at array (len=nelem) are the atomic mass absorption coefficients
* (cm^2/g) at the fluorescence energy (ie used to calc the fluorescence attenuation)
*    mu_at[k] = mass absorption coefficient of element k (cm^2/atom)
*
* Note these are related to the atomic cross sections and scattering factors via:
*    sig_at = 2*r_e*lambda*(fpp+coh+incoh) --> cm^2/atom
*    mu_at = (Na/A)*sig_at
*
* Therefore for a given layer (density rho (g/cm^3)) and formula weight (MW):
*   mu_t =  10^-8 *(rho/MW) * Sum( A[k] * mu_at[k] * x[k])
*   att  = exp(-mu_t*thickness)
*
* The calculations also allow the inclusion of interface roughness through
* a simple Debye-Waller type model.  Note that in the above model there are
* 3 interfaces.  Therefore the 'sigma' array should have three values, each
* being the rms roughness value in angstroms:
*   sigma[0] = layer0/layer1 interface roughness
*   sigma[1] = layer1/layer2 interface roughness
*   sigma[2] = layer2/layer3 interface roughness
*
* Note that the Debye-Waller model fails for high values of
* sigma (e.g. > 20 angstroms).  For rougher interfaces, or for
* rapidly varying compsition distributions, the best way to model these
* is by generating a discretized model with many homogeneous layers
* that approximates the continuous distribution (e.g. density or composition
* profiles that follow an erf distribution)
*
* Refs
* 1. Bartels et. al., Acta Cryst (1986) A42 539-545.
* 2. de Boer, PRB (1991) V44 498-511
* 3. Krol et al, PRB (1988) V38 8579-8592
* 4. Vidal and Vincent, Applied Optics (1984) V23 1794-1801
* 5. Trainor, Templeton and Eng, J Electron Spec Rel Phenom (2006) V150, 66-85
*
* For information about scattering factors see (e.g.)
* Chantler, et al. (2005), X-Ray Form Factor,
* Attenuation and Scattering Tables (version 2.1):
* http://physics.nist.gov/ffast
* Chantler, C.T., J. Phys. Chem. Ref. Data V29, 597-1048 (2000)
* Chantler, C.T., J. Phys. Chem. Ref. Data V24, 71-643 (1995).
*
* Todo
* ----
* - Fix bug in attenuation calculations
* - Look at improving roughness calcs.  ie add more models...
* - At the top of the reflectivity function, parse out the
*   calc params into local variables so we dont have to
*   remember how they are indexed in other section of the code
* - Check convolution calculation!
* - Try out the GSL integration routine
* - We should return the penetration depth for plotting
* - Make an interface to the Icalc to generate a matrix for plotting...
* - TEST, CLEANUP, COMMENT (see matlab code, ie copy stuff here)
* - Include incident and exit beam attenuation -- transmision cell geometry
*   ie |Ai| is reduced given the cell path length
*   and we can attenuate R and FY????  Tricky geom...
* - Look at modifying the book keep so that we keep track of some
*   calculated values indexed by interface rather than layer.  e.g. r[0], t[0] could
*   be the reflection and transmission coefficients for the 1/0 interface
*   (ie the 1st interface) therefore our loops could go over interface number
*   rather than layer....
*
******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "utils.h"
#include "numfcns.h"
#include "xrr.h"
#include "complex.h"

//local function prototypes
double calc_FY(double theta, ref_model *ref, angle_calc *ang_c, layer_calc *lay_c);
int    calc_X(double theta, ref_model *ref, angle_calc *ang_c, layer_calc *lay_c);
int    calc_A(double theta, ref_model *ref, angle_calc *ang_c, layer_calc *lay_c);
double calc_I(int layer_idx, double z, ref_model *ref, angle_calc *ang_c, layer_calc *lay_c);
double calc_penetration_depth(int layer_idx, ref_model *ref, angle_calc *ang_c);
int    calc_del_bet_mu(ref_model *ref, layer_calc *lay_c);
double calc_atten(int layer_idx, double z, ref_model *ref, layer_calc *lay_c);
double calc_area(double theta, double b_vert, double b_horz, double xtal_len);
double calc_spilloff(double theta, double b_vert, double xtal_len);

/******************************************************************************
* xref()
* Master function for reflectivity and reflection xsw calculations
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
* take all raw arrays, use static copies of structs and assign pointers
* then pass along.  This is the main public function
*
* Note theta array and results
* arrays are dim nthet
* int     nthet;
* double *theta_arr;
* double *R;
* double *Y;
******************************************************************************/
int xref( int nlayer, int nelem, int nthet, double *calc_params,
          double *d, double *rho, double *sigma, double **comp,
          double *elem_z, double *fp, double *fpp, double *amu, double *mu_at,
          double *theta_arr, double *R, double *Y,
          double *del, double *bet, double *amu_t, double *mu_t,
          double *Re_X, double *Im_X, double *Re_Ai, double *Im_Ai,
          double *Re_Ar, double *Im_Ar, double *Re_g, double *Im_g)
{

    int         j, ret, fy_idx;
    double      energy, ref_scale;
    double      theta, k, lam, y, wconv;
    ref_model   ref;
    layer_calc  lay_c;
    angle_calc  ang_c;
    double complex Xtop;

    //calc k
    energy = calc_params[0];
    lam = 12398.0 / (energy);
    k = 2*M_PI/lam;

    // get fy_idx
    fy_idx = (int) calc_params[6];

    //ref scale factor
    ref_scale = calc_params[13];

    //fill in the data structure
    //ref_model
    ref.calc_params = calc_params;
    ref.k           = k;
    ref.nlayer      = nlayer;
    ref.d           = d;
    ref.rho         = rho;
    ref.comp        = comp;
    ref.sigma       = sigma;
    ref.nelem       = nelem;
    ref.elem_z      = elem_z;
    ref.fp          = fp;
    ref.fpp         = fpp;
    ref.amu         = amu;
    ref.mu_at       = mu_at;

    //layer_calc
    lay_c.nlayer    = nlayer;
    lay_c.del       = del;
    lay_c.bet       = bet;
    lay_c.mu_t      = mu_t;
    lay_c.amu_t     = amu_t;

    //angle_calc
    ang_c.nlayer    = nlayer;
    ang_c.Re_X      = Re_X;
    ang_c.Im_X      = Im_X;
    ang_c.Re_Ai     = Re_Ai;
    ang_c.Im_Ai     = Im_Ai;
    ang_c.Re_Ar     = Re_Ar;
    ang_c.Im_Ar     = Im_Ar;
    ang_c.Re_g      = Re_g;
    ang_c.Im_g      = Im_g;

    // get layer_calc data
    // this stuff doesnt depend on
    // incidence angle
    ret = calc_del_bet_mu(&ref, &lay_c);
    if (ret == FAILURE) return(FAILURE);

    // loop over theta
    for (j=0;j<nthet;j++){

        theta = theta_arr[j];

        // calc X's and compute reflectivity
        ret = calc_X(theta, &ref, &ang_c, &lay_c);
        if (ret == FAILURE) return(FAILURE);

        Xtop = ang_c.Re_X[nlayer-1] + ang_c.Im_X[nlayer-1] * I;
        R[j] = cabs(Xtop);

        if (calc_params[5] > 0.0){
            R[j] = R[j]*calc_spilloff(theta,calc_params[3],calc_params[2]);
        }
        if (ref_scale > 0.0){
            R[j] = R[j]*ref_scale;
        }

        if (fy_idx >= 0) {
            // calculate A's
            ret = calc_A(theta, &ref, &ang_c, &lay_c);
            if (ret == FAILURE) return(FAILURE);

            // calc FY
            y = calc_FY(theta, &ref, &ang_c, &lay_c);

            // mult by area
            if (calc_params[5] > 0.0){
                y = y*calc_area(theta,calc_params[3],calc_params[4],calc_params[2]);
                y = y*calc_spilloff(theta,calc_params[3],calc_params[2]);
            }
            Y[j] = y;
        }
    }

    // convolve and normalize
    wconv = calc_params[1];

    // Note convolve is not padded!
    if (wconv > 0.0) {
        convolve(theta_arr, R, nthet, wconv);
    }

    if (fy_idx >= 0){
        //convolve yield
        if (wconv > 0.0){
            convolve(theta_arr, Y, nthet, wconv);
        }
        // normalize yield
        norm_array(theta_arr, Y, nthet, calc_params[9], 1.0);
    }

    return(SUCCESS);
}

/******************************************************************************
* calc_FY()
* Calculate the fluorescent yield for a specied element at a given theta
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
* We assume here that calc_X and calc_A (and calc_del_bet_mu)
* have already been called!
*
* This computes the FY by first computing the element concentration
* in atoms/cm^3
*
* Note that ref->comp are mole fractions/stoichiometric coefficients
* not number density.  Therefore to compute number density:
*    N   = (Na/amu_t)*rho*x --> atoms/cm^3
*
* To get the number of atoms in a volume slice:
*    N*dV = N*del_z*A = atoms
*
* We can assume the A = 1cm^2 (unit area).  Therefore N*A -> atoms/cm
* and N*A*10^-8 = atoms/angstrom.  Collecting all the conversion factors:
*    N*A = (con/amu_t)*rho*x
*    con = Na*10^-8
*
* Note since FY is normalized, the application of constants
* and unit conversions is arbitrary, ie can set con = 1 if we want
*
******************************************************************************/
double calc_FY(double theta, ref_model *ref, angle_calc *ang_c, layer_calc *lay_c)
{
    int    j, fy_idx;
    double Intens, atten, del_z, pdepth;
    double d, z, znext, y, N, con, delta;

    //con   = 6.0221417899999999e15;
    con    = 1.0;
    fy_idx = (int) ref->calc_params[6];
    if (fy_idx < 0){
        return (0.0);
    }
    del_z  = ref->calc_params[11];
    pdepth = fabs(ref->calc_params[12]);

    // loop over all layer
    y = 0.0;
    for (j = 0; j < ref->nlayer; j++){

        // Compute the element concetration in layer j
        N = con*(ref->comp[fy_idx][j])*(ref->rho[j]) / lay_c->amu_t[j] ;

        //printf("**Layer=%d, Comp=%6.6f, rho=%6.6f, MW=%6.6f\n",j,
        //ref->comp[fy_idx][j],ref->rho[j],lay_c->amu_t[j]);

        // only loop through z if N > 0
        if (N > 0.0){
            d = ref->d[j];

            // for j==0 we are in the base layer.
            if (j == 0) {
                //if (d > 0.){
                if (pdepth > 0.0){
                    d = pdepth*calc_penetration_depth(0,ref,ang_c);
                }
            }
            // Note make sure allways d > 0
            // ie for j == 0 we compute into
            // the substrate, while the others
            // we go towards the top
            if (d < 0.0) d = -1.0* d;

            // loop over z within a layer
            // first segment is integral 0<->delta
            // second segment is del_z <-> 2*delta
            // etc.  the field intensity is calc
            // at the midpoint of each segment == z
            // note if del_z > d then we just do one
            // slice with thickness delta = d
            // otherwise integral vol element delta = del_z
            if (del_z>=d){
                delta = d;
            } else {
                delta = del_z;
            }
            znext = delta/2.0;
            while (znext < d){
                z = znext;
                //increment and check end of range
                znext = z+delta;
                if (znext > d){
                    z = z - 0.5*delta;
                    delta = d - z;
                    z = z + delta/2.0;
                }
                // get I(z). note if j == 0, calc_I handles
                // converting z to < 0
                Intens = calc_I(j, z, ref, ang_c, lay_c);
                // calc atten
                atten = calc_atten(j, z, ref, lay_c);
                //atten = 1.0;

                //printf("Layer=%d, N=%6.6f, d=%6.6f, z=%6.6f, delta=%6.6f, I=%6.6f,
                //atten=%6.6f\n",j,N,d,z,delta,I,atten);

                //compute integral element
                y = y + Intens*N*delta*atten;
            }
        }
    }

    return(y);
}

/******************************************************************************
* calc_X()
* This function calcs the X's and g's for a given theta.
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
* This function computes
*   X[j] = Ar[j]/Ai[j]
* where Ar[j] and Ai[j] are the field amplitdues within layer j at
* the j/j-1 interface. For the base layer X[0] = 0,
* ie there is no reflected field in the substrate. While
* X[nlayers-1] is the field intensity at the top of the multilayer,
* therefore |X[nlayers-1]|^2 = R
*
* The X's are calculated using the optical recursion formlua, starting at the bottom
* where X[0] = 0 and calc X[1], and proceding to calc of X[nlayer-1]:
*
*   X[j] = ( r[j] + X[j-1]*phase) / (1 + r[j]*X[j-1]*phase)
*
* The r[j] values are the reflection coefficients for the j/j-1 interface,
* and calculated from:
*
*   r[j] = (g[j]-g[j-1])/(g[j]+g[j-1])
*
* This function also includes an interfacial roughness term mult into
* the reflection coefficient for the j/j-1 interface:
*
*   r[j] = r[j]*exp(-1.0(q*sig[j])^2)
*
* The phase term accounts for the phase shift and attentuation of the field
* traversing the j-1 layer (ie from the j-1/j-2 interface to the j/j-1 interface)
*
*   phase = exp( -i* 2*k*d[j-1] * g[j-1])
*
* Note if j == 1 then the phase for layer j-1 is 0.0,
* ie infinitley thick base layer
*
* Note:  g[j] = sqrt( n[j]^2 - cos(theta)^2 )
* where n[j] is the layers index of refraction: n[j] = 1-del[j]-i*bet[j]
* Hence, each layer has a unique g value. We assume here
* that the del[j] and bet[j] have already been calculated
* (see calc_del_bet(), this should be called before this function).
*
*
* Note: double check how doing the phase calc below, is that the best way??
* e.g. could we use:
*    phase = gsl_complex_exp( gsl_complex_mul_real(g1,-1.0*k*d[j]));
* instead of
*    phase = gsl_complex_exp( gsl_complex_mul(phase,g1));
*
* Note we should move the g and r calc to a new function, possibly calc the
* t values as well and store these so we dont have to recalc in calc_A() (???)
*
******************************************************************************/
int calc_X(double theta, ref_model *ref, angle_calc *ang_c, layer_calc *lay_c){
    int      i, j;
    double   cos_sqr_thet, k, q, dw;
    double  *d, *sig;
    double complex g1, g2, num, den, r, X2, X1, n1, n2, phase;

    // convenience vars
    cos_sqr_thet = square(cos(theta*M_PI/180.));
    k   = ref->k;
    q   = 2.0*k*sin((theta*M_PI/180.));
    d   = ref->d;
    sig = ref->sigma;

    // start at bottom
    // ie interface 0
    ang_c->Re_X[0] = ang_c->Im_X[0] = 0.0;
    X1 = 0.0 + 0.0 * I;

    // calc g for base layer (g[0])
    // n = 1 - del -i*bet
    // g = sqrt( n^2 - cos(theta)^2 )

    n1 =  (1.0 - lay_c->del[0]) + (-1.0*lay_c->bet[0] * I);

    g1 = csqrt( cmul(n1, n1) - (cos_sqr_thet+0.0*I));

    // store g for base layer
    ang_c->Re_g[0] = creal(g1);
    ang_c->Im_g[0] = cimag(g1);

    // loop over all layers
    // note j is upper layer and i is the interface number
    // the loop starts considering the first interface (i=0)
    // which is the interface between the base layer (j=0)
    // and the first layer (j=1) ie the 1/0 interface
    // note there are nlayer-1 interfaces...
    i = 0;
    for ( j = 1; j < ref->nlayer; j++){
        // calc g for upper layer
        // n = 1 - del -i*bet
        // g2 = sqrt( n^2 - cos(theta)^2 )
        n2 = (1.0 - lay_c->del[j]) + (-1.0 * lay_c->bet[j] * I);
        g2 = csqrt( cmul(n2, n2) - (cos_sqr_thet+0.0*I));

        // store g for layer j
        ang_c->Re_g[j] = creal(g2);
        ang_c->Im_g[j] = cimag(g2);

        // calculate r for the j/j-1 interface
        // r_j = (g2-g1)/(g2+g1)
        num = csub(g2,g1);
        den = cadd(g2,g1);
        r   = cdiv(num,den);

        // calc r including roughness
        // simple dw roughness model
        // r = r*exp(-1.0(q*sig)^2)
        // sigma[i]
        if (ref->sigma[j] > 0.0){
            //dw = exp(-1.0*square(q*ref->sigma[j]));
            dw = exp(-1.0*square(q*ref->sigma[i]));
            r  = cmul(r, (dw+0*I));
        }

        //calc phase shift and attentuation for
        //field traversing the j-1 layer
        //phase = exp( -i* 2*k*d[j-1] * g1)
        // if j == 1 then the phase for layer j-1
        // is 0.0, ie infinitley thick base layer
        if (j == 1){
 	    phase = 0.0 + 0.0*I;
        } else {
            // check here
            phase = 0.0 - 2.0*k*d[j-1]*I;
            phase = cexp( cmul(phase, g1));
        }

        // calc Xj
        // X2 = ( r + X1*phase) / (1 + r*X1*phase)

        num = cadd(r, cmul(X1, phase));
        den = cadd(cmul(r, cmul(X1,phase)), 1.0+0.0*I);
        X2 = cdiv(num, den);
        if (fabs(creal(X2)) < 1e-10){
            X2 = 0.0 + cimag(X2) * I;
        }
        if (fabs(cimag(X2)) < 1e-10){
	    X2 = creal(X2) + 0.0*I;
        }
        // store values of Xj
        ang_c->Re_X[j] = creal(X2);
        ang_c->Im_X[j] = cimag(X2);

        // use top layer as bottom for next time through
        X1 = X2;
        // do the same for g's
        g1 = g2;
        //increment the interface number
        i++;
    }
    return (SUCCESS);
}

/******************************************************************************
* calc_A()
* This function calcs the A's for a given theta.
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
* This function assumes that the X[j] and g[j] values have already
* been calculated (see calc_X())
*
* Note:  X[j] = Ar[j]/Ai[j]
* And we can derive:
*
*   Ai[j] = phase*Ai[j+1]*t[j+1] / ( 1 + phase^2 * X[j] * r[j+1])
*   Ar[j] = Ai[j]*X[j]
*
* where Ar[j] and Ai[j] are the field amplitudes within layer j at
* the j/j-1 interface.
*
* For the top layer Ai[nlayer-1] = 1.  ie everything is referenced to
* unit intensity (Io=|Ai|^2) within the top layer.  Therefore Ar[0] = X[0].
*
* Given the known values for the top layer we can procede to calculate each Ai and Ar
* value starting at the top.
*
* For the base layer X[0] = 0, ie Ar[0] = 0
*
* The values of r[j+1] and t[j+1] are the reflection and transmission coefficients
* for the j+1/j interface, respectively, and calculated from:
*
*   r[j+1] = (g[j+1]-g[j])/(g[j+1]+g[j])
*   t[j+1] =  2*g[j+1]/(g[j+1] + g[j]) = 1 + r[j+1]
*
* This function also includes an interfacial roughness term mult into
* the reflection coefficient for the j+1/j interface:
*
*   r[j+1] = r[j+1]*exp(-1.0(q*sig[j+1])^2)
*
* Note we include roughness in the t calc using the below:
*   t[j+1] = 1 + r[j+1]*exp(-1.0(q*sig[j+1])^2)
*
* Therefore as r-> 0, t-> 1.  This is only done if calc_params[10] > 0.
* Otherwise t is calc explicitley from the g's
*
* The phase term accounts for the phase shift and attentuation of the field
* traversing the j layer (ie from the j+1/j interface to the j/j-1 interface),
* ie the Ai[j] and Ar[j] are the field amplitudes at the j/j-1 interface:
*
*   phase = exp( -i*k*d[j] * g[j])
*
* Note if j == 0 then the phase for layer j is 0.0, and X[0] = 0
* ie in the infinitley thick base layer we assumer the is no reflected
* field (Ar[0] = 0.0) and we are calculating the the field amplitude at the
* top of layer 0 rather than the bottom of layer 0.
*
*
* Note: double check how doing the phase calc below, is that the best way??
* e.g. could we use:
*  phase = gsl_complex_exp( gsl_complex_mul_real(g1,-1.0*k*d[j]));
* instead of
*  phase = gsl_complex_exp( gsl_complex_mul(phase,g1));
*
*
******************************************************************************/
int calc_A(double theta, ref_model *ref, angle_calc *ang_c, layer_calc *lay_c){
    int     i, j, idx;
    double  cos_sqr_thet, k, q, dw;
    double  *d, *sig;
    double complex g1, g2, num, den, r, t, Ai_2, Ar_2, Ai_1, Ar_1, X1, phase, phase_sqr;
    //gsl_complex n1, n2;

    // convenience vars
    cos_sqr_thet = square(cos(theta*M_PI/180.));
    k = ref->k;
    q   = 2.0*k*sin((theta*M_PI/180.));
    d   = ref->d;
    sig = ref->sigma;

    // start at top
    idx = ref->nlayer - 1;

    // get g for upper layer, ie air/vaccum slab
    // ie del ~ bet ~ 0.
    // n = 1 - del -i*bet ~ 1
    // g2 = sqrt( n^2 - cos(theta)^2 )

    //g2 = gsl_complex_sqrt( gsl_complex_sub_real(gsl_complex_mul(n2, n2), cos_sqr_thet) );
    g2 = ang_c->Re_g[idx] + ang_c->Im_g[idx]*I;

    // for the top (j=nlayer-1) layer, Ai =1.0
    // therefore Ar = Ai*X = X
    // store these in A arrays
    Ai_2 = 1.0 + 0.0*I;
    Ar_2 = ang_c->Re_X[idx] + ang_c->Im_X[idx]*I;
    ang_c->Re_Ai[idx] = creal(Ai_2) ;
    ang_c->Im_Ai[idx] = cimag(Ai_2) ;
    ang_c->Re_Ar[idx] = creal(Ar_2) ;
    ang_c->Im_Ar[idx] = cimag(Ar_2) ;

    // loop over all layers, start at top - 1
    // ie loop from j = nlayer-2 --> j = 0
    // note j is lower layer in the r and t calcs
    // therefore we are considering the j+1/j interface
    // i is the interface number, ie which interface
    // are the r's and t's calc for
    i = idx - 1;
    for ( j = idx - 1; j > -1 ; j--){
        // calc g and phase for lower layer
        // n = 1 - del -i*bet
        // g1 = sqrt( n^2 - cos(theta)^2 )

        //creal(n1) = 1.0 - xpar->del[j];
        //cimag(n1) = -1.0 * xpar->bet[j];
        //g1 = gsl_complex_sqrt( gsl_complex_sub_real(gsl_complex_mul(n1, n1), cos_sqr_thet) );
        g1 = ang_c->Re_g[j] + ang_c->Im_g[j] * I;

        //calc phase and attenuation
        // phase = exp( -i*k*d[j] * g1)
        // phase_sqr = (phase)^2
        if (j == 0){
   	    phase = phase_sqr = 1.0 + 0.0*I;
        } else {
            phase = 0.0 - 1.0*k*d[j]*I;
            phase = cexp(cmul(phase, g1));
            phase_sqr = 0.0 - 2.0*k*d[j]*I;
            phase_sqr = cexp(cmul(phase_sqr, g1));
        }

        // calc r for upper layer
        // r = (g2-g1)/(g2+g1)
        num = csub(g2,g1);
        den = cadd(g2,g1);
        r   = cdiv(num, den);

        // include simple dw type roughness
        // r = r*exp(-1.0(q*sig)^2)
        if (ref->sigma[j+1] > 0.0){
            //dw = exp(-1.0*square(q*ref->sigma[j+1]));
            dw = exp(-1.0*square(q*ref->sigma[i]));
            r  = cmul(r, (dw+0.0*I));
        }

        // calc t for upper layer
        // t = 2*g2/(g2 + g1)
        if (ref->calc_params[10] == 0.0){
	    num = cmul(g2, (2.0+0.0*I));
            t = cdiv(num,den);
        } else {
            // note as as r->0,  t->1
            // ie t = 1+r, so we could just calc t from r?
             t = cadd(r, 1.0+0.0*I);
        }

        // calc Ai_1 and Ar_1
        // these are field mag in the lower slab (layer j)
        // at the bottom of the slab, relative to Ai in the top layer (j=nlayer-1)
        // Ai_1 = phase*Ai_2*t / ( 1 + phase^2 * X1 * r)
        // Ar_1 = Ai_1*X1
        X1 = ang_c->Re_X[j] + ang_c->Im_X[j] * I;
        num = cmul(t, cmul(Ai_2,phase));
        den = cadd(cmul(phase_sqr,cmul(X1,r)), 1.0+0.0*I);
        Ai_1 = cdiv(num,den);
        if (fabs(creal(Ai_1)) < 1.0e-12){
	    Ai_1 = 0.0 + cimag(Ai_1);
        }
        if (fabs(cimag(Ai_1)) < 1.0e-12){
	    Ai_1 = creal(Ai_1) + 0.0*I;
        }
        Ar_1 = cmul(Ai_1,X1);

        //store results
        ang_c->Re_Ai[j] = creal(Ai_1);
        ang_c->Im_Ai[j] = cimag(Ai_1);
        ang_c->Re_Ar[j] = creal(Ar_1);
        ang_c->Im_Ar[j] = cimag(Ar_1);

        // this is the field intensity at the bottom of layer j
        // except for layer j = 0, this is the intensity at the top
        // of the base layer (bc of how phase is calc above)

        // use bottom layer as top for next time through
        Ai_2 = Ai_1;
        Ar_2 = Ar_1;
        // do the same for g's
        g2 = g1;

        // increment interface number
        i--;
    }
    return (SUCCESS);
}

/******************************************************************************
* calc_I()
* Calculate the field intensity at an arbitrary point within layer
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
* Note assume that 0 <= z <= d[j] (but we dont check that here)
* except for the base layer (j==0), in that case z <= 0
*
******************************************************************************/
double calc_I(int layer_idx, double z, ref_model *ref, angle_calc *ang_c, layer_calc *lay_c){
    double  Intens, k;
    double complex Ai, Ar, Ei, Er, g, phase_i, phase_r;

    k = ref->k;

    Ai = ang_c->Re_Ai[layer_idx] + ang_c->Im_Ai[layer_idx] * I;
    Ar = ang_c->Re_Ar[layer_idx] + ang_c->Im_Ar[layer_idx] * I;
    g  = ang_c->Re_g[layer_idx]  + ang_c->Im_g[layer_idx]  * I;

    // calculate Ei at z within layer j
    // make sure z is neg for base layer
    if (layer_idx == 0){
        if (z > 0) z = -1.0*z;
    } else {
        if (z < 0) z = -1.0*z;
    }
    phase_i = 0.0 + 1.0*k*z * I;
    phase_i = cexp( cmul(phase_i,g));
    Ei = cmul(Ai,phase_i);

    // calculate Er at z within layer j
    if (layer_idx > 0) {
        phase_r = 0.0 - 1.0*k*z * I;
        phase_r = cexp( cmul(phase_r,g));
        Er = cmul(Ar,phase_r);
    } else {
      phase_r = Er = 0.0 + 0.0 *I;
    }

    // I = gsl_complex_abs2(gsl_complex_add(Ei,Er));
    Intens = cadd(Ei, Er);
    Intens = ( creal(Intens) * creal(Intens) +
	       cimag(Intens) * cimag(Intens) ) ;
    return Intens;
}

/******************************************************************************
* calc_penetration_depth()
* Calculate the 1/e length for the given material properties
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
* Assume that k is in angstroms^-1, therefore the penetration
* depth is angstroms. Note this always well defined for the base
* layer (layer_idx = 0), but only well defined for other layers
* when the E-field is evanescent.
******************************************************************************/
double calc_penetration_depth(int layer_idx, ref_model *ref, angle_calc *ang_c){
    double  ze, k, gpp;

    k   = ref->k;
    gpp = ang_c->Im_g[layer_idx];
    ze = 1./(2*k*gpp);
    ze = fabs(ze);
    return (ze);
}

/******************************************************************************
* calc_del_bet_mu()
* Get del and bet and mu for FY.
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
* Note these calcs are independant of incidence angle.
* (ie only need to do once at the begining)
*
* del and bet are the components of the index of refraction:
*   n = 1 - del - i*bet
*   del = (r_e*lambda^2*Na*rho/(2pi*MW))*Sum(x(Z+f'))
*   bet = (r_e*lambda^2*Na*rho/(2pi*MW))*Sum(xf'')
*
* For the mu calc we assume here that mu_at[j] = atomic mass
* absorption coefficients (cm^2/g) computed at the fluorescent energy.
* These are related to the atomic cross sections via:
*    sig_at = 2*r_e*lambda*(fpp+coh+incoh) --> cm^2/atom
*    mu_at = (Na/A)*sig_at
*
* Therefore for a given layer (density rho (g/cm^3)) and formula weight (MW):
*   mu_t =  10^-8 *(rho/MW) * Sum( A[j] * mu_at[j] * x[j])
*
*     and
*
*  att = exp(-mu_t*thickness)
*
* Where thickness is in angstroms and the 10^-8 above converts
* mu from cm^-1 to angstroms^-1
*
* This algo also computes and stores (in lay_c) the molecular wieght
* (amu_t) of the given layer
*
******************************************************************************/
int calc_del_bet_mu(ref_model *ref, layer_calc *lay_c){

    double fpt, fppt, amut, zt, con, mu_t;
    double x, energy;
    double TINY = 1.0e-20;
    int j,k, calc_mu;

    energy  = ref->calc_params[0];
    calc_mu = (int) ref->calc_params[6];

    for (j=0;j<ref->nlayer;j++){
        fpt  = 0.0;
        fppt = 0.0;
        amut = 0.0;
        zt   = 0.0;
        mu_t = 0.0;

        for (k = 0; k< ref->nelem; k++){
            x    = ref->comp[k][j];
            fpt  = fpt  + ref->fp[k]  * x;
            fppt = fppt + ref->fpp[k] * x;
            amut = amut + ref->amu[k] * x;
            zt   = zt   + ref->elem_z[k] * x;
            if (calc_mu >= 0){
                mu_t = mu_t + ref->amu[k] * ref->mu_at[k] * x;
            }
        }
        // make sure amut is > 0.  ie 1/amut occurs in various places.
        // note that amut = 0 only if there is nothing in the layer
        // ie vaccum!
        if (amut <= 0.0){ amut = TINY; }

        // n = 1 - del - i*bet
        con = (415.181 * ref->rho[j]) / ( square(energy) * amut);
        lay_c->del[j]   = con * (zt + fpt);
        lay_c->bet[j]   = con * fppt;
        lay_c->amu_t[j] = amut;

        // absorption coefficient for FY
        if (calc_mu >= 0){
            con = 1.0e-8;
            con = con * ( ref->rho[j]  / lay_c->amu_t[j] );
            lay_c->mu_t[j] = con * mu_t;
        } else {
            lay_c->mu_t[j] = 0.0;
        }
    }
    return (SUCCESS);
}

/******************************************************************************
* calc_atten()
* Calculate the attenuation
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
* Here for a given z within a given layer, calc the total attenuation
*  for the ray from z to the det.
*
* Note this may be faster if had stored exp(mu_fy_j*dj) for each layer
* then just have to correct the first term for z and then sum the rest
*
* Note if layer_idx == 0, then z is the depth below the interface, so
* we just compute from z
*
******************************************************************************/
double calc_atten(int layer_idx, double z, ref_model *ref, layer_calc *lay_c){

    double atten, sin_thetdet, d;
    int j;

    d = ref->d[layer_idx];
    // make sure z and d positive
    if (d < 0.0) d = -1.0* d;
    if (z < 0.0) z = -1.0* z;

    sin_thetdet = sin((ref->calc_params[8])*M_PI/180.0);
    atten = 0.0;
    if (layer_idx == 0){
        atten = lay_c->mu_t[layer_idx] * ( z )  / sin_thetdet ;
    }else{
        atten = lay_c->mu_t[layer_idx] * ( d - z )  / sin_thetdet ;
    }
    for ( j = layer_idx+1; j < ref->nlayer-1; j++){
        d = ref->d[j];
        if (d < 0.0) d = -1.0* d;
        atten = atten + lay_c->mu_t[j] *  d  / sin_thetdet  ;
    }
    atten = exp(-1.0*atten);
    return(atten);
}

/******************************************************************************
* calc_area()
* Calculate the area of the incident beam foot-print in cm^2
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
* b_vert, b_horz and xtal_len are in mm
*
******************************************************************************/
double calc_area(double theta, double b_vert, double b_horz, double xtal_len){
    double a, sin_thet;

    if (xtal_len <= 0.) return(1.0);

    sin_thet = sin(theta*M_PI/180.0);
    if (sin_thet < 1e-10){
        a = xtal_len*b_horz;
    }else{
        a = b_vert/sin_thet;
        if ( a > xtal_len){
            a = xtal_len;
        }
        a = a*b_horz;
    }
    //convert to cm^2
    a = 0.01 * a;
    return(a);
}

/******************************************************************************
* calc_spilloff()
* Calculate the spill off correction factor.
*
* Parameters
* ---------
*
* Returns
* -------
*
* Notes
* -----
* Io' = f*Io where f is the fraction of the projected
* beam footprint that contains the sample
*
* Useage: y_calc = y_calc * f
*
******************************************************************************/
double calc_spilloff(double theta, double b_vert, double xtal_len){
    double f, sin_thet;

    if (xtal_len <= 0.) return(1.0);

    sin_thet = sin(theta*M_PI/180.0);
    f = (xtal_len/b_vert)*sin_thet;
    if (f>1.0){
        f = 1.0;
    }
    //printf("theta = %6.5f,  f is %g\n",theta, f);
    return(f);
}
/*****************************************************************************/
