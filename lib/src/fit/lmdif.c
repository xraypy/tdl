/******************************************************************************
* lmdif.c
* Levenberg-marquardt code
*
* Authors/modifications
* ---------------------
* This code is based on the minpack routines from:
* Argonne national laboratory. minpack project. march 1980.
* Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
*
* Error analysis code supplied by Matt Newville, University of Washington
* and Unviersity of Chicago.  
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lmdif.h"

/*****************************************************************************/

double MACHEP = 1.2e-14; /* resolution of arithmetic */
//extern double MACHEP;

double DWARF = 1.0e-31; /* smallest nonzero number */
//extern double DWARF;

char *infmsg[] = {
    "improper input parameters",
    "the relative error in the sum of squares is at most tol",
    "the relative error between variables and the solution is at most tol",
    "conditions for info = 1 and info = 2 both hold",
    "residual is orthogonal to the columns of the jacobian to machine precision",
    "number of calls to fcn has reached or exceeded 200*(n+1)",
    "tol is too small. no further reduction in the sum of squares is possible",
    "tol too small. no further improvement in approximate solution possible"
};

// Local function prototypes
void lmpar(int, double *, int, int *, double *, double *, double, 
           double *, double *, double *, double *, double *);
void qrfac(int, int, double *, int, int, int *, int, 
           double *, double  *, double *);
void qrsolv(int, double *, int, int *, double *, double *, 
            double *, double *, double *);
double enorm(int, double *);
void fdjac2(int (*)(int, int, double *, double *, void *), 
            int, int, double *, double *, void *,
            double *, int, int *, double, double *);
int gaussj(double *, int);
double dmax1(double, double);
double dmin1(double, double);

/******************************************************************************
* lm_info()
* Return string corresponding to lmdif info code
*
* Parameters
* ---------
* - info: integer index returned from lmdif
*
* Returns
* -------
* - string message
*
* Notes
* -----
*
******************************************************************************/
char *lm_info(int info) 
{
    if ( (info < 0) || (info > 7) ) return NULL;
    return infmsg[info];
}

/******************************************************************************
* lmfit()
* Main user interface to non-linear least squares optimization
* based on  subroutine lmdif1 from minpack 
* (see authorship and copyright notice below)
*
* Parameters
* ---------
* - (*fcn): pointer to a function, with argument list
*   iflag = (*fcn)(int m, int n, double *x, double *fvec, void *farg) 
*   - fcn is the user supplied function that calculates the 
*     array fvec (the residual) which is minimized with respect 
*     to the array x (the variables), where m = length of fvec, 
*     and n = length of x (number of variables).  farg is
*     a pointer to additional function arguments
*   - return iflag is an integer. if iflag < 0 this terminates the routine
*
* - m = length of fvec (length of residual array)
* - n = length of x (number of variables)
* - x = array of variables to be minimized
* - fvec = residual array
* - farg = pointer to additional function arguments (passed to fcn)
* - tol  = fitting tolerance
* - info = return code (see lmdif)
*
* Returns
* -------
* - number of function evaluations
*
* Notes
* -----
* The purpose of lmfit is to minimize the sum of the squares of
* m nonlinear functions in n variables by a modification of
* the levenberg-marquardt algorithm (lmdif). The user must provide a
* function which calculates the functions (ie the residual array).
*
******************************************************************************/
int lmfit(int (*fcn)(int m, int n, double *x, double *fvec, void *farg),
          int m, int n, double *x, double *fvec, void *farg, double tol, int *info) 
{
    double epsfcn = MACHEP;
    double xxtol  = dmax1(tol, MACHEP);
    double ftol   = xxtol;
    double xtol   = xxtol;
    double gtol   = 0.0;
    double factor = 0.1;
    int nfev      = 0;
    int ldfjac    = m;
    int mode      = 1;
    int nprint    = 0;
    double *diag, *fjac, *qtf,  *wa1, *wa2, *wa3, *wa4;
    int *ipvt;
    
    int maxfev = 200 * (n+10);
    
    diag  = (double *) calloc(n+10, sizeof(double));
    qtf   = (double *) calloc(n+10, sizeof(double));
    wa1   = (double *) calloc(n+10, sizeof(double));
    wa2   = (double *) calloc(n+10, sizeof(double));
    wa3   = (double *) calloc(n+10, sizeof(double));
    wa4   = (double *) calloc(m+10, sizeof(double));
    fjac  = (double *) calloc((n+10) * (m+10), sizeof(double));
    ipvt  = (int *)    calloc(m+10, sizeof(int));
    
    *info = 0; 
    nfev  = lmdif(fcn, m, n, x, fvec, farg, ftol, xtol, gtol, maxfev, epsfcn, 
                  diag, mode, factor, nprint, info, fjac, ldfjac, ipvt, qtf,
                  wa1, wa2, wa3, wa4);
    if (*info >= 8) *info = 4;
    
    free(diag);  free(qtf);  free(wa1);  free(wa2);  free(wa3);  free(wa4);
    free(fjac);  free(ipvt);
    return nfev;
}

/******************************************************************************
* lmdif()
* Least squares optimization
*
* Parameters
* ---------
* - fcn is a pointer to the user-supplied funtion. 
*   - iflag = fcn(m,n,x,fvec,farg)
*   - calculate the function at x (variables) and
*     return the resulting vector in fvec.
*     The values of x are optimized to minimize
*     the sum of squares of fvec.
*   - fcn should return an integer (iflag).  if the return 
*     is less than zero this routine terminates.
* 
* - m is a positive integer input variable set to the number 
*   of functions (length of fvec).
*
* - n is a positive integer input variable set to the number
*   of variables (lngth of x). n must not exceed m.
*
* - x is an array of length n. on input x must contain
*   an initial estimate of the solution vector. on output x
*   contains the final estimate of the solution vector.
*
* - fvec is an output array of length m which contains
*   the functions evaluated at the output x.
*
* - farg is a pointer to additional function arguments (passed to fcn)
*
* - ftol is a nonnegative input variable. termination
*   occurs when both the actual and predicted relative
*   reductions in the sum of squares are at most ftol.
*   therefore, ftol measures the relative error desired
*   in the sum of squares.
*
* - xtol is a nonnegative input variable. termination
*   occurs when the relative error between two consecutive
*   iterates is at most xtol. therefore, xtol measures the
*   relative error desired in the approximate solution.
*
* - gtol is a nonnegative input variable. termination
*   occurs when the cosine of the angle between fvec and
*   any column of the jacobian is at most gtol in absolute
*   value. therefore, gtol measures the orthogonality
*   desired between the function vector and the columns
*   of the jacobian.
*
* - maxfev is a positive integer input variable. termination
*   occurs when the number of calls to fcn is at least
*   maxfev by the end of an iteration.
*
* - epsfcn is an input variable used in determining a suitable
*   step length for the forward-difference approximation. this
*   approximation assumes that the relative errors in the
*   functions are of the order of epsfcn. if epsfcn is less
*   than the machine precision, it is assumed that the relative
*   errors in the functions are of the order of the machine
*   precision.
*
* - diag is an array of length n. if mode = 1 (see
*   below), diag is internally set. if mode = 2, diag
*   must contain positive entries that serve as
*   multiplicative scale factors for the variables.
*
* - mode is an integer input variable. if mode = 1, the
*   variables will be scaled internally. if mode = 2,
*   the scaling is specified by the input diag. other
*   values of mode are equivalent to mode = 1.
*
* - factor is a positive input variable used in determining the
*   initial step bound. this bound is set to the product of
*   factor and the euclidean norm of diag*x if nonzero, or else
*   to factor itself. in most cases factor should lie in the
*   interval (.1,100.). 100. is a generally recommended value.
*
* - nprint is an integer input variable that enables controlled
*   printing of iterates if it is positive. 
*
* - info is an integer output variable. if the user has
*   terminated execution, info is set to the (negative)
*   value of iflag. see description of fcn. otherwise,
*   info is set as follows.
*    - info = 0  improper input parameters.
*    - info = 1  both actual and predicted relative reductions
*      in the sum of squares are at most ftol.
*    - info = 2  relative error between two consecutive iterates
*      is at most xtol.
*    - info = 3  conditions for info = 1 and info = 2 both hold.
*    - info = 4  the cosine of the angle between fvec and any
*      column of the jacobian is at most gtol in absolute value.
*    - info = 5  number of calls to fcn has reached or exceeded maxfev.
*    - info = 6  ftol is too small. no further reduction in
*      the sum of squares is possible.
*    - info = 7  xtol is too small. no further improvement in
*      the approximate solution x is possible.
*    - info = 8  gtol is too small. fvec is orthogonal to the
*      columns of the jacobian to machine precision.
*
* - fjac is an output m by n array. the upper n by n submatrix
*   of fjac contains an upper triangular matrix r with
*   diagonal elements of nonincreasing magnitude such that
*        t     t       t
*       p *(jac *jac)*p = r *r,
*   where p is a permutation matrix and jac is the final
*   calculated jacobian. column j of p is column ipvt(j)
*   (see below) of the identity matrix. the lower trapezoidal
*   part of fjac contains information generated during
*   the computation of r.
*
* - ldfjac is a positive integer input variable not less than m
*   which specifies the leading dimension of the array fjac.
*
* - ipvt is an integer output array of length n. ipvt
*   defines a permutation matrix p such that jac*p = q*r,
*   where jac is the final calculated jacobian, q is
*   orthogonal (not stored), and r is upper triangular
*   with diagonal elements of nonincreasing magnitude.
*   column j of p is column ipvt(j) of the identity matrix.
*
* - qtf is an output array of length n which contains
*   the first n elements of the vector (q transpose)*fvec.
*
* - wa1, wa2, and wa3 are work arrays of length n.
*
* - wa4 is a work array of length m.
*
* Returns
* -------
* - number of calls to fcn.
*
* Notes
* -----
* The purpose of lmdif is to minimize the sum of the squares of
* m nonlinear functions in n variables by a modification of
* the levenberg-marquardt algorithm. The user must provide a
* function which calculates the fvec. the jacobian is
* then calculated by a forward-difference approximation.
*
* argonne national laboratory. minpack project. march 1980.
* burton s. garbow, kenneth e. hillstrom, jorge j. more
*
******************************************************************************/
int lmdif(int (*fcn)(int m, int n, double *x, double *fvec,  void *farg),
          int m, int n, double *x, double *fvec,  void *farg,
          double ftol, double xtol, double gtol, 
          int maxfev, double epsfcn, double *diag, int mode, 
          double factor, int nprint, int *info, double *fjac,
          int ldfjac, int *ipvt, double *qtf, 
          double *wa1, double *wa2, double *wa3,double *wa4) 
{ 
    int i,iflag,ij,jj,iter,j,l;
    double actred,delta,dirder,fnorm,fnorm1,gnorm;
    double par,pnorm,prered,ratio;
    double sum,temp,temp1,temp2,temp3,xnorm;
    static double one = 1.0;
    static double p1  = 0.1;
    static double p5  = 0.5;
    static double p25 = 0.25;
    static double p75 = 0.75;
    static double p0001 = 1.0e-4;
    static double zero = 0.0;
    int nfev = 1;
    delta = zero;
    xnorm = zero;
    *info = 0;
    iflag = 0;

    /* check the input parameters for errors. */
    if( (n <= 0) || (m < n) || (ldfjac < m) || (ftol < zero)
        || (xtol < zero) || (gtol < zero) || (maxfev <= 0)
        || (factor <= zero) ) {  goto L300; }
    
    if( mode == 2 ){ 
        /* scaling by diag[] */
        for( j=0; j<n; j++ ) {
            if( diag[j] <= 0.0 ) { goto L300;}
        }   
    }

    /* evaluate the function at the starting point
    *  and calculate its norm.  */
    iflag = (*fcn)(m,n,x,fvec,farg);
    if(iflag < 0) {goto L300;}
    fnorm = enorm(m,fvec);
    /*  initialize levenberg-marquardt parameter and iteration counter.*/
    par = zero;
    iter = 1;

    /*  beginning of the outer loop. */
L30:
    /* calculate the jacobian matrix. */
    iflag = 2;
    fdjac2(fcn, m,n,x,fvec,farg,fjac,ldfjac,&iflag,epsfcn,wa4);
    nfev += n;
    if(iflag < 0) { goto L300;}

    /* if requested, call fcn to enable printing of iterates.*/
    if( nprint > 0) {
        iflag = 0;
        if(((iter-1) % nprint) == 0) {
            iflag = (*fcn)(m,n,x,fvec,farg);
            if(iflag < 0) { goto L300;}
            printf( "fnorm %.15e\n", enorm(m,fvec) );
            //buff_print( "fnorm %.15e\n", enorm(m,fvec) ); 
        }
    }
    /* compute the qr factorization of the jacobian.*/
    qrfac(m,n,fjac,ldfjac,1,ipvt,n,wa1,wa2,wa3);

    /* on the first iteration and if mode is 1, scale according
    *  to the norms of the columns of the initial jacobian.*/
    if(iter == 1)   {
        if(mode != 2) {
            for( j=0; j<n; j++ ) {
                diag[j] = wa2[j];
                if( wa2[j] == zero ) {diag[j] = one;}
            }
        }
        /* on the first iteration, calculate the norm of the scaled x
        *  and initialize the step bound delta.    */
        for( j=0; j<n; j++ ) {  wa3[j] = diag[j] * x[j];}
        xnorm = enorm(n,wa3);
        delta = factor*xnorm;
        if(delta == zero) { delta = factor;}
    }
    /* form (q transpose)*fvec and store the first n components in qtf.*/
    for( i=0; i<m; i++ ) {  wa4[i] = fvec[i];}
    jj = 0;
    for( j=0; j<n; j++) {   
        temp3 = fjac[jj];
        if(temp3 != zero) {
            sum  = zero;
            ij   = jj;
            for(i=j; i<m; i++, ij++) {  sum += fjac[ij] * wa4[i]; }
            temp = -sum/temp3;
            ij   = jj;
            for( i=j; i<m; i++, ij++) { wa4[i] += fjac[ij] * temp; }
        }
        fjac[jj] = wa1[j];
        jj += m+1;  /* fjac[j+m*j] */
        qtf[j] = wa4[j];
    }
    /* compute the norm of the scaled gradient.*/
    gnorm = zero;
    if(fnorm != zero) {
        jj = 0;
        for( j=0; j<n; j++) {
            l = ipvt[j];
            if(wa2[l] != zero) {
                sum = zero;
                ij = jj;
                for( i=0; i<=j; i++, ij++) {
                    sum += fjac[ij]*(qtf[i]/fnorm);
                }
                gnorm = dmax1(gnorm,fabs(sum/wa2[l]));
            }
            jj += m;
        }
    }
    /*   test for convergence of the gradient norm. */
    if(gnorm <= gtol) {  *info = 4;}
    if( *info != 0) {goto L300;}
    /* rescale if necessary. */
    if(mode != 2) {
        for( j=0; j<n; j++ ) { diag[j] = dmax1(diag[j],wa2[j]);}
    }
    
    /* beginning of the inner loop.*/
L200:
    /*   determine the levenberg-marquardt parameter.*/
    lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,&par,wa1,wa2,wa3,wa4);
    /*   store the direction p and x + p. calculate the norm of p. */
    for( j=0; j<n; j++ ) {
        wa1[j] = -wa1[j];
        wa2[j] = x[j] + wa1[j];
        wa3[j] = diag[j]*wa1[j];
    }
    pnorm = enorm(n,wa3);
    /* on the first iteration, adjust the initial step bound. */
    if(iter == 1) {   delta = dmin1(delta,pnorm); }
    /* evaluate the function at x + p and calculate its norm. */
    iflag = (*fcn)(m,n,wa2,wa4,farg);
    nfev += 1;
    if(iflag < 0) { goto L300;}
    fnorm1 = enorm(m,wa4);
    /* compute the scaled actual reduction. */
    actred = -one;
    if( (p1*fnorm1) < fnorm) {
        temp = fnorm1/fnorm;
        actred = one - temp * temp;
    }
    /* compute the scaled predicted reduction and
    *  the scaled directional derivative.   */
    jj = 0;
    for( j=0; j<n; j++ ) {
        wa3[j] = zero;
        l = ipvt[j];
        temp = wa1[l];
        ij = jj;
        for( i=0; i<=j; i++,ij++ )  {
            wa3[i] += fjac[ij]*temp;
        }
        jj += m;
    }
    temp1  = enorm(n,wa3)/fnorm;
    temp2  = (sqrt(par)*pnorm)/fnorm;
    prered = temp1*temp1 + (temp2*temp2)/p5;
    dirder = -(temp1*temp1 + temp2*temp2);
    /* compute the ratio of the actual to the predicted
    *  reduction.*/
    ratio = zero;
    if(prered != zero) {    ratio = actred/prered; }
    /* update the step bound. */
    if(ratio <= p25) {
        temp  = (actred >= zero) ?  p5 : p5*dirder/(dirder + p5*actred);
        if( ((p1*fnorm1) >= fnorm) || (temp < p1) )  {  temp = p1;} 
        delta = temp*dmin1(delta,pnorm/p1);
        par   = par/temp;
    } else {
        if( (par == zero) || (ratio >= p75) ) {
            delta = pnorm/p5;
            par   = p5*par;
        }
    }
    /* test for successful iteration. */
    if(ratio >= p0001) {
        /* successful iteration. update x, fvec, and their norms. */
        for( j=0; j<n; j++ ) {
            x[j]   = wa2[j];
            wa2[j] = diag[j]*x[j];
        }
        for( i=0; i<m; i++) { fvec[i] = wa4[i];}
        xnorm = enorm(n,wa2);
        fnorm = fnorm1;
        iter += 1;
    }

    /*  tests for convergence.*/
    if((fabs(actred) <= ftol) && (prered <= ftol) && (p5*ratio <= one)) {*info = 1;}
    if(delta <= xtol*xnorm) {*info = 2;}
    if( (fabs(actred) <= ftol) && (prered <= ftol) && (p5*ratio <= one)
        && ( *info == 2) ) {*info = 3;}
    if( *info != 0) { goto L300;}
    /*  tests for termination and stringent tolerances.*/
    if( nfev >= maxfev) { *info = 5;}
    if((fabs(actred)<= MACHEP) && (prered<= MACHEP) && (p5*ratio <= one)) {*info = 6;}
    if(delta <= MACHEP*xnorm) {  *info = 7;}
    if(gnorm <= MACHEP) { *info = 8;}
    if( *info != 0) { goto L300;}
    /* end of the inner loop. repeat if iteration unsuccessful.*/
    if(ratio < p0001) {goto L200;}
    /* end of the outer loop.*/
    goto L30;
    
L300:
    /* termination, either normal or user imposed. */
    if(iflag < 0) { *info = iflag;}
    return nfev;
} 

/******************************************************************************
* lmpar()
*
* Parameters
* ---------
* - n is a positive integer input variable set to the order of r.
*
* - r is an n by n array. on input the full upper triangle
*   must contain the full upper triangle of the matrix r.
*   on output the full upper triangle is unaltered, and the
*   strict lower triangle contains the strict upper triangle
*   (transposed) of the upper triangular matrix s.
*
* - ldr is a positive integer input variable not less than n
*   which specifies the leading dimension of the array r.
*
* - ipvt is an integer input array of length n which defines the
*   permutation matrix p such that a*p = q*r. column j of p
*   is column ipvt(j) of the identity matrix.
*
* - diag is an input array of length n which must contain the
*   diagonal elements of the matrix d.
*
* - qtb is an input array of length n which must contain the first
*   n elements of the vector (q transpose)*b.
*
* - delta is a positive input variable which specifies an upper
*   bound on the euclidean norm of d*x.
*
* - par is a nonnegative variable. on input par contains an
*   initial estimate of the levenberg-marquardt parameter.
*   on output par contains the final estimate.
*
* - x is an output array of length n which contains the least
*   squares solution of the system a*x = b, sqrt(par)*d*x = 0,
*   for the output par.
*
* - sdiag is an output array of length n which contains the
*   diagonal elements of the upper triangular matrix s.
*
* - wa1 and wa2 are work arrays of length n.
*
*
* Notes
* -----
* Given an m by n matrix a, an n by n nonsingular diagonal
* matrix d, an m-vector b, and a positive number delta,
* the problem is to determine a value for the parameter
* par such that if x solves the system
*
*    a*x = b ,     sqrt(par)*d*x = 0 ,
*
* in the least squares sense, and dxnorm is the euclidean
* norm of d*x, then either par is zero and
*
*   (dxnorm-delta) .le. 0.1*delta ,
*
* or par is positive and
*
*    abs(dxnorm-delta) .le. 0.1*delta .
*
* This subroutine completes the solution of the problem
* if it is provided with the necessary information from the
* qr factorization, with column pivoting, of a. That is, if
* a*p = q*r, where p is a permutation matrix, q has orthogonal
* columns, and r is an upper triangular matrix with diagonal
* elements of nonincreasing magnitude, then lmpar expects
* the full upper triangle of r, the permutation matrix p,
* and the first n components of (q transpose)*b. On output
* lmpar also provides an upper triangular matrix s such that
*
*     t   t           t
*    p *(a *a + par*d*d)*p = s*s .
*
* s is employed within lmpar and may be of separate interest.
*
* Only a few iterations are generally needed for convergence
* of the algorithm. if, however, the limit of 10 iterations
* is reached, then the output par will contain the best
* value obtained so far.
*
* argonne national laboratory. minpack project. march 1980.
* burton s. garbow, kenneth e. hillstrom, jorge j. more
*
******************************************************************************/
void lmpar(int n, double *r, int ldr, int *ipvt, double *diag, 
           double *qtb, double delta, double *par, double *x, 
           double *sdiag, double *wa1, double *wa2) 
{
    int i,iter,ij,jj,j,jm1,jp1,k,l,nsing;
    double dxnorm,fp,gnorm,parc,parl,paru;
    double sum,temp;
    static double zero = 0.0;
    static double p1 = 0.1;
    static double p001 = 0.001;
    extern double DWARF;
    
    /* compute and store in x the gauss-newton direction. if the
    *  jacobian is rank-deficient, obtain a least squares solution. */
    nsing = n;
    jj = 0;
    for( j=0; j<n; j++) {
        wa1[j] = qtb[j];
        if( (r[jj] == zero) && (nsing == n)) {nsing = j;}
        if(nsing < n) { wa1[j] = zero;}
        jj += ldr+1; /* [j+ldr*j] */
    }
    
    if(nsing >= 1) {
        for( k=0; k<nsing; k++) {
            j = nsing - k - 1;
            wa1[j] = wa1[j]/r[j+ldr*j];
            temp = wa1[j];
            jm1 = j - 1;
            if(jm1 >= 0) {
                ij = ldr * j;
                for( i=0; i<=jm1; i++ )  {
                    wa1[i] -= r[ij]*temp;
                    ij += 1;
                }
            }
        }
    }
    for( j=0; j<n; j++) {
        l = ipvt[j];
        x[l] = wa1[j];
    }

    /* initialize the iteration counter.  evaluate the function at the
    * origin, and test  for acceptance of the gauss-newton direction. */
    iter = 0;
    for( j=0; j<n; j++ ) {  wa2[j] = diag[j]*x[j]; }
    dxnorm = enorm(n,wa2);
    fp = dxnorm - delta;
    if(fp <= p1*delta)  {goto L220;}

    /* if the jacobian is not rank deficient, the newton
    *  step provides a lower bound, parl, for the zero of
    *  the function. otherwise set this bound to zero. */
    parl = zero;
    if(nsing >= n) {
        for( j=0; j<n; j++) {
            l = ipvt[j];
            wa1[j] = diag[l]*(wa2[l]/dxnorm);
        }
        jj = 0;
        for( j=0; j<n; j++) {
            sum = zero;
            jm1 = j - 1;
            if(jm1 >= 0) {
                ij = jj;
                for( i=0; i<=jm1; i++, ij++) { sum += r[ij]*wa1[i]; }
            }
            wa1[j] = (wa1[j] - sum)/r[j+ldr*j];
            jj += ldr; /* [i+ldr*j] */
        }
        temp = enorm(n,wa1);
        parl = ((fp/delta)/temp)/temp;
    }
    /* calculate an upper bound, paru, for the zero of the function. */
    jj = 0;
    for( j=0; j<n; j++) {
        sum = zero;
        ij = jj;
        for( i=0; i<=j; i++, ij++) {
            sum += r[ij]*qtb[i];
        }
        l = ipvt[j];
        wa1[j] = sum/diag[l];
        jj += ldr; /* [i+ldr*j] */
    }
    gnorm = enorm(n,wa1);
    paru = gnorm/delta;
    if(paru == zero) { paru = DWARF/dmin1(delta,p1);}
    /* if the input par lies outside of the interval (parl,paru),
    *  set par to the closer endpoint.   */
    *par = dmax1( *par,parl);
    *par = dmin1( *par,paru);
    if( *par == zero) { *par = gnorm/dxnorm;}
    /* beginning of an iteration. */
L150:
    iter += 1;
    /* evaluate the function at the current value of par.*/
    if( *par == zero) {  *par = dmax1(DWARF,p001*paru);}
    temp = sqrt( *par );
    for( j=0; j<n; j++) {wa1[j] = temp*diag[j];}
    qrsolv(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2);
    for( j=0; j<n; j++ ) {wa2[j] = diag[j]*x[j];}
    dxnorm = enorm(n,wa2);
    temp = fp;
    fp = dxnorm - delta;
    /* if the function is small enough, accept the current value
    *  of par. also test for the exceptional cases where parl
    *  is zero or the number of iterations has reached 10. */
    if(!( (fabs(fp) <= p1*delta) || 
        ((parl == zero) && (fp <= temp) && (temp < zero)) || (iter == 10) )) {
        /* compute the newton correction.*/
        for( j=0; j<n; j++) {
            l = ipvt[j];
            wa1[j] = diag[l]*(wa2[l]/dxnorm);
        }
        jj = 0;
        for( j=0; j<n; j++) {
            wa1[j] = wa1[j]/sdiag[j];
            temp = wa1[j];
            jp1 = j + 1;
            if(jp1 < n) {
                ij = jp1 + jj;
                for( i=jp1; i<n; i++, ij++) {
                    wa1[i] -= r[ij]*temp;
                }
            }
            jj += ldr; /* ldr*j */
        }
        temp = enorm(n,wa1);
        parc = ((fp/delta)/temp)/temp;
        /* depending on the sign of the function, update parl or paru. */
        if(fp > zero) { parl = dmax1(parl, *par);}
        if(fp < zero) { paru = dmin1(paru, *par);}
        /* compute an improved estimate for par. */
        *par = dmax1(parl, *par + parc);
        /* end of an iteration.*/
        goto L150;
    }
L220:
    /* termination. */
    if(iter == 0) { *par = zero;}
    return;
} 

/******************************************************************************
* qrfac()
*
* Parameters
* ---------
* - m is a positive integer input variable set to the number
*   of rows of a.
*
* - n is a positive integer input variable set to the number
*   of columns of a.
*
* - a is an m by n array. on input a contains the matrix for
*   which the qr factorization is to be computed. on output
*   the strict upper trapezoidal part of a contains the strict
*   upper trapezoidal part of r, and the lower trapezoidal
*   part of a contains a factored form of q (the non-trivial
*   elements of the u vectors described above).
*
* - lda is a positive integer input variable not less than m
*   which specifies the leading dimension of the array a.
*
* - pivot is a logical input variable. if pivot is set true,
*   then column pivoting is enforced. if pivot is set false,
*   then no column pivoting is done.
*
* - ipvt is an integer output array of length lipvt. ipvt
*   defines the permutation matrix p such that a*p = q*r.
*   column j of p is column ipvt(j) of the identity matrix.
*   if pivot is false, ipvt is not referenced.
*
* - lipvt is a positive integer input variable. if pivot is false,
*   then lipvt may be as small as 1. if pivot is true, then
*   lipvt must be at least n.
*
* - rdiag is an output array of length n which contains the
*   diagonal elements of r.
*
* - acnorm is an output array of length n which contains the
*   norms of the corresponding columns of the input matrix a.
*   if this information is not needed, then acnorm can coincide
*   with rdiag.
*
* - wa is a work array of length n. if pivot is false, then wa
*   can coincide with rdiag.
*
* Notes
* -----
* This subroutine uses householder transformations with column
* pivoting (optional) to compute a qr factorization of the
* m by n matrix a. That is, qrfac determines an orthogonal
* matrix q, a permutation matrix p, and an upper trapezoidal
* matrix r with diagonal elements of nonincreasing magnitude,
* such that a*p = q*r. the householder transformation for
* column k, k = 1,2,...,min(m,n), is of the form
*
*            t
*    i - (1/u(k))*u*u
*
* where u has zeros in the first k-1 positions. the form of
* this transformation and the method of pivoting first
* appeared in the corresponding linpack subroutine.
*
* argonne national laboratory. minpack project. march 1980.
* burton s. garbow, kenneth e. hillstrom, jorge j. more
*
******************************************************************************/
void qrfac(int m, int n, double *a, int lda, int pivot, int *ipvt,
           int lipvt, double *rdiag, double  *acnorm, double *wa) 
{
    int i,ij,jj,j,jp1,k,kmax,minmn;
    double ajnorm,sum,temp;
    static double zero = 0.0;
    static double one = 1.0;
    static double p05 = 0.05;
    extern double MACHEP;

    /* compute the initial column norms and initialize several arrays. */
    ij = 0;
    for( j=0; j<n; j++) {
        acnorm[j] = enorm(m,&a[ij]);
        rdiag[j] = acnorm[j];
        wa[j] = rdiag[j];
        if(pivot != 0) { ipvt[j] = j;}
        ij += m; /* m*j */
    }
    /* reduce a to r with householder transformations. */
    minmn = ( m <= n ) ? m : n;
    for( j=0; j<minmn; j++ ){
        if(pivot != 0) {
            /* bring the column of largest norm into the pivot position. */
            kmax = j;
            for( k=j; k<n; k++) {
                if(rdiag[k] > rdiag[kmax]) { kmax = k;}
            }
            if(kmax != j) { 
                ij = m * j;
                jj = m * kmax;
                for( i=0; i<m; i++, ij++, jj++) {
                    temp = a[ij]; /* [i+m*j] */
                    a[ij] = a[jj]; /* [i+m*kmax] */
                    a[jj] = temp;
                }
                rdiag[kmax] = rdiag[j];
                wa[kmax] = wa[j];
                k = ipvt[j];
                ipvt[j] = ipvt[kmax];
                ipvt[kmax] = k;
            }
        }
        /* compute the householder transformation to reduce the
        *  j-th column of a to a multiple of the j-th unit vector. */
        jj = j + m*j;
        ajnorm = enorm(m-j,&a[jj]);
        if(ajnorm != zero) {    
            if(a[jj] < zero) { ajnorm = -ajnorm; }
            ij = jj;
            for( i=j; i<m; i++, ij++) {  a[ij] /= ajnorm; }
            a[jj] += one;
            /* apply the transformation to the remaining columns
            *  and update the norms. */
            jp1 = j + 1;
            if(jp1 < n ) {
                for( k=jp1; k<n; k++) {
                    sum = zero;
                    ij = j + m*k;
                    jj = j + m*j;
                    for( i=j; i<m; i++, ij++, jj++) { sum += a[jj]*a[ij];  }
                    temp = sum/a[j+m*j];
                    ij = j + m*k;
                    jj = j + m*j;
                    for( i=j; i<m; i++, ij++, jj++) { a[ij] -= temp*a[jj]; }
                    if( (pivot != 0) && (rdiag[k] != zero)) {
                        temp = a[j+m*k]/rdiag[k];
                        temp = dmax1( zero, one-temp*temp );
                        rdiag[k] *= sqrt(temp);
                        temp = rdiag[k]/wa[k];
                        if( (p05*temp*temp) <= MACHEP)  {
                            rdiag[k] = enorm(m-j-1,&a[jp1+m*k]);
                            wa[k] = rdiag[k];
                        }
                    }
                }
            }
        }
        rdiag[j] = -ajnorm;
    }
    return;
} 


/******************************************************************************
* qrsolv()
*
* Parameters
* ---------
* - n is a positive integer input variable set to the order of r.
*
* - r is an n by n array. on input the full upper triangle
*   must contain the full upper triangle of the matrix r.
*   on output the full upper triangle is unaltered, and the
*   strict lower triangle contains the strict upper triangle
*   (transposed) of the upper triangular matrix s.
*
* - ldr is a positive integer input variable not less than n
*   which specifies the leading dimension of the array r.
*
* - ipvt is an integer input array of length n which defines the
*   permutation matrix p such that a*p = q*r. column j of p
*   is column ipvt(j) of the identity matrix.
*
* - diag is an input array of length n which must contain the
*   diagonal elements of the matrix d.
*
* - qtb is an input array of length n which must contain the first
*   n elements of the vector (q transpose)*b.
*
* - x is an output array of length n which contains the least
*   squares solution of the system a*x = b, d*x = 0.
*
* - sdiag is an output array of length n which contains the
*   diagonal elements of the upper triangular matrix s.
*
* - wa is a work array of length n.
*
* Notes
* -----
* Given an m by n matrix a, an n by n diagonal matrix d,
* and an m-vector b, the problem is to determine an x which
* solves the system
*
*   a*x = b ,     d*x = 0 ,
*
* in the least squares sense.
*
* This subroutine completes the solution of the problem
* if it is provided with the necessary information from the
* qr factorization, with column pivoting, of a. That is, if
* a*p = q*r, where p is a permutation matrix, q has orthogonal
* columns, and r is an upper triangular matrix with diagonal
* elements of nonincreasing magnitude, then qrsolv expects
* the full upper triangle of r, the permutation matrix p,
* and the first n components of (q transpose)*b. the system
* a*x = b, d*x = 0, is then equivalent to
*
*       t       t
*    r*z = q *b ,  p *d*p*z = 0 ,
*
* where x = p*z. if this system does not have full rank,
* then a least squares solution is obtained. on output qrsolv
* also provides an upper triangular matrix s such that
*
*     t   t       t
*    p *(a *a + d*d)*p = s *s .
*
* s is computed within qrsolv and may be of separate interest.
*
* argonne national laboratory. minpack project. march 1980.
* burton s. garbow, kenneth e. hillstrom, jorge j. more
*
******************************************************************************/
void qrsolv(int n, double *r, int ldr,int *ipvt, double *diag, 
            double *qtb, double *x, double *sdiag, double *wa) 
{
    int i,ij,ik,kk,j,jp1,k,kp1,l,nsing;
    double cos,cotan,qtbpj,sin,sum,tan,temp;
    static double zero = 0.0;
    static double p25 = 0.25;
    static double p5 = 0.5;
    
    /* copy r and (q transpose)*b to preserve input and initialize s.
    *  in particular, save the diagonal elements of r in x.*/
    kk = 0;
    for( j=0; j<n; j++) {
        ij = kk;
        ik = kk;
        for( i=j; i<n; i++ , ij++) {
            r[ij] = r[ik];
            ik += ldr; /* [j+ldr*i] */
        }
        x[j] = r[kk];
        wa[j] = qtb[j];
        kk += ldr+1; /* j+ldr*j */
    }
    /* eliminate the diagonal matrix d using a givens rotation. */
    for( j=0; j<n; j++ ) {
        /* prepare the row of d to be eliminated, locating the
        *  diagonal element using p from the qr factorization.*/
        l = ipvt[j];
        if ( diag[l] != zero) { 
            for( k=j; k<n; k++ ) {sdiag[k] = zero;}
            sdiag[j] = diag[l];
            /* the transformations to eliminate the row of d
            *  modify only a single element of (q transpose)*b
            *  beyond the first n, which is initially zero.     */
            qtbpj = zero;
            for( k=j; k<n; k++) {
                /* determine a givens rotation which eliminates the
                * appropriate element in the current row of d. */
                if(sdiag[k] == zero) {   continue;  }
                kk = k + ldr * k;
                if(fabs(r[kk]) < fabs(sdiag[k])) {
                    cotan = r[kk]/sdiag[k];
                    sin = p5/sqrt(p25+p25*cotan*cotan);
                    cos = sin*cotan;
                } else {
                    tan = sdiag[k]/r[kk];
                    cos = p5/sqrt(p25+p25*tan*tan);
                    sin = cos*tan;
                }
                /* compute the modified diagonal element of r and
                *  the modified element of ((q transpose)*b,0).*/
                r[kk] = cos*r[kk] + sin*sdiag[k];
                temp  = cos*wa[k] + sin*qtbpj;
                qtbpj = -sin*wa[k] + cos*qtbpj;
                wa[k] = temp;
                /* accumulate the tranformation in the row of s. */
                kp1 = k + 1;
                if( n > kp1 ) {
                    ik = kk + 1;
                    for( i=kp1; i<n; i++, ik++ ) {
                        temp = cos*r[ik] + sin*sdiag[i];
                        sdiag[i] = -sin*r[ik] + cos*sdiag[i];
                        r[ik] = temp;
                    }
                }
            }
        }
        /* store the diagonal element of s and restore
        *  the corresponding diagonal element of r. */
        kk = j + ldr*j;
        sdiag[j] = r[kk];
        r[kk] = x[j];
    }
    /* solve the triangular system for z. if the system is
    *  singular, then obtain a least squares solution.*/
    nsing = n;
    for( j=0; j<n; j++) {
        if( (sdiag[j] == zero) && (nsing == n) ) {nsing = j;}
        if(nsing < n) {wa[j] = zero;}
    }
    if(nsing >= 1) {
        for( k=0; k<nsing; k++ ) {
            j = nsing - k - 1;
            sum = zero;
            jp1 = j + 1;
            if(nsing > jp1) {
                ij = jp1 + ldr * j;
                for( i=jp1; i<nsing; i++ , ij++) {
                    sum += r[ij]*wa[i];
                }
            }
            wa[j] = (wa[j] - sum)/sdiag[j];
        }
    }
    /* permute the components of z back to components of x. */
    for( j=0; j<n; j++ ) {
        l = ipvt[j];
        x[l] = wa[j];
    }
    return;
}

/******************************************************************************
* enorm()
* Euclidina norm of a vector
*
* Parameters
* ---------
* - n is a positive integer input variable.
*
* - x is an input array of length n.
*
* Notes
* -----
* Given an n-vector x, this function calculates the
* euclidean norm of x.
*
* The euclidean norm is computed by accumulating the sum of
* squares in three different sums. The sums of squares for the
* small and large components are scaled so that no overflows
* occur. Non-destructive underflows are permitted. Underflows
* and overflows do not occur in the computation of the unscaled
* sum of squares for the intermediate components.
* The definitions of small, intermediate and large components
* depend on two constants, rdwarf and rgiant. The main
* restrictions on these constants are that rdwarf**2 not
* underflow and rgiant**2 not overflow. The constants
* given here are suitable for every known computer.
*
* argonne national laboratory. minpack project. march 1980.
* burton s. garbow, kenneth e. hillstrom, jorge j. more
*
******************************************************************************/
double enorm(int n, double *x) 
{
    int i;
    double agiant, floatn, s1, s2, s3, xabs, x1max, x3max, ans, temp;
    double GIANT;
    extern double DWARF;

    GIANT = 1./DWARF;
    s1    = 0;  s2    = 0;  s3    = 0.;
    x1max = 0;  x3max = 0;
    floatn = n;
    agiant = GIANT/floatn;
    for(i=0; i<n; i++) { 
        /* sum components */
        xabs = fabs(x[i]);
        if( (xabs > DWARF) && (xabs < GIANT))  {  /* intermediate components.*/
            s2 += xabs*xabs;
            continue;
        }
        if(xabs > DWARF) {  
            /* large components. */
            if(xabs > x1max) {
                temp  = x1max/xabs;
                s1    = 1 + s1*temp*temp;
                x1max = xabs;
            } else {
                temp = xabs/x1max;
                s1 += temp*temp;
            }
            continue;
        }
        if(xabs > x3max) {   
            /* small components. */
            temp  = x3max/xabs;
            s3    = 1 + s3*temp*temp;
            x3max = xabs;
        } else {
            if(xabs != 0) {
                temp = xabs/x3max;
                s3 += temp*temp;
            }
        }
    }
    /* calculation of norm. */
    if(s1 != 0) {
        temp = s1 + (s2/x1max)/x1max;
        ans = x1max*sqrt(temp);
        return(ans);
    }
    if(s2 != 0) {
        if(s2 >= x3max) { 
            temp = s2*(1+(x3max/s2)*(x3max*s3));
        } else {
            temp = x3max*((s2/x3max)+(x3max*s3));
        } 
        ans = sqrt(temp);
    } else {
        ans = x3max*sqrt(s3);
    }
    return(ans);
} 

/******************************************************************************
* fdjac2()
*
* Parameters
* ---------
* - fcn is the name of the user-supplied funtion. 
*   - iflag = fcn(m,n,x,fvec,farg)
*   - calculate the function at x (variables) and
*     return the resulting vector in fvec.
*     The values of x are optimized to minimize
*     the sum of squares of fvec.
*   - fcn should return an integer (iflag).  if the return 
*     is less than zero this routine terminates.
*
* - m is a positive integer input variable set to the number 
*   of functions (length of fvec).
*
* - n is a positive integer input variable set to the number
*   of variables (lngth of x). n must not exceed m.
*
* - x is an array of length n. on input x must contain
*   an initial estimate of the solution vector. on output x
*   contains the final estimate of the solution vector.
*
* - fvec is an output array of length m which contains
*   the functions evaluated at the output x.
*
* - fjac is an output m by n array which contains the
*   approximation to the jacobian matrix evaluated at x.
*
* - ldfjac is a positive integer input variable not less than m
*   which specifies the leading dimension of the array fjac.
*
* - iflag is an integer variable which can be used to terminate
*   the execution of fdjac2. see description of fcn.
*
* - epsfcn is an input variable used in determining a suitable
*   step length for the forward-difference approximation. this
*   approximation assumes that the relative errors in the
*   functions are of the order of epsfcn. if epsfcn is less
*   than the machine precision, it is assumed that the relative
*   errors in the functions are of the order of the machine
*   precision.
*
* - wa is a work array of length m.
*
* Notes
* -----
* This subroutine computes a forward-difference approximation
* to the m by n jacobian matrix associated with a specified
* problem of m functions in n variables.
*
* argonne national laboratory. minpack project. march 1980.
* burton s. garbow, kenneth e. hillstrom, jorge j. more
*
******************************************************************************/
void fdjac2(int (*fcn)(int m, int n, double *x, double *fvec, void *farg), 
            int m, int n, double *x, double *fvec,  void *farg, double *fjac,
            int ldfjac,int *iflag, double epsfcn, double *wa) 
{
    int i,j,k;
    double eps,h,temp;
    extern double MACHEP;
    
    temp = dmax1(epsfcn,MACHEP);
    eps  = sqrt(temp);
    k    = 0;
    for(j=0; j<n; j++) {
        temp = x[j];
        h    = eps * fabs(temp);
        if(h == 0.) { h = eps;}
        x[j] = temp + h;
        *iflag = (*fcn)(m,n,x,wa,farg);
        if( *iflag < 0) return;
        x[j] = temp;
        for(i=0; i<m; i++, k++ ) {
            fjac[k] = (wa[i] - fvec[i])/h;
        }
    }
    return;
}

/******************************************************************************
* fiterr()
* Error analysis for a fit using the minpack routines 
*
* Parameters
* ---------
* - fcn is a pointer to the user-supplied funtion. 
*   - iflag = fcn(nfit, nvar, x, ftemp, farg)
*   - calculate the function at x (variables) and
*     return the resulting vector in ftemp.
*     The values of x are are the optimized variables
*     that produced the min of the sum of squares of ftemp.
*   - nfit is the length of ftemp
*   - nvar is the length of x
*   - farg is a pointer to additional function arguments
*   - fcn should return an integer (iflag).  if the return 
*     is less than zero this routine terminates.
* - nfit is number of function evaluations (length of ftemp and fbest) 
* - nvar is the number of variables (length of x)
* - fbest is the array of fit residual for best fit 
* - x is the array of best fit values for variables 
* - delta is the array of uncertainties for the variables   (computed here) 
* - correl is the array of two-variable correlations  (nvar,nvar)  (computed here) 
* - iflag is an integer array (length nvar) whose elements are 1 if the  
*   corresponding variable is suspected of causing the failure of the inversion 
*   of the curvature matrix. these may be null variables. 
*
* Return
* ------
* - 0 if succesful
*   1 if the matrix could not be inverted
*   <0 is an error code from fcn
*
* Notes
* -----
* Given a subroutine, *fcn*, to generate a fitting function 
* with *nfit* evaluations from a set of *nvar* variables, 
* with best-fit values *x* and residuals *fbest* determined, 
* this will return the uncertainties in *x* in *delta*, and 
* the correlations between the variables in *correl*.    
*
* The algorithm here is to construct and invert the nvar x nvar 
* curvature matrix  alpha, whose elements are found by summing 
* over the elements of the jacobian matrix, fjac: 
*
*    fjac(i,j) = dfvect(i) / dx(j)   (i = 1, nfit; j = 1, nvar) 
*
* where fvect is the residual array for the fit and dx is a small 
* change in one variable away from the best-fit solution. then 
*
*    alpha(j,k) = alpha(k,j) 
*               = sum_(i=1)^nfit (fjac(i,j) * fjac(i,k)) 
*
* The inverse of alpha gives the curvature matrix, whose diagonal 
* elements are used as the uncertainties and whose off-diagonal 
* elements give the correlations between variables. 
*
* Copyright 1994 university of washington  matt newville 
*
******************************************************************************/            
int fiterr(int (*fcn)(int nfit, int nvar, double *x, double *ftemp, void *farg), 
           int nfit, int nvar, double *x,  double *fbest, void *farg, 
           double *delta,  double *correl, int *iflag) 
{
    int i, j, k, ier;
    double tempx, eps, sum, delx, d;
    static double tiny = 1.e-10;
    double *ftemp, *fjac, *alpha;
    extern int gaussj(double *, int);
      
    ftemp = (double *) calloc(nfit,      sizeof(double));
    fjac  = (double *) calloc(nfit*nvar, sizeof(double));
    alpha = (double *) calloc(nvar*nvar, sizeof(double));

    /* construct jacobian using the best possible guess for the */
    /* relative error in each variable to evaluate the derivatives. */
    /* if not available, use 1% of the value for the variable. */
    eps   = .01;
    for (j = 0; j < nvar; j++) {
        tempx = x[j];
        d     = eps * fabs(tempx);
        delx  = dmax1(tiny,d);
        x[j]  = tempx + delx;
        ier   = (*fcn)(nfit, nvar, x, ftemp, farg);
        x[j]  = tempx;
        if (ier < 0) {
            free(ftemp);
            free(fjac);
            free(alpha);
            return(ier);
        }
        for (i = 0; i < nfit; i++) {
            fjac[i + j * nfit] = (fbest[i] - ftemp[i]) / delx;
        }
    }
    /* re-evaluate best-fit to restore any common block stuff */
    ier = (*fcn)(nfit, nvar, x, ftemp, farg);

    /* collect the symmetric curvature matrix, store in alpha */
    for (j = 0; j < nvar; j++) {
        for (k = 0; k <= j; k++) {
            sum = 0.;
            for (i = 0; i < nfit; i++) {
                sum += fjac[i + j * nfit] * fjac[i + k * nfit];
            }
            alpha[j + k * nvar] = sum;
            if (k != j) { alpha[k + j * nvar] = sum; }
            /*buff_print( "#==  %i  %i :  %15.7g\n",  j, k, sum);      */
        }
    }
    
    /* invert alpha to give the covariance matrix.  gaussj does */
    /* gauss-jordan elimination which will die if the matrix is */
    /* singular. although more efficient versions of this method */
    /* exist, in the event of a singular matrix, this one will */
    /* preserve the original matrix and set ier to 1. */
    ier = gaussj(alpha, nvar);
    
    /* if alpha could not be inverted, flag those variables with */
    /* small diagonal components of alpha - these are the likely */
    /* null variables that caused the matrix inversion to fail. */
    if (ier >= 1) {
        for (i = 0; i < nvar; i++) {
            iflag[i] = 0;
            if ((fabs( alpha[i + i * nvar])) <= tiny) { iflag[i] = 1; }
        }
        free(ftemp);
        free(fjac);
        free(alpha);
        return(1);
    }
    
    /* alpha now contains the covariance matrix, and is easily */
    /* converted into delta, the uncertainty for each variable, */
    /* and correl, the two-variable correlation matrix. */
    for (i = 0; i < nvar; i++) {
        delta[i] = dmax1(tiny,sqrt(fabs(alpha[i + i * nvar])));        
        for (j = 0; j <= i; j++) {
            correl[j + i * nvar] =  alpha[j + i * nvar]/(delta[i] * delta[j]);
            correl[i + j * nvar] =  correl[j + i * nvar];
        }
    }
    //buff_print( "#==  %i  %i : %15.7g %15.7g\n",i, j, correl[i+j*nvar], alpha[i+j*nvar]);
    free(ftemp);
    free(fjac);
    free(alpha);
    return (0);
}

/******************************************************************************
* gaussj()
* Gauss-jordan elimination to invert a matrix. 
*
* Parameters
* ---------
* - a    matrix to invert / solution on output  
* - na   number of elements in a to use  
* 
* Return
* ------
* - 0 on success / 1  on error
*
* Notes
* -----
* If the matrix cannot be inverted, a contains garbage 
*
* Copyright (c) 1998 matt newville 
*
******************************************************************************/
int gaussj(double *a, int na) 
{
    int  icol, irow, i, j, k;
    int  *col, *row, *pivot;
    double piv, amax, tmp;
    
    col   = (int *) calloc(na,  sizeof(int));
    row   = (int *) calloc(na,  sizeof(int));
    pivot = (int *) calloc(na,  sizeof(int));
    irow   = 0;  icol   = 0;

    /* main loop over the columns to be reduced */
    for (i = 0; i < na; i++) {
        amax = 0.;
        /* find biggest element to use as pivot */
        for (j = 0; j < na; j++) {
            if (pivot[j] != 1) {
                for (k = 0; k < na; k++) {
                    if (pivot[k] == 0) {
                        if (fabs(a[j + k * na]) >= amax) {
                            amax = fabs(a[j + k * na]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }
        pivot[icol]++;
        row[i] = irow;
        col[i] = icol;
        /* use this pivot to switch row / column elements */
        if (irow != icol) {
            for (j = 0; j < na; j++) {
                tmp              = a[irow + j * na];
                a[irow + j * na] = a[icol + j * na];
                a[icol + j * na] = tmp;
            }
        }
        /* divide the pivot row by the pivot element */
        if (a[icol + icol*na] == 0.) {
            free(row);  free(col);  free(pivot);
            return (1);
        }
        piv               = 1. / a[icol + icol*na];
        a[icol + icol*na] = 1.;
        for (j = 0; j < na; j++) { a[icol + j * na] *= piv;}
        /* reduce the rows except for the pivot one */
        for (j = 0; j < na; j++) {
            if (j != icol) {
                tmp            = a[j + icol*na];
                a[j + icol*na] = 0.;
                for (k = 0; k < na; k++) {
                    a[j + k*na] -= a[icol + k*na]*tmp;
                }
            }
        }
    }
    /* unpack: interchange column pairs in the reverse */
    /* order of the above  permutation */
    for (i = na-1; i >= 0; i--) {
        if (row[i] != col[i]) {
            for (j = 0; j < na; ++j) {
                tmp                = a[j + row[i] * na];
                a[j + row[i] * na] = a[j + col[i] * na];
                a[j + col[i] * na] = tmp;
            }
        }
    }
    free(row);  free(col);  free(pivot);
    return (0);
}

/******************************************************************************
* sum_squares()
* Returns the sum of squares of an array 
******************************************************************************/
double sum_squares (double *a, int n) {
    double s = 0.;
    int i ;
    for (i = 0 ; i < n ; i++ ) { s +=  a[i]*a[i]; }
    return s;
}

/******************************************************************************
* dmax1(), dmin1()
* Find max or min of two values
******************************************************************************/
double dmax1(double a, double b) { return (a >= b ) ? a : b; }

double dmin1(double a, double b) { return (a <= b ) ? a : b; }

/*****************************************************************************/