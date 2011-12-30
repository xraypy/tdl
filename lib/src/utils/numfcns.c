/******************************************************************************
* Math utility functions 
*
* Authors/modifications
* ---------------------
* T. Trainor, 4-06-03
*
*
* Notes
* -----
* The complex functions were originally taken from the public
* domain nurmerical recipies stite:
* http://www.nr.com/pubdom/complex.c.txt
*
* Todo
* ----
* - Put complex fcns in here to remove gsl dependency
* - Work on convolutions (add fft method)
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <utils/utils.h>
#include <utils/numfcns.h>

/******************************************************************************
* square()
* Calc the square of a number
* 
* Parameters
* ---------
* - a single double value
* 
* Returns
* -------
* - the square of the number
*
* Notes
* -----
*
******************************************************************************/
double square(double a){
    return(a*a);
}

/******************************************************************************
* d2i
* Double to integer conversion
*
* Parameters
* ---------
* - a single double value
* 
* Returns
* -------
* - corresponding integer value
*
* Notes
* -----
* This does not round, it just returns the integer part.  This is 
* essentially the same result as using the type cast: i = (int) a;
******************************************************************************/
int d2i(double a){
    int result;
    //result = (int) a;
    if (a > 0.0){
       result = (int) floor(a);
    } else {
       result = (int) ceil(a);
    }
    return(result);
}

/******************************************************************************
* line()
* Calc point on a line
* 
* Parameters
* ---------
* - x
* - offset
* - slope
* 
* Returns
* -------
* - y
*
* Notes
* -----
* Returns: y = offset + slope * x
*
******************************************************************************/
double line(double x, double offset, double slope ){
    double  y;
    y =   slope * x + offset ;
    return(y);
}

/******************************************************************************
* gauss()
* Calc point on a gaussian distribution 
*
* Parameters
* ---------
* - x
* - xcen
* - fwhm
* - mag
* 
* Returns
* -------
* - y
*
* Notes
* -----
* Returns: y = mag*exp( - (x-xcen)^2 / sig^2)
* where the standard deviation of a normal 
* distribution sig = 0.600561 * fwhm 
* 
* Note fwhm = 0.0 is safe, returns 1 for x = xcen (within 
* a small differential)
******************************************************************************/
double gauss(double x, double xcen, double fwhm, double mag) {

    double a, y;
    // note this makes it safe 
    // and makes it perform like a 
    // Kroniker-delta function
    if (fwhm == 0.0){
        if ( fabs(x - xcen) < 1e-13 ){
            return (1.0);
        } else {
            return(0.0);
        }
    }
    a =   ( x - xcen )  / ( 0.600561 * fwhm ) ;        
    y =   mag *  exp( -1 * square(a) ) ;
    return(y);
}

/******************************************************************************
* lor()
* Calc point on a lorentzian distr
* 
* Parameters
* ---------
* - x
* - xcen
* - fwhm
* - mag
* 
* Returns
* -------
* - y
*
* Notes
* -----
* Returns: y = mag/ ( 1 + (x-xcen)^2 / sig^2)
* where the standard deviation of a Lorentzian 
* distribution sig = 0.5 * fwhm 
*
******************************************************************************/
double lor(double x, double xcen, double fwhm, double mag){

    double a, y;
    if (fwhm == 0.0) return(0);
    a =   ( x - xcen )  / ( 0.5 * fwhm ) ;
    y = mag / ( 1 + square( a ) ) ;
    return(y);
}

/******************************************************************************
* trapz_int()
* Trapezoidal integration
*
* Parameters
* ---------
* - x array of length n
* - y array of length n
* - n
*
* Returns
* -------
* - area
*
* Notes
* -----
* The arrays must have at least 2 points to compute the integral
******************************************************************************/
double trapz_int(double *x, double *y, int n){
    int j;
    double a=0;
    
    if (n < 2 ) return(0);
    a = 0.0;
    for (j=0; j < n-2; j++){
        a = a +  0.5 * ( y[j+1] + y[j] )  / ( x[j+1] - x[j]  ) ;
    }    
    //buff_print("a=%f\n",a);    
    return(a);
}

/******************************************************************************
* convolve()
* Simple gaussian convolution (unpadded)
*
* Parameters
* ---------
* - x array of length num
* - y array of length num
* - wconv is the convolution width (gaussian fwhm)
*
* Returns
* -------
* - SUCCESS/FAILURE
*
* Notes
* -----
* Note y is the input 'data' array and is recomputed to 
* hold the convolution.  ie this overwrites the y data
*
* Note this has end point problems since the arrays are not padded.
*
******************************************************************************/
int convolve(double *x, double *y, int num, double wconv){

    double  *yp, *temp, xcen, n, g;
    int j, k;

    yp = new_array(num, double);
    temp = new_array(num, double);

    for (j=0;j<num;j++){
        xcen = x[j];
        n = 0.0;
        // compute gaussian scaled amplitude
        // of each point.  n is the gaussian 
        // weight used later to normalize
        for (k=0;k<num;k++){
            g     = gauss(x[k], xcen, wconv, 1.0);
            n     = n + g;
            yp[k] = y[k]*g;
        }
        //temp[j] = trapz_int(x,yp, num);
        temp[j] = 0.0;
        for (k=0;k<num;k++){
            temp[j] = temp[j] + yp[k];
        }
        temp[j] = temp[j] / n;
    }
    for (j=0;j<num;j++){
        y[j] = temp[j];
    }
    free(yp);
    free(temp);
    return(SUCCESS);
}

/******************************************************************************
* convolve_pad()
* Simple gaussian convolution (padded)
*
* - x array of length num
* - y array of length num
* - wconv is the convolution width (gaussian fwhm)
*
* Returns
* -------
* - SUCCESS/FAILURE
*
* Returns
* -------
* - SUCCESS/FAILURE
*
* Notes
* -----
* Note y is the input 'data' array and is recomputed to 
* hold the convolution.  ie this overwrites the y data
*
* This routine pads the arrays to get rid of end point probs
* Note this is not normalized like the above one,
* Therefore it adds wierd stuff!  (NEEDS WORK) 
*
* Have to add normalization stuff to loop as did above
* and should work ok.  Should really be doing this with fft!
*
******************************************************************************/
int convolve_pad(double *x, double *y, int num, double wconv){
    double  *yp,*xp,*ypp, *temp, del, xcen;
    int j, k, np;

    np = 3*num - 2;
    xp = new_array(np, double);
    yp = new_array(np, double);
    ypp = new_array(np, double);
    temp = new_array(num, double);
    del = fabs(x[num-1] - x[0])/(num-1);

    for (j=0;j<np;j++){
        if (j < num - 1){
            xp[j] = x[0] - del*(num-1) + del*j;
            yp[j] = y[0] ;
        } else if (j > 2*(num - 1) ) {
            xp[j] = x[num-1] + del*( j - 2*(num - 1) );
            yp[j] = y[num-1] ;
        } else {
            xp[j] = x[j-num + 1];
            yp[j] = y[j-num + 1];
        }
    }

    for (j=0;j<num;j++){
        xcen = x[j];
        for (k=0;k<np;k++){
            ypp[k] = yp[k]*gauss(xp[k], xcen, wconv, 1.0);
        }
        temp[j] = trapz_int(xp,ypp,np);
    }
    for (j=0;j<num;j++){
        y[j] = temp[j];
    }
    free(xp);
    free(yp);
    free(ypp);
    free(temp);
    return(SUCCESS);
}

/******************************************************************************
* normalize()
* Normalize an array to a specified value at a point given by the ordinate
*
* Parameters
* ---------
* - x array of length num
* - y array of length num
* - xnorm is the closest x position for normalization
* - norm_val is the value to normalize the data to at xnorm
*
* Returns
* -------
* - SUCCESS/FAILURE
* 
* Notes
* -----
* This modifies the values in the y-array
*
******************************************************************************/
int norm_array(double *x, double *y, int num, double xnorm, double norm_val){
    double  norm;
    int j;

    j = 0;
    while (j < num-1){
        if (x[j] < xnorm){
            j++;
        } else {
            break;
        }
    }
    if (y[j] == 0.0) return(FAILURE);
    norm = norm_val/y[j];
    for (j=0;j<num;j++){
        y[j] = y[j] * norm;
    }
    return (SUCCESS);
}

/******************************************************************************
* print_martrix()
* Display a matrix
******************************************************************************/
//void print_matrix(int m, int n, double *y) 
//{
//  int i, j, k;
//  k = 0;
//  for( i=0; i<m; i++) {
//    for( j=0; j<n; j++) { 
//      buff_print( "%14.8g", y[k++]);
//    }
//    buff_print( "\n" );
//  }
//  return;
//}
/******************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//  Complex fncs
///////////////////////////////////////////////////////////////////////////////

//typedef struct FCOMPLEX {float r,i;} fcomplex;

fcomplex Cadd(fcomplex a, fcomplex b)
{
    fcomplex c;
    c.r=a.r+b.r;
    c.i=a.i+b.i;
    return c;
}

fcomplex Csub(fcomplex a, fcomplex b)
{
    fcomplex c;
    c.r=a.r-b.r;
    c.i=a.i-b.i;
    return c;
}


fcomplex Cmul(fcomplex a, fcomplex b)
{
    fcomplex c;
    c.r=a.r*b.r-a.i*b.i;
    c.i=a.i*b.r+a.r*b.i;
    return c;
}

fcomplex Complex(float re, float im)
{
    fcomplex c;
    c.r=re;
    c.i=im;
    return c;
}

fcomplex Conjg(fcomplex z)
{
    fcomplex c;
    c.r=z.r;
    c.i = -z.i;
    return c;
}

fcomplex Cdiv(fcomplex a, fcomplex b)
{
    fcomplex c;
    float r,den;
    if (fabs(b.r) >= fabs(b.i)) {
        r=b.i/b.r;
        den=b.r+r*b.i;
        c.r=(a.r+r*a.i)/den;
        c.i=(a.i-r*a.r)/den;
    } else {
        r=b.r/b.i;
        den=b.i+r*b.r;
        c.r=(a.r*r+a.i)/den;
        c.i=(a.i*r-a.r)/den;
    }
    return c;
}

float Cabs(fcomplex z)
{
    float x,y,ans,temp;
    x=fabs(z.r);
    y=fabs(z.i);
    if (x == 0.0)
        ans=y;
    else if (y == 0.0)
        ans=x;
    else if (x > y) {
        temp=y/x;
        ans=x*sqrt(1.0+temp*temp);
    } else {
        temp=x/y;
        ans=y*sqrt(1.0+temp*temp);
    }
    return ans;
}

fcomplex Csqrt(fcomplex z)
{
    fcomplex c;
    float x,y,w,r;
    if ((z.r == 0.0) && (z.i == 0.0)) {
        c.r=0.0;
        c.i=0.0;
        return c;
    } else {
        x=fabs(z.r);
        y=fabs(z.i);
        if (x >= y) {
            r=y/x;
            w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
        } else {
            r=x/y;
            w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
        }
        if (z.r >= 0.0) {
            c.r=w;
            c.i=z.i/(2.0*w);
        } else {
            c.i=(z.i >= 0) ? w : -w;
            c.r=z.i/(2.0*c.i);
        }
        return c;
    }
}

fcomplex RCmul(float x, fcomplex a)
{
    fcomplex c;
    c.r=x*a.r;
    c.i=x*a.i;
    return c;
}

/******************************************************************************/

