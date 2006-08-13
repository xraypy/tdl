"""
This module contains miscellaneous math functions

Author:
   Mark Rivers

Created: 
   Sept. 16, 2002

Modifications:
   Sept. 26, 2002  MLR
      - Added newton function from scipy.optimize.  Put here so users don't
        need scipy
"""

import Numeric
import LinearAlgebra

############################################################
def polyfitw(x, y, w, ndegree, return_fit=0):
   """
   Performs a weighted least-squares polynomial fit with optional error estimates.

   Inputs:
      x: 
         The independent variable vector.

      y: 
         The dependent variable vector.  This vector should be the same 
         length as X.

      w: 
         The vector of weights.  This vector should be same length as 
         X and Y.

      ndegree: 
         The degree of polynomial to fit.

   Outputs:
      If return_fit==0 (the default) then polyfitw returns only C, a vector of 
      coefficients of length ndegree+1.
      If return_fit!=0 then polyfitw returns a tuple (c, yfit, yband, sigma, a)
         yfit:  
            The vector of calculated Y's.  Has an error of + or - Yband.

         yband: 
            Error estimate for each point = 1 sigma.

         sigma: 
            The standard deviation in Y units.

         a: 
            Correlation matrix of the coefficients.

   Written by:   George Lawrence, LASP, University of Colorado,
                 December, 1981 in IDL.
                 Weights added, April, 1987,  G. Lawrence
                 Fixed bug with checking number of params, November, 1998, 
                 Mark Rivers.  
                 Python version, May 2002, Mark Rivers
   """
   n = min(len(x), len(y)) # size = smaller of x,y
   m = ndegree + 1         # number of elements in coeff vector
   a = Numeric.zeros((m,m),Numeric.Float)  # least square matrix, weighted matrix
   b = Numeric.zeros(m,Numeric.Float)    # will contain sum w*y*x^j
   z = Numeric.ones(n,Numeric.Float)     # basis vector for constant term

   a[0,0] = Numeric.sum(w)
   b[0] = Numeric.sum(w*y)

   for p in range(1, 2*ndegree+1):     # power loop
      z = z*x   # z is now x^p
      if (p < m):  b[p] = Numeric.sum(w*y*z)   # b is sum w*y*x^j
      sum = Numeric.sum(w*z)
      for j in range(max(0,(p-ndegree)), min(ndegree,p)+1):
         a[j,p-j] = sum

   a = LinearAlgebra.inverse(a)
   c = Numeric.matrixmultiply(b, a)
   if (return_fit == 0):
      return c     # exit if only fit coefficients are wanted

   # compute optional output parameters.
   yfit = Numeric.zeros(n,Numeric.Float)+c[0]   # one-sigma error estimates, init
   for k in range(1, ndegree +1):
      yfit = yfit + c[k]*(x**k)  # sum basis vectors
   var = Numeric.sum((yfit-y)**2 )/(n-m)  # variance estimate, unbiased
   sigma = Numeric.sqrt(var)
   yband = Numeric.zeros(n,Numeric.Float) + a[0,0]
   z = Numeric.ones(n,Numeric.Float)
   for p in range(1,2*ndegree+1):     # compute correlated error estimates on y
      z = z*x		# z is now x^p
      sum = 0.
      for j in range(max(0, (p - ndegree)), min(ndegree, p)+1):
         sum = sum + a[j,p-j]
      yband = yband + sum * z      # add in all the error sources
   yband = yband*var
   yband = Numeric.sqrt(yband)
   return c, yfit, yband, sigma, a


############################################################
def fit_gaussian(chans, counts):
   """
   Fits a peak to a Gaussian using a linearizing method
   Returns (amplitude, centroid, fwhm).
   Inputs:
      chans:
         An array of x-axis coordinates, typically channel numbers

      counts:
         An array of y-axis coordinates, typically counts or intensity

   Outputs:
      Returns a tuple(amplitude, centroid, fwhm)
      amplitude:
         The peak height of the Gaussian in y-axis units
      centroid:
         The centroid of the gaussian in x-axis units

      fwhm:
         The Full Width Half Maximum (FWHM) of the fitted peak

   Method:
      Takes advantage of the fact that the logarithm of a Gaussian peak is a
      parabola.  Fits the coefficients of the parabola using linear least
      squares.  This means that the input Y values (counts)  must not be 
      negative.  Any values less than 1 are replaced by 1.
   """
   center = (chans[0] + chans[-1])/2.
   x = Numeric.asarray(chans, Numeric.Float)-center
   y = Numeric.log(Numeric.clip(counts, 1, max(counts)))
   w = Numeric.asarray(counts, Numeric.Float)**2
   w = Numeric.clip(w, 1., max(w))
   fic = polyfitw(x, y, w, 2)
   fic[2] = min(fic[2], -.001)  # Protect against divide by 0
   amplitude = Numeric.exp(fic[0] - fic[1]**2/(4.*fic[2]))
   centroid  = center - fic[1]/(2.*fic[2])
   sigma     = Numeric.sqrt(-1/(2.*fic[2]))
   fwhm      = 2.35482 * sigma
   return amplitude, centroid, fwhm


############################################################
def compress_array(array, compress):
   """
   Compresses an 1-D array by the integer factor "compress".  
   Temporary fix until the equivalent of IDL's 'rebin' is found.
   """

   l = len(array)
   if ((l % compress) != 0):
      print 'Compression must be integer divisor of array length'
      return array

   temp = Numeric.resize(array, (l/compress, compress))
   return Numeric.sum(temp, 1)/compress

############################################################
def expand_array(array, expand, sample=0):
   """
   Expands an 1-D array by the integer factor "expand".  
   if 'sample' is 1 the new array is created with sampling, if 1 then
   the new array is created via interpolation (default)
   Temporary fix until the equivalent of IDL's 'rebin' is found.
   """

   l = len(array)
   if (expand == 1): return array
   if (sample == 1): return Numeric.repeat(array, expand)

   kernel = Numeric.ones(expand, Numeric.Float)/expand
   # The following mimic the behavior of IDL's rebin when expanding
   temp = Numeric.convolve(Numeric.repeat(array, expand), kernel, mode=2)
   # Discard the first "expand-1" entries
   temp = temp[expand-1:]
   # Replace the last "expand" entries with the last entry of original
   for i in range(1,expand): temp[-i]=array[-1]
   return temp

##############################################################
def newton(func, x0, fprime=None, args=(), tol=1.48e-8, maxiter=50):
    """
    newton is taken directly from scipy.optimize.minpack.py.  
    I have extracted it here so that users of my software don't have to install 
    scipy.  This may change in the future as I use more scipy features, and 
    scipy will be required

    Given a function of a single variable and a starting point,
    find a nearby zero using Newton-Raphson.
    fprime is the derivative of the function.  If not given, the
    Secant method is used.
    """
    if fprime is not None:
        p0 = x0
        for iter in range(maxiter):
            myargs = (p0,)+args
            fval = func(*myargs)
            fpval = fprime(*myargs)
            if fpval == 0:
                print "Warning: zero-derivative encountered."
                return p0
            p = p0 - func(*myargs)/fprime(*myargs)
            if abs(p-p0) < tol:
                return p
            p0 = p
    else: # Secant method
        p0 = x0
        p1 = x0*(1+1e-4)
        q0 = apply(func,(p0,)+args)
        q1 = apply(func,(p1,)+args)
        for iter in range(maxiter):
            try:                
                p = p1 - q1*(p1-p0)/(q1-q0)
            except ZeroDivisionError:
                if p1 != p0:
                    print "Tolerance of %g reached" % (p1-p0)
                return (p1+p0)/2.0
            if abs(p-p0) < tol:
                return p
            p0 = p1
            q0 = q1
            p1 = p
            q1 = apply(func,(p1,)+args)
    raise RuntimeError, "Failed to converge after %d iterations, value is %f" % (maxiter,p)

