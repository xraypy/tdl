"""
Simple math addons and wrappers

Authors/Modifications:
----------------------
* T. Trainor (tptrainor@alaska.edu)
  M. Newville (newville@cars.uchicago.edu)
* minimize and random from original tdl

Todo:
-----
* peak fit
 
"""
#######################################################################

import types
import numpy as num
import scipy

#######################################################################
def ave(x):
    """
    average of an array
    """
    #return (sum(x)/float(len(x)))
    return num.ave(x)

def std(x):
    """
    standard deviation of an array
    """
    #x_ave = self.ave(x)
    #return( num.sqrt( sum( (x-x_ave)**2 ) / float(len(x)) ) )
    return num.std(x)

def line(x, offset, slope ):
    """
    calculation of a line
    """
    y =   slope * x + offset
    return(y)

def square(a):
    """
    square of two numbers
    """
    return a*a

def cosd(x):
    """
    num.cos(x), x in degrees
    """
    return num.cos(num.radians(x))

def sind(x):
    """
    num.sin(x), x in degrees
    """
    return num.sin(num.radians(x))

def tand(x):
    """
    num.tan(x), x in degrees
    """
    return num.tan(num.radians(x))

def arccosd(x):
    """
    num.arccos(x), result returned in degrees
    """
    return num.degrees(num.arccos(x))

def arcsind(x):
    """
    num.arcsin(x), result returned in degrees
    """
    return num.degrees(num.arcsin(x))

def arctand(x):
    """
    num.arctan(x), result returned in degrees
    """
    return num.degrees(num.arctan(x))

def cartesian_mag(v):
    """
    Calculate the norm of a vector defined in
    a cartesian basis.
    
    This should give same as num.linalg.norm
    """
    m = num.sqrt(num.dot(v,v))
    return m

def cartesian_angle(u,v):
    """
    Calculate angle between two vectors defined in
    a cartesian basis.

    Result is always between 0 and 180 degrees
    """
    uv = num.dot(u,v)
    um = cartesian_mag(u)
    vm = cartesian_mag(v)
    denom = (um*vm)
    if denom == 0: return 0.
    arg = uv/denom
    if num.fabs(arg) > 1.0:
        arg = arg / num.fabs(arg)
    alpha = arccosd(arg)
    return alpha

#######################################################################
def minimize(f,x,y,params,*args,**kws):
    """
    Simple wrapper around scipy.optimize.leastsq

    Parameters:
    -----------
    * f is the function to be optimized
    * x is a vector of independant varibles (floats) - the abscissa
    * y is the corresponding vector of known/dependant values - the ordinate
    * params is a tuple of doubles which are to be optimized.
    * args and kws are additional arguments for f

    Notes:
    ------
    >>params = minimize(f,x,y,params,*args,**kw)

    where
        ycalc = f(x,*args,**kw)
    and
        args should be all single valued (floats)

    Examples:
    ---------
    # Define a function and optimize (a,b)
    >>def fun(x,a,b,c,d=1):
    >>   ...calc y... 
    >>   return y
    >>(a,b) = minimize(f,x,yobs,(a,b),c,d=10)
    
    """
    from scipy.optimize import leastsq
    XX    = x
    YY    = y
    FUNC  = f
    ###########################################    
    def _residual(parameters,*arguments):
        """
        if the last arg is a dictionary assume
        its the kw args for the function
        """
        kw = {}
        if len(arguments) > 0:
            if type(arguments[-1]) == types.DictionaryType:
                kw  = arguments[-1]
                arguments = arguments[0:-1]
        # Now combine all parameters into a single tuple
        parameters = tuple(parameters) + tuple(arguments)
        # calculate theory
        yc = FUNC(XX,*parameters,**kw)
        #return residual
        return (YY-yc)
    ###########################################
    # make sure params is a tuple
    params = tuple(params)
    args   = args + (kws,)
    test   = _residual(params,*args)
    if len(test) != len(x):
        print 'cannot minimize function '
    
    result = leastsq(_residual,params,args=args)
    return  result[0]

#######################################################################
def random_seed(x=None):
    """
    wrapper for numpy random seed
    Seeds the random number generator
    """
    if x is None:
        return num.random.seed()
    else:
        try:
            return num.random.seed([x])
        except:
            return num.random.seed()

def random(a=1,b=1,c=1,npts=1,distribution='normal',**kw):
    """
    wrapper for numpy random distributions

    Parameters:
    -----------
    * a,b,c are default arguments for the dist functions
      e.g. NR.normal a = mean, b = stdev of the distrobution
    * npts is the number of points
    
    Outputs:
    --------
    returns npts random numbers.
    """ 
    NR = num.random
    if   distribution == 'binomial':        return NR.binomial(a,b,size=npts)
    elif distribution == 'geometric':       return NR.geometric(a,size=npts)    
    elif distribution == 'poisson':         return NR.poisson(a,size=npts)    
    elif distribution == 'zipf':            return NR.zipf(a,size=npts)    
    elif distribution == 'beta':            return NR.beta(a,b,size=npts)    
    elif distribution == 'chisquare':       return NR.chisquare(a,size=npts)    
    elif distribution == 'exponential':     return NR.exponential(a,size=npts)
    elif distribution == 'gamma':           return NR.gamma(a,b,size=npts)
    elif distribution == 'gumbel':          return NR.gumbel(a,b,size=npts)
    elif distribution == 'laplace':         return NR.laplace(a,b,size=npts)
    elif distribution == 'lognormal':       return NR.lognormal(a,b,size=npts)    
    elif distribution == 'logistic':        return NR.logistic(a,b,size=npts)    
    elif distribution == 'multivariate_normal':  return NR.multivariate_normal(a,b,size=npts)
    elif distribution == 'noncentral_chisquare': return NR.noncentral_chisquare(a,b,size=npts)
    elif distribution == 'noncentral_f':    return NR.noncentral_f(a,b,c,size=npts)
    elif distribution == 'normal':          return NR.normal(a,b,size=npts)
    elif distribution == 'pareto':          return NR.pareto(a,size=npts)
    elif distribution == 'power':           return NR.power(a,size=npts)
    elif distribution == 'randint':         return NR.randint(a,b,size=npts)
    elif distribution == 'random_integers': return NR.random_integers(a,b,size=npts)
    elif distribution == 'rayleigh':        return NR.rayleigh(a,size=npts)
    elif distribution == 'standard_cauchy': return NR.standard_cauchy(size=npts)
    elif distribution == 'standard_exponential': return NR.standard_exponential(size=npts)
    elif distribution == 'standard_gamma':  return NR.standard_gamma(a,size=npts)
    elif distribution == 'standard_normal': return NR.standard_normal(size=npts)
    elif distribution == 'standard_t':      return NR.standard_t(a,size=npts)
    elif distribution == 'uniform':         return NR.uniform(a,b,size=npts)
    elif distribution == 'wald':            return NR.wald(a,b,size=npts)
    elif distribution == 'weibull':         return NR.weibull(a,b,size=npts)

#######################################################################
