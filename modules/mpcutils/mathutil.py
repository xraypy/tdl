#######################################################################
"""
T. Trainor (fftpt@uaf.edu)
Simple math addons and wrappers

Modifications:
--------------
- minimize and random from original tdl, written by Matt Newville

"""
#######################################################################
"""
Todo

 - peak fit
 
"""
#######################################################################

import types
import numpy as num
import scipy

#######################################################################
# some simple stuff
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

#######################################################################
def minimize(f,x,y,params,*args,**kws):
    """
    Simple wrapper around scipy.optimize.leastsq
      f is the function to be optimized
      x is a vector of independant varibles (floats) - the abscissa
      y is the corresponding vector of known/dependant values - the ordinate
      params is a tuple of doubles which are to be optimized.
      args and kws are additional arguments for f

    use:
         params = minimize(f,x,y,params,*args,**kw)
    where
         ycalc = f(x,*args,**kw)
    and
         args should be all single valued (floats)

    E.g. Define a function and optimize (a,b)
    def fun(x,a,b,c,d=1):
       ...calc y... 
       return y
    (a,b) = minimize(f,x,yobs,(a,b),c,d=10)
    
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
    wrap numpy random seed
    Seds the random number generator
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
    returns npts random numbers.
    a,b,c are default arguments for the dist functions
    e.g. NR.normal a = mean, b = stdev of the distrobution
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

#######################################################################
# Peaks
#######################################################################

#######################################################################
def gauss(x, xcen, fwhm, mag):
    """
    calculate a gaussian profile:
       y = mag * exp( (x-xcen)**2 / (2 * sig**2))
    note the fwhm is related to sigma by:
       sigma = fwhm/2.35482
         or
       2*sig**2 = (sqrt(2)*fwhm/2.35482) = 0.600561*fwhm
    """
    # note this makes it safe 
    # and makes it perform like a 
    # Kroniker-delta function
    if fwhm == 0.0:
        xx  = num.fabs(x-xcen)
        idx = num.where(xx == xx.min())
        idx = idx[0]
        y = num.zeros(len(x))
        y[idx] = mag
        return y
    a =   ( x - xcen )  / ( 0.600561 * fwhm ) 
    y =   mag *  num.exp( -1. * (a**2.) ) 

    return(y)

#######################################################################
def lor( x, xcen, fwhm, mag):
    """
    calculate a lorentzian profile:
       y = mag * sigma**2 / (sigma**2 + (x-xcen)**2))
       y = mag / (1 + ( (x-xcen)/sigma )**2 )
    """
    if (fwhm == 0.0): return(0.0)
    a =   ( x - xcen )  / ( 0.5 * fwhm ) 
    y = mag / ( 1 + a**2. ) 
    return(y)

#######################################################################
def voigt( x, xcen, fwhm, mag, *args):
#def voigt( x, xcen, fwhm, mag, flor):
    """
    Calculate a psuedo-voigt profile:
       y = flor*lor + (1-flor)*gauss
    This approximates the voigt-profile (which is a convolution of
    a gaussian and lorentzian)

    If flor is not passed as an argument it is set to zero
    """
    if len(args) > 0:
        flor = args[0]
    else:
        flor = 0.0
    g = gauss( x, xcen, fwhm, mag )
    l = lor( x, xcen, fwhm, mag )
    p = flor * l  +  ( 1 - flor ) * g
    return (p)    

#######################################################################
class Peak:
    def __init__(self,npeaks=0):
        """
        self.bgr_params[0],   lin bgr offset
        self.bgr_params[1],   lin bgr slope
        self.pk_params[0][0], center
        self.pk_params[0][1], fwhm
        self.pk_params[0][2], mag
        self.pk_params[0][3], frac lor
        etc..
        """
        self.npeaks      = npeaks
        self.bgr_include = 1
        self.bgr_params  = [0.0, 0.0]
        self.pk_params   = []
        self.pk_include  = []
        for j in range(npeaks):
            self.pk_include.append(1)
            self.pk_params.append([0.0, 0.0, 0.0, 0.0]) 

    def set_bgr(self,offset=None,slope=None):
        if self.bgr_include == 0:
            self.bgr_include = 1
        if offset != None: self.bgr_params[0] = offset
        if slope  != None: self.bgr_params[1] = slope

    def set_peak(self,idx=-1,cen=None,fwhm=None,mag=None,flor=None,include=None):
        if idx < 0:
            self.pk_include.append(1)
            self.pk_params.append([0.0,0.0,0.0,0.0])
            idx = len(self.pk_params) - 1

        self.npeaks = len(self.pk_params)
        if idx > self.npeaks - 1:
            idx = self.npeaks - 1
        
        if cen     != None: self.pk_params[idx][0] = cen
        if fwhm    != None: self.pk_params[idx][1] = fwhm
        if mag     != None: self.pk_params[idx][2] = mag
        if flor    != None: self.pk_params[idx][3] = flor
        if include != None: self.pk_include[idx] = include

    def calc(self,x):
        x = num.array(x)
        y = num.zeros(len(x))
        # bgr        
        if self.bgr_include == 1:
            offset = self.bgr_params[0]
            slope  = self.bgr_params[1]
            y = line(x, offset, slope)
        # pks
        self.npeaks = len(self.pk_params)
        for j in range(self.npeaks):
            if self.pk_include[j] == 1:
                xcen = self.pk_params[j][0]
                fwhm = self.pk_params[j][1]
                mag  = self.pk_params[j][2]
                flor = self.pk_params[j][3]
                print j, xcen, fwhm, mag, flor
                y = y + voigt(x, xcen, fwhm, mag, flor)
        return y

    def fit(self,x,yobs):
        pass

#######################################################################
class LinReg:
    """
    Equations for calculating linear regression
      y = mx + b
    x and y should be simple double arrays
    call as:
    >>lr = LinReg(x,y,plot=True)

    ###
    Note also see the linear regression funciton
    in scipy (=> stats.linregress), for example:
    >>from scipy.stats import linregress
    >>m = 1.4
    >>b = 10.0
    >>x = num.linespace(1,10)
    >>y = m*x + b
    >>(m_s,b_s,r,tt,stderr)=linregress(x,y)
    >>print('Linear regression using stats.linregress')
    >>print('parameters: m=%.2f   b=%.2f \n
            regression: m_s=%.2f b_s=%.2f, std error= %.3f' %
            (m,b,m_s,b_s,stderr))
    """
    ########################################################
    def __init__(self,x,y,plot=False):
        self.fit(x,y,plot=plot)
        
    ########################################################
    def fit(self, x, y, plot=False):
        """
        fit lr
        """
        # store x 
        self.x = x
        
        # calc stats
        n     = len(x)
        x_ave = sum(x) / n
        y_ave = sum(y) / n
        Sxx   = sum( (x-x_ave)**2)
        Syy   = sum( (y-y_ave)**2)
        Sxy   = sum( (x-x_ave)*(y-y_ave))

        # slope and intercept
        m = Sxy/Sxx
        b = y_ave - m*x_ave

        # Residual and residual standard deviation
        resid = (y - (b + m*x) )
        SS_resid = sum(resid**2)
        s_r = num.sqrt(SS_resid/(n-2))

        # std dev of slope and intercept
        s_m = num.sqrt(s_r**2/Sxx)
        s_b = s_r * num.sqrt( 1 /(n - (sum(x))**2 / sum(x**2) ) )

        self.m=m
        self.b=b
        self.s_r=s_r
        self.Sxx=Sxx
        self.Syy=Syy
        self.Sxy=Sxy
        self.y_ave=y_ave
        self.x_ave=x_ave
        self.n=n
        self.s_m=s_m
        self.s_b=s_b
        
        if plot:
            self.plot(x,y)

    ########################################################
    def plot(self,x=None,y=None):
        import pylab
        print "slope     = ", self.m, " +/- ", self.s_m
        print "intercept = ", self.b, " +/- ", self.s_b
        if x == None: x = self.x
        pylab.clf()
        ycalc = self.calc_y(x)
        pylab.plot(x,ycalc)
        if y != None:
            pylab.plot(x,y,'bo')
            resid = y - ycalc
            pylab.plot(x,resid, 'ro')
            pylab.plot([x[0],x[n-1]],[0,0],'k-')
    
    ########################################################
    def calc_y(self,x=None):
        """
        calulate y given x
        y = mx + b
        """
        if x == None:
            x = self.x
        m = self.m
        b = self.b
        return (m*x + b)

    ########################################################
    def calc_x(self,y):
        """
         given cal data calc x given y
         x = (y-b)/m
        """
        m = self.m
        b = self.b
        return ( (y-b)/m )

    ########################################################
    def calc_x_err(self,y):
        """
        calc x given y
         x = (y-b)/m
         this function also computes the error
         and allows for y to be from several measurements
        """
        m     = self.m
        b     = self.b
        s_r   = self.s_r
        Sxx   = self.Sxx
        nc    = self.n
        yc_ave = self.y_ave

        n     = len(y)
        y_ave = num.sum(y)/float(n)

        x = self.calc_x(y_ave)

        s_x = (s_r/m) * num.sqrt( (1/n) + (1/nc) + (y_ave - yc_ave)**2 / (Sxx * m**2 ) )

        return ([x,s_x])

#######################################################################
#######################################################################
def test_lr():
    """
    test linear regression
    """
    n      = 100
    slope  = 7.0
    offset = 0.23
    x   = num.linspace(1,10,num=n)
    yo  = offset + slope*(x)
    y   = offset + slope*(x+num.random.randn(n)/5.)
    #lr = LinReg(x,y,plot=True)
    line = LinReg(x,y,plot=True).calc_y()
    print line
    
def test_peak():
    """
    test peak
    """
    n      = 100
    slope  = 7.0
    offset = 0.23
    x = num.linspace(1,10,num=n)
    p = Peak()
    p.set_bgr(slope=slope,offset=offset)
    p.set_peak(cen=3.,fwhm=1.2,mag=20.,flor=.3)
    p.set_peak(cen=6.,fwhm=0,mag=35.,flor=.3)
    y = p.calc(x)
    import pylab
    pylab.clf()
    pylab.plot(x,y,'m.-')
    
def test_minimize():
    """
    test minimize
    """
    n      = 100
    slope  = 7.0
    offset = 0.23
    x   = num.linspace(1,10,num=n)
    y   = offset + slope*(x+num.random.randn(n)/5.)

    # function    
    def fun(x,a,b,c,d=1):
        yc = x*a + b
        return yc

    # parameter set
    a = 1.
    b = 2.
    c = 3.
    d = 4.
    # these are initial guess of the two to be optimized
    p = (a,b)
    p = minimize(fun,x,y,(a,b),c,d=d)
    print p
    a,b = p
    yc = fun(x,a,b,c,d=d)
    
    import pylab
    pylab.clf()
    pylab.plot(x,y,'o')
    pylab.plot(x,yc,'-')

#######################################################################
#######################################################################
if __name__ == "__main__":
    print "Test lr"
    test_lr()
    raw_input('hit enter')
    print "Test peak"
    test_peak()
    raw_input('hit enter')
    print "Test minimize"
    test_minimize()
    
