"""
Fcns for peak fitting

Authors/Modifications:
----------------------
* T. Trainor (tptrainor@alaska.edu)

Todo:
-----
* add background
* more peak shape profiles
* numerical integration and error analysis
 
"""
#######################################################################

import types
import numpy as num
import scipy

#######################################################################
def gauss(x, xcen, fwhm, mag):
    """
    Calculate a gaussian profile:

    Notes:
    ------
    Computes the following
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
    calculate a lorentzian profile

    Notes:
    ------
    y = mag * sigma**2 / (sigma**2 + (x-xcen)**2))
    y = mag / (1 + ( (x-xcen)/sigma )**2 )
    """
    if (fwhm == 0.0): return(0.0)
    a =   ( x - xcen )  / ( 0.5 * fwhm ) 
    y = mag / ( 1 + a**2. ) 
    return(y)

#######################################################################
def voigt( x, xcen, fwhm, mag, flor=0.0):
    """
    Calculate a psuedo-voigt profile:

    Notes:
    ------
    y = flor*lor + (1-flor)*gauss
    
    This approximates the voigt-profile (which is a convolution of
    a gaussian and lorentzian)

    """
    g = gauss( x, xcen, fwhm, mag )
    l = lor( x, xcen, fwhm, mag )
    p = flor * l  +  ( 1 - flor ) * g
    return (p)    

#######################################################################
class Peak:
    """
    Peak class
    """
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

    Notes:
    ------
    Regression model is based on equation of line
      y = mx + b
    x and y should be simple double arrays

    Examples:
    ---------
    >>lr = LinReg(x,y,plot=True)

    Note also see the linear regression function
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
        from matplotlib import pyplot
        print "slope     = ", self.m, " +/- ", self.s_m
        print "intercept = ", self.b, " +/- ", self.s_b
        if x == None: x = self.x
        pyplot.clf()
        ycalc = self.calc_y(x)
        pyplot.plot(x,ycalc)
        if y != None:
            pyplot.plot(x,y,'bo')
            resid = y - ycalc
            pyplot.plot(x,resid, 'ro')
            pyplot.plot([x[0],x[n-1]],[0,0],'k-')
    
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
        nc    = float(self.n)
        yc_ave = self.y_ave

        n     = float(len(y))
        y_ave = num.sum(y)/n

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
    from matplotlib import pyplot
    pyplot.clf()
    pyplot.plot(x,y,'m.-')
    
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
    
    from matplotlib import pyplot
    pyplot.clf()
    pyplot.plot(x,y,'o')
    pyplot.plot(x,yc,'-')

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
    
