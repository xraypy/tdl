#######################################################################
"""
Tom Trainor (tptrainor@alaska.edu)

Methods to handle generic backgorund determination

Modifications:
--------------


"""
#######################################################################
"""
Todo
   
"""
#######################################################################

import types
import copy
import numpy as num
from matplotlib import pyplot

from  mpcutils.mathutil import LinReg

#######################################################################
def linear_background(data,nbgr=0):
    """
    Calculate a linear background of the data based on the endpoints

    Note we assume that the data is on an evenly spaced grid so no
    abscica is used in calculating the background line
    """
    ndat = len(data)
    if nbgr <= 0:
        return num.zeros(ndat)
    if ndat < 2*nbgr + 1:
        return num.zeros(ndat)
    # calc linear bgr from end points
    xlin = num.arange(0,nbgr,1,dtype=float)
    xlin = num.append(xlin,num.arange(ndat-nbgr,ndat,1))
    ylin = num.array(data[0:nbgr],dtype=float)
    ylin = num.append(ylin,data[ndat-nbgr:])
    lr   = LinReg(xlin,ylin,plot=False)
    linbgr = lr.calc_y(num.arange(ndat))
    return linbgr

#######################################################################
def background(data,nbgr=0,width=0,pow=0.5,tangent=False,debug=False):
    """
    Calculate the background for a given line, y.

    This is a simplified form of the algorithm used for fitting
    XRF backgrounds (see Kajfosz and Kwiatek (1987) Nuc.
    Instr. Meth. in Phys. Res., B22, 78-81)

    At each data point (i) a calculated polynomial is brought up from underneath
    untill its first contact with any data point within the polynomial range
    (+/-width).  The distance between the apex of the polynomial (which occurs
    at point i) and zero is taken as the background for point i.

    This algorithm also allows for the inclusion of a linear background
    based on the end points.

    Note we assume that data is on a uniform
    grid, ie delta_x steps between data values are all the same,
    therefore no abscissa is passed in here.
    (therefore you might need to spline your data to a
    uniform grid for this algo to work!)

    * data is the data (ie y-data or ordinate).  

    * nbgr is the number of end points to use for calculation of a linear
      background.  This part of the background is removed from the data
      before polynomials are adjusted.  We then add it back to the polynomial
      sum so the total background returned from this function includes both
      
    * width is the polynomial (half) width in units of steps
      in y.  If pow = 0.5, width is the radius of a circle.
      Note that the extent of the polynomial used in bgr fitting
      for a given point y[i] is given by 2*width+1.  If width=0
      then we just calculate and return a linear bgr based on
      the end points

    * pow is the power of the polynomial (pow>=0).
      default value of 0.5 gives a circle
      flatter polynomials will result with  pow < 0.5 (pow = 0 are linear)
      steeper polynomials will result with  pow  > 0.5

    * tangent is a flag (True/False) to indicate if we add the
      average local linear slope of the data to the polynomial
      (this effectively removes the local slope of the data while
      adjusting the position of the polynomials)
    
    * debug is a flag (True/False) to indicate if additional debug arrays
      should be calculated
    
    Note given the polynomial:
        poly    = (r**2. - delx**2)**pow  - (r**2.)**pow
    calc
        d^2(poly)/dx^2 = 0
    to determine where this function will have the max slope.
    The result is:
       delx_max_slope = r/sqrt(2*pow-1)
    If we assume that the 'width' argument corresponds to this
    delx_max_slope, then we can calc r for the polynomial as:
       r = width*sqrt(2*pow-1)
    (if pow<=1/2 just use r = width)

    """
    # make sure pow is positive
    # and create some debug stuff
    if pow < 0.:
        print "Warning power is less than 0, changing it to positive"
        pow = -1.*pow
    if debug:
        p = []
        d = []
    
    # linear bgr subtract data
    linbgr = linear_background(data,nbgr=nbgr)
    if width <= 0.: return linbgr
    y = data - linbgr

    # create bgr array
    ndat = len(y)
    bgr  = num.zeros(ndat)
    
    # calc polynomial
    npoly   = int(2*width) + 1
    pdelx   = num.array(range(npoly),dtype=float) - float((npoly-1.)/2.)
    if pow <= 0.5:
        r = float(width)
    else:
        r = float(width)*sqrt(2.*pow-1.)
    poly = (r**2. - pdelx**2.)**pow  - (r**2.)**pow
    # renorm poly
    pmax = num.max(num.fabs(poly))
    poly = ((width*pow)/pmax) * poly
    
    # loop through each point
    #delta = num.zeros(len(poly))
    n = (npoly-1)/2
    for j in range(ndat):
        # data and polynomial indicies
        dlidx = num.max((0,j-n))
        dridx = num.min((ndat,j+n+1))
        plidx = num.max((0,n-j))
        pridx = num.min((npoly,ndat-j+n))
        if tangent:
            # calc avg val to l and r of center
            # and use to calc avg slope
            nl    = len(y[dlidx:j])
            if nl == 0:
                lyave = 0.0
                lxave = 0.0
            else:
                lyave = num.sum(y[dlidx:j])/nl
                lxave = num.sum(num.arange(dlidx,j))
            nr    = len(y[j+1:dridx])
            if nr == 0:
                ryave = 0.0
                rxave = 0.0
            else:
                ryave = num.sum(y[j+1:dridx])/nr
                rxave = num.sum(num.arange(j+1,dridx))
            slope = (ryave - lyave)/ num.abs(rxave - lxave)
            line = slope*num.arange(-1*nl,nr+1) 
            #print line
        else:
            line = num.zeros(len(y[dlidx:dridx]))
        delta  = y[dlidx:dridx] - (y[j] + poly[plidx:pridx]+line) 
        bgr[j] = y[j] + num.min((0,num.min(delta)))
        # debug arrays
        if debug:
            p.append((y[j] + poly[plidx:pridx] + line))
            d.append(delta)

    # add back linbgr
    bgr = bgr + linbgr
    if debug:
        return (bgr,p,d,linbgr)
    else:
        return bgr

############################################################################
def plot_bgr(data,nbgr=0,width=0,pow=0.5,tangent=False,debug=False):
    """
    make a fancy background plot
    """
    #
    if debug:
        (bgr,p,d,l) = background(data,nbgr=nbgr,width=width,pow=pow,
                                 tangent=tangent,debug=debug)
    else:
        bgr = background(data,nbgr=nbgr,width=width,pow=pow,
                         tangent=tangent,debug=debug)
        
    # plot data and bgr
    pyplot.figure(1)
    pyplot.clf()
    npts = len(data)
    pyplot.subplot(1,1,1)
    pyplot.plot(data,'k-o',label='data')
    pyplot.plot(bgr,'r-*',label='bgr')
    pyplot.plot(data-bgr,'g-',label='data-bgr')
    pyplot.plot(num.zeros(npts),'k-')
    pyplot.legend(loc=2)

    # plot data and polynomials
    # and for each polynomial
    # plot the diff between data
    # and the polynomial
    if debug == False: return
    
    pyplot.figure(2)
    pyplot.clf()
    pyplot.subplot(2,1,1)
    pyplot.plot(data-l,'k-o')
    pyplot.subplot(2,1,2)
    pyplot.plot(data-l,'k-o')
    pyplot.plot(num.zeros(npts),'k-')
    pyplot.plot(bgr-l,'k--')

    for j in range(len(p)):
        n = len(p[j])
        if j < int(width):
            s = 0
        else:
            s = j - int(width)
        xx = range(s,s+n)
        pyplot.subplot(2,1,1)
        pyplot.plot(xx,p[j],'*-')
        pyplot.subplot(2,1,2)
        #pyplot.plot(xx,d[j],'*-')
        dd = num.min((0,num.min(d[j]))) 
        pyplot.plot(xx,p[j]+dd,'*-')
    
################################################################################
################################################################################
if __name__ == '__main__':
    from matplotlib import pyplot
    from mpcutils.mathutil import gauss
    # generate a curve
    npts = 35
    x = num.array(range(npts))
    g1  = gauss(x, npts/2., 5., 300)
    g2  = gauss(x, npts/2., 30., 100)
    r = num.random.normal(size=npts)
    r = 10.*r/num.max(r)
    y = (1.0*r+ 10.*x) + g1 + g2
    # plot bgr
    width=5
    plot_bgr(y,nbgr=3,width=width,pow=8.,tangent=True,debug=True)
    
