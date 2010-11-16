"""
Methods to handle generic background determination

Authors/Modifications:
----------------------
*  Tom Trainor (tptrainor@alaska.edu)
*  The polynomial background model follows closley 
   the XRF background code by Mark Rivers

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

    Parameters:
    -----------
    * data is the data to be fit.  We assume that the data is on an
      evenly spaced grid so no abscica is used in calculating the
      background line. If your data has uneven spacing then you
      should use a spline to regrid it before calling this routine
    * ngr is the number of end points to use in the fit.

    Example:
    --------
    >>bgr = linear_background(data,nbgr=3)
    
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
def background(data,nbgr=0,width=0,pow=0.5,tangent=False,
               compress=1,debug=False):
    """
    Calculate the background under a curve.

    Parameters:
    -----------
    * data is the data (ie y-data or ordinate) to be fit.
      We assume that the data is on an
      evenly spaced grid so no abscica is used in calculating the
      background line. If your data has uneven spacing then you
      should use a spline to regrid it before calling this routine

    * nbgr is the number of end points to use for calculation of a linear
      background.  This part of the background is removed from the data
      before polynomials are adjusted.  We then add it back to the polynomial
      sum so the total background returned from this function includes both
      if nbgr is zero no linear background is applied.  
      
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

    * compress is a compression factor for the array.  It should help
      speed up the fitting and smooth the background
      
    * debug is a flag (True/False) to indicate if additional debug arrays
      should be calculated
    
    Notes:
    ------
    This is a simplified form of the algorithm used for fitting
    XRF backgrounds (see Kajfosz and Kwiatek (1987) Nuc.
    Instr. Meth. in Phys. Res., B22, 78-81)

    At each data point (i) a calculated polynomial is brought up from underneath
    untill its first contact with any data point within the polynomial range
    (+/-width).  The distance between the apex of the polynomial (which occurs
    at point i) and zero is taken as the background for point i.

    This algorithm also allows for the inclusion of a linear background
    based on the end points.

    Note on computing the polynomial based on the width:
    Given the polynomial:
        poly    = (r**2. - delx**2)**pow  - (r**2.)**pow
    calc
        d^2(poly)/dx^2 = 0
    to determine where this function will have the max slope.
    The result is:
        delx_max_slope = r/sqrt(2*pow-1)
    If we assume that the 'width' argument corresponds to this
    delx_max_slope, then we can calc r for the polynomial as:
        r = width*sqrt(2*pow-1)
    If pow<=1/2 the we just use r = width

    Note we should rename pow, since pow is a builtin...
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
    if width <= 0. or pow == 0.: return linbgr
    y = data - linbgr

    # Compression
    # note we compress after linear bgr because
    # (1) linear fit is fast and (2) if we've
    # removed most of the slope then we limit
    # end point errors that result from trunction
    # if the array length is not an integer divisor
    # of the compression factor
    if compress > 1:
        (y,rem) = compress_array(y,compress)
        width = int(width/compress)
        if width == 0: width = 1 

    # create bgr array
    ndat = len(y)
    bgr  = num.zeros(ndat)
    
    # calc polynomial
    """
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
    """
    ## edits
    # make sure npoly is odd
    #npoly = num.min(int(6*width)+1, 2*int(ndat/2.)+1)
    #pdelx = num.array(range(npoly),dtype=float) - float((npoly-1.)/2.)
    npoly = num.min(int(3*width), int(ndat/2.))
    pdelx = range(-npoly,-1)
    pdelx.append(1)
    pdelx.extend(range(2,npoly+1))
    #print pdelx
    pdelx = num.array(pdelx,dtype=float)
    npoly = len(pdelx)
    r = float(width)
    poly = -1.*(pdelx/r)**(2.*pow)
    # renorm poly 
    pnorm = num.fabs(poly[int(npoly/2.) + int(width)])
    poly  = poly/pnorm
    ## end edits
    
    # loop through each point
    # NOTE this loop is the bottleneck
    # in the background calculations!
    # We should speed up this loop.  
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

    # do another linbgr to get residual
    # note this seem important since the polynomials
    # come up from the bottom, they will always tend to
    # underfit.  So having an additional linear part
    # helps to limit the residual. 
    linbgr2 = linear_background(y-bgr,nbgr=nbgr)
    bgr = bgr + linbgr2

    # Compression
    if compress > 1:
        bgr = expand_array(bgr,compress)
        if rem > 0:
            temp = bgr[-1]*num.ones(rem,dtype=bgr.dtype)
            bgr = num.append(bgr,temp)

    # Add back the original linear background / slope
    bgr = bgr + linbgr
    
    if debug:
        return (bgr,p,d,linbgr)
    else:
        return bgr

############################################################################
def plot_bgr(data,nbgr=0,width=0,pow=0.5,tangent=False,compress=1,debug=False):
    """
    make a fancy background plot
    """
    #
    if debug:
        (bgr,p,d,l) = background(data,nbgr=nbgr,width=width,pow=pow,
                                 tangent=tangent,compress=compress,debug=debug)
    else:
        bgr = background(data,nbgr=nbgr,width=width,pow=pow,
                         tangent=tangent,compress=compress,debug=debug)
        
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
        if j < int(n/2):
            s = 0
        else:
            s = j - int(n/2)
        xx = range(s,s+n)
        pyplot.subplot(2,1,1)
        pyplot.plot(xx,p[j],'*-')
        pyplot.subplot(2,1,2)
        #pyplot.plot(xx,d[j],'*-')
        dd = num.min((0,num.min(d[j]))) 
        pyplot.plot(xx,p[j]+dd,'*-')
    dma = num.max(data)
    pyplot.subplot(2,1,1)
    mi,ma = pyplot.ylim()
    mi = num.max((mi,-dma))
    pyplot.ylim(mi,1.1*dma)
    pyplot.subplot(2,1,2)
    mi,ma = pyplot.ylim()
    mi = num.max((mi,-dma))
    pyplot.ylim(mi,1.1*dma)
    
############################################################
def compress_array(array, compress):
    """
    Compresses a 1-D array by the integer factor "compress".

    Returns:
    --------
    * newarray is the compressed array
    * rem is the number of end points that were not
      compressed (ie compression has to be by integer divisor).

    Notes:
    ------
    Because compression needs to be an integer factor, not all
    points in the array will be compressed.  We simply resize the
    input array by chopping off the end so that its length becomes
    integer divisible by the compress factor.
    """
    compress = int(compress)
    alen = len(array)
    nlen = int(len(array)/compress)
    rem  = alen % compress

    temp = num.resize(array, (nlen, compress))
    newarray = num.sum(temp, 1)/compress
    #ra = array[alen-rem:]
    return (newarray,rem)

############################################################
def expand_array(array, expand, sample=0, rem=0):
    """
    Expands an 1-D array by the integer factor "expand".

    Parameters:
    -----------
    * array is the array to be expanded
    * expand is the expansion factor
    * sample is the sampling flag. if 'sample' is 1 the new array is
      created with sampling (ie no interpolation), if 0 then the new
      array is created via interpolation (default)
    * rem is not used... 
    """

    alen = len(array)
    if (expand == 1): return array
    if (sample == 1): return num.repeat(array, expand)

    kernel = num.ones(expand)/float(expand)
    temp = num.convolve(num.repeat(array, expand), kernel, mode=2)
    # Discard the first "expand-1" entries
    temp = temp[expand-1:]
    # Replace the last "expand" entries with the last entry of original
    for i in range(1,expand): temp[-i]=array[-1]
    if temp.dtype != array.dtype:
        temp = num.array(temp,dtype=array.dtype)
    return temp
    
################################################################################
################################################################################
if __name__ == '__main__':
    from matplotlib import pyplot
    from mpcutils.mathutil import gauss
    # generate a curve
    npts = 21
    x   = num.array(range(npts))
    g1  = gauss(x, npts/2., 5., 300)
    g2  = gauss(x, npts/2., 30., 100)
    r = num.random.normal(size=npts)
    r = 10.*r/num.max(r)
    y = (1.0*r+ 10.*x) + g1 + g2
    # plot bgr
    width=1
    plot_bgr(y,nbgr=3,width=width,pow=3.,tangent=False,
             compress=1,debug=True)
    
