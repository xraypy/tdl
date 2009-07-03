##########################################################################
"""
Tom Trainor (fftpt@uaf.edu)
Simple wrapper functions and classes
around some of matplotlib.pylab

Modifications:
--------------

"""
################################################################################

import sys
import types
import time
import numpy as Num

from   pds.util import PrintExceptErr

################################################################################
def newplotter(*arg,**kw):
    """
    create a new plot window
    """
    import pylab
    pylab.figure()
    #pylab.clf()
    if arg:
        plotter(*arg,**kw)
        
################################################################################
def plotter(x,y=None,fmt='k-',xerr=None,yerr=None,xscale=1.,yscale=1.,
            xnorm=False,ynorm=False,xoff=0.,yoff=0.,xcen=False,
            ycen=False,xlog=False,ylog=False,nr=1,nc=1,np=1,hold=True,**kw):
    """
    plotting with many options
    """
    import pylab

    x = Num.array(x)
    if len(x) == 0:
        print "No data"
        return
    
    if y == None:
        y = x
        x = Num.arange(len(y))
    else:
        y = Num.array(y)
        if len(y) == 0:
            print "No data"
            return

    # calc x and y
    if xcen == True:
        ave_x = (Num.min(x) + Num.max(x))/2.
    else:
        ave_x = 0.
    if ycen == True:
        ave_y = (Num.min(y) + Num.max(y))/2.
    else:
        ave_y = 0.
        
    if xnorm == True:
        xscale = xscale/Num.max(x-ave_x)
    if ynorm == True:
        yscale = yscale/Num.max(y-ave_y)
    try:
        x_plot = (x-ave_x)*xscale + xoff
        y_plot = (y-ave_y)*yscale + yoff
    except:
        print "Error with plot data"
        return

    # check hold and set subplot (double check this...)
    # if clf: pylab.clf()
    #pylab.hold(hold)
    pylab.subplot(nr,nc,np)
    if hold == False:
        pylab.cla()
    else:
        pylab.hold(True)
    
    # if errbars, make an errorbar plot
    if (xerr!=None) or (yerr != None):
        if xerr != None:
            xerr_plot = xerr*xscale
        else:
            xerr_plot = None
        if yerr != None:
            yerr_plot = yerr*yscale
        else:
            yerr_plot = None
        try:
            pylab.errorbar(x_plot,y_plot,xerr=xerr_plot,yerr=yerr_plot,fmt=fmt,**kw)
        except:
            PrintExceptErr('Error bar plot failed')
    else:
        try:
            pylab.plot(x_plot,y_plot,fmt,**kw)
        except:
            PrintExceptErr('Plot failed')

    # check logs
    if xlog:
        try:
            pylab.semilogx()
        except:
            PrintExceptErr('X log failed')
    if ylog:
        try:
            pylab.semilogy()
        except:
            PrintExceptErr('X log failed')
    
    return

################################################################################
class PlotClick:
    def __init__(self,fig=None,verbose=False):
        import pylab
        self.backend = pylab.get_backend()
        #
        if self.backend == "WXAgg":
            import wx
            self.GUI = wx
        elif self.backend == "TkAgg":
            self.GUI = pylab.plot_root
        else:
            self.GUI = None
        #
        self.verbose = verbose
        self.x = 0.
        self.y = 0.
        self.clicked = False
        if fig == None:
            self.fig = pylab.figure()
        else:
            self.fig = pylab.figure(fig)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)

    def __repr__(self):
        l = 'X=%g, Y=%g' % (self.x, self.y)
        return l

    def on_click(self,event):
        (x, y) = (event.x, event.y)
        if event.button==1:
            if event.inaxes is not None:
                if self.verbose:
                    print 'X=', event.xdata, 'Y=', event.ydata
                self.x = event.xdata
                self.y = event.ydata
                self.clicked = True

    def get_click(self,msg="Click a point"):
        """
        get the x and y coords
        """
        print msg
        self.clicked = False
        while self.clicked == False:
            time.sleep(.025)
            if self.backend == "WXAgg":
                #wx.YieldIfNeeded()
                self.GUI.YieldIfNeeded()
            elif self.backend == "TkAgg":
                self.GUI.update()
            else:
                print "Sorry dont know how to deal with this backend"
                return
        return (self.x,self.y)

    def get_lasso(self):
        """
        get ((x1,y1),(x2,y2)) lasso
        http://matplotlib.sourceforge.net/examples/event_handling/lasso_demo.html
        Note could also just use zoom and then...
        [x1,x2,y1,y2] = pylab.axis()
        """
        pass
        
################################################################################
def cursor(fig=None,verbose=False):
    """
    Allows interactive clicks...
    cursor()
    """
    #import pylab
    #pylab.connect('button_press_event', click.on_click)
    click = PlotClick(fig=fig,verbose=verbose)
    return click

################################################################################



