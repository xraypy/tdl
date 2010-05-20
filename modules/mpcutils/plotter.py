"""
Simple wrapper functions and classes
around some of matplotlib.pyplot

Authors/Modifications:
----------------------
Tom Trainor (tptrainor@alaska.edu)

"""
################################################################################

import sys
import types
import time
import numpy as num

################################################################################
def newplotter(*arg,**kw):
    """
    create a new plot window
    """
    from matplotlib import pyplot
    pyplot.figure()
    #pyplot.clf()
    if arg:
        plotter(*arg,**kw)
        
################################################################################
def plotter(x,y=None,fmt='k-',xerr=None,yerr=None,xscale=1.,yscale=1.,
            xnorm=False,ynorm=False,xoff=0.,yoff=0.,xcen=False,
            ycen=False,xlog=False,ylog=False,nr=1,nc=1,np=1,hold=True,**kw):
    """
    plotting with many options
    """
    from matplotlib import pyplot

    x = num.array(x)
    if len(x) == 0:
        print "No data"
        return
    
    if y == None:
        y = x
        x = num.arange(len(y))
    else:
        y = num.array(y)
        if len(y) == 0:
            print "No data"
            return

    # calc x and y
    if xcen == True:
        ave_x = (num.min(x) + num.max(x))/2.
    else:
        ave_x = 0.
    if ycen == True:
        ave_y = (num.min(y) + num.max(y))/2.
    else:
        ave_y = 0.
        
    if xnorm == True:
        xscale = xscale/num.max(x-ave_x)
    if ynorm == True:
        yscale = yscale/num.max(y-ave_y)
    try:
        x_plot = (x-ave_x)*xscale + xoff
        y_plot = (y-ave_y)*yscale + yoff
    except:
        print "Error with plot data"
        return

    # check hold and set subplot (double check this...)
    # if clf: pyplot.clf()
    #pyplot.hold(hold)
    pyplot.subplot(nr,nc,np)
    if hold == False:
        pyplot.cla()
    else:
        pyplot.hold(True)
    
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
            pyplot.errorbar(x_plot,y_plot,xerr=xerr_plot,yerr=yerr_plot,fmt=fmt,**kw)
        except:
            print 'Error bar plot failed'
    else:
        try:
            pyplot.plot(x_plot,y_plot,fmt,**kw)
        except:
            print 'Plot failed'

    # check logs
    if xlog:
        try:
            pyplot.semilogx()
        except:
            print 'X log failed'
    if ylog:
        try:
            pyplot.semilogy()
        except:
            print 'X log failed'
    
    return

################################################################################
class PlotClick:
    """
    Handle plot click's 
    """
    ##########################################################################
    def __init__(self,fig=None,verbose=False):
        from matplotlib import pyplot
        self.backend = pyplot.get_backend()
        #
        if self.backend == "WXAgg":
            import wx
            self.GUI = wx
        elif self.backend == "TkAgg":
            self.GUI = pyplot.plot_root
        else:
            self.GUI = None
        #
        self.verbose = verbose
        self.clicked = False
        self.x = 0.
        self.y = 0.
        self.subplot=-1
        self.zoom=[]
        if fig == None:
            self.fig = pyplot.figure()
        else:
            self.fig = pyplot.figure(fig)
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

    ##########################################################################
    def _disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cid)

    ##########################################################################
    def __repr__(self):
        l = 'X=%g, Y=%g' % (self.x, self.y)
        return l

    ##########################################################################
    def on_click(self,event):
        (x, y) = (event.x, event.y)
        if event.button==1:
            if event.inaxes is not None:
                self.clicked = True
                self.x = event.xdata
                self.y = event.ydata
                #self.axes = event.inaxes
                self._subplot(event.inaxes)
                self._zoom()
                if self.verbose:
                    a = self.fig.gca()
                    xlbl = a.xaxis.label.get_text()
                    if len(xlbl) == 0: xlbl = 'X'
                    ylbl = a.yaxis.label.get_text()
                    if len(ylbl) == 0: ylbl = 'Y'
                    #print 'X=', event.xdata, 'Y=', event.ydata
                    #print 'X=', self.x, 'Y=', self.y
                    #print 'Subplot=', self.subplot, ', Zoom=',self.zoom
                    print "%s=%s, %s=%s" % (xlbl,str(self.x),ylbl,str(self.y))
                    
    ##########################################################################
    def _subplot(self,axes):
        a = self.fig.get_axes()
        n = len(a)
        for j in range(n):
            if axes == a[j]:
                self.fig.sca(a[j])
                self.subplot = j
                return
        self.subplot = 0
        #self.fig.sca(a[0])
        return

    ##########################################################################
    def _zoom(self):
        """
        get ((x1,y1),(x2,y2))
        http://matplotlib.sourceforge.net/examples/event_handling/lasso_demo.html
        Note could also just use zoom and then...
        #z = pyplot.axis()
        #z = [x1,x2,y1,y2]
        """
        a = self.fig.gca()
        x = a.get_xlim()
        y = a.get_ylim()
        self.zoom = [(x[0],y[0]),(x[1],y[1])]
        return
    
    ##########################################################################
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
        
################################################################################
def cursor(fig=None,verbose=False):
    """
    Allows interactive clicks...
    cursor()
    """
    #from matplotlib import pyplot
    #pyplot.connect('button_press_event', click.on_click)
    click = PlotClick(fig=fig,verbose=verbose)
    return click

################################################################################


