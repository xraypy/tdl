############################################################################
# Make pylab plotting functions available to tdl
# T. Trainor
############################################################################
# imports
import sys
import types
from Util import PrintExceptErr
import numpy as Num

# setup for matplotlib with Tk as default
DEFAULT_BACKEND = "TkAgg"
PLOT_ROOT = None

#####################################################################################
def plotter(x,y=None,fmt='k-',xerr=None,yerr=None,xscale=1,yscale=1,xnorm=False,ynorm=False,
            xoff=0,yoff=0,xcen=False,ycen=False,xlog=False,ylog=False,nr=1,nc=1,np=1,hold=True,**kw):
    """
    plotting with many options
    """
    import pylab
    if y == None:
        y = x
        x = Num.arange(len(y))

    # calc x and y
    if xcen == True:
        ave_x = (Num.min(x) + Num.max(x))/2.
    else:
        ave_x = 0
    if ycen == True:
        ave_y = (Num.min(y) + Num.max(y))/2.
    else:
        ave_y = 0
        
    if xnorm == True:
        xscale = xscale/Num.max(x-ave_x)
    if ynorm == True:
        yscale = yscale/Num.max(y-ave_y)
        
    x_plot = (x-ave_x)*xscale + xoff
    y_plot = (y-ave_y)*yscale + yoff

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


########################################################
def on_click(event):
    # get the x and y coords, flip y from top to bottom
    x, y = event.x, event.y
    if event.button==1:
        if event.inaxes is not None:
            print 'data coords', event.xdata, event.ydata

#########################################################################
def _init_pylab(tdl):
    backend = tdl.symbolTable.getSymbolValue('_builtin.GUI')

    if backend == None:
        backend = DEFAULT_BACKEND
    try:
        if backend == "TkAgg":
            import Tkinter
            PLOT_ROOT = Tkinter.Tk()
            PLOT_ROOT.iconify()
        elif backend == "WXAgg":
            import wx
            PLOT_ROOT = None
    except:
        raise "Error assigning plot backend"
    
    txt = "    **INIT MATPLOTLIB, backend = %s\n" % backend
    sys.__stdout__.write(txt)
    
    import matplotlib
    matplotlib.use(backend)
    import matplotlib.pylab as pylab

    pylab.show._needmain=False
    #matplotlib.interactive(True)
    pylab.ion()

    if backend=="TkAgg":
        try:
            pylab.plot_root = PLOT_ROOT                
        except:
            PrintExceptErr("Error in plotter")
            return

    # use this to print clicks to screen
    def cursor():
        pylab.connect('button_press_event', on_click)
    
    #################################
    _func_ = {'pylab.axes':pylab.axes,
          'pylab.axis':pylab.axis,
          'pylab.bar':pylab.bar,
          'pylab.barh':pylab.barh,
          'pylab.boxplot':pylab.boxplot,
          'pylab.cla':pylab.cla,
          'pylab.clf':pylab.clf,
          'pylab.close':pylab.close,
          'pylab.colorbar':pylab.colorbar,
          'pylab.cohere':pylab.cohere,
          'pylab.csd':pylab.csd,
          'pylab.cursor':cursor,
          'pylab.draw':pylab.draw,
          'pylab.errorbar':pylab.errorbar,
          'pylab.figlegend':pylab.figlegend,
          'pylab.figtext':pylab.figtext,
          'pylab.figimage':pylab.figimage,
          'pylab.figure':pylab.figure,
          'pylab.fill':pylab.fill,
          'pylab.gca':pylab.gca,
          'pylab.gcf':pylab.gcf,
          'pylab.gci':pylab.gci,
          'pylab.get':pylab.get,
          'pylab.gray':pylab.gray,
          'pylab.grid':pylab.grid,
          'pylab.barh':pylab.barh,
          'pylab.jet':pylab.jet,
          'pylab.hist':pylab.hist,
          'pylab.hold':pylab.hold,
          'pylab.imread':pylab.imread,
          'pylab.imshow':pylab.imshow,
          'pylab.ioff':pylab.ioff,
          'pylab.ion':pylab.ion,
          'pylab.isinteractive':pylab.isinteractive,
          'pylab.legend':pylab.legend,
          'pylab.loglog':pylab.loglog,
          'pylab.quiver':pylab.quiver,
          'pylab.rc':pylab.rc,
          'pylab.pcolor':pylab.pcolor,
          'pylab.pcolormesh':pylab.pcolormesh,
          'pylab.plot':pylab.plot,
          'pylab.psd':pylab.psd,
          'pylab.savefig':pylab.savefig,
          'pylab.scatter':pylab.scatter,
          'pylab.setp':pylab.setp,
          'pylab.semilogx':pylab.semilogx,
          'pylab.semilogy':pylab.semilogy,
          'pylab.show':pylab.show,
          'pylab.specgram':pylab.specgram,
          'pylab.stem':pylab.stem,
          'pylab.subplot':pylab.subplot,
          'pylab.subplots_adjust':pylab.subplots_adjust,
          'pylab.table':pylab.table,
          'pylab.text':pylab.text,
          'pylab.title':pylab.title,
          'pylab.xlabel':pylab.xlabel,
          'pylab.xticks':pylab.xticks,
          'pylab.ylabel':pylab.ylabel,
          'pylab.yticks':pylab.yticks,
          'pylab.pie':pylab.pie,
          'pylab.polar':pylab.polar,
          'pylab.plotter':plotter}
    for nam,val in _func_.items():
        cmdOut = None
        asCmd  = True
        func   = val
        #print nam
        if type(val) == types.TupleType:
            func = val[0]
            if len(val) > 1: cmdOut = val[1]
            if len(val) > 2: asCmd  = val[3]
        x =tdl.symbolTable.addFunction(nam,func,cmd_out=cmdOut,as_cmd=asCmd)

    return

##########################################################################

_groups_ = [('pylab',True)]
#_var_    = {'pylab.var':None}
_func_ = {'pylab._init':_init_pylab}

# code to run on initialization (no args, but will get a 'tdl reference')
#_init_ = _init_pylab


"""
#import pylab
_func_ = {'pylab._init':_init_pylab,
          'pylab.axes':pylab.axes,
          'pylab.axis':pylab.axis,
          'pylab.bar':pylab.bar,
          'pylab.boxplot':pylab.boxplot,
          'pylab.cla':pylab.cla,
          'pylab.clf':pylab.clf,
          'pylab.close':pylab.close,
          'pylab.colorbar':pylab.colorbar,
          'pylab.cohere':pylab.cohere,
          'pylab.csd':pylab.csd,
          'pylab.draw':pylab.draw,
          'pylab.errorbar':pylab.errorbar,
          'pylab.figlegend':pylab.figlegend,
          'pylab.figtext':pylab.figtext,
          'pylab.figimage':pylab.figimage,
          'pylab.figure':pylab.figure,
          'pylab.fill':pylab.fill,
          'pylab.gca':pylab.gca,
          'pylab.gcf':pylab.gcf,
          'pylab.gci':pylab.gci,
          'pylab.get':pylab.get,
          'pylab.gray':pylab.gray,
          'pylab.barh':pylab.barh,
          'pylab.jet':pylab.jet,
          'pylab.hist':pylab.hist,
          'pylab.hold':pylab.hold,
          'pylab.imread':pylab.imread,
          'pylab.imshow':pylab.imshow,
          'pylab.ioff':pylab.ioff,
          'pylab.ion':pylab.ion,
          'pylab.isinteractive':pylab.isinteractive,
          'pylab.legend':pylab.legend,
          'pylab.loglog':pylab.loglog,
          'pylab.quiver':pylab.quiver,
          'pylab.rc':pylab.rc,
          'pylab.pcolor':pylab.pcolor,
          'pylab.pcolormesh':pylab.pcolormesh,
          'pylab.plot':pylab.plot,
          'pylab.psd':pylab.psd,
          'pylab.savefig':pylab.savefig,
          'pylab.scatter':pylab.scatter,
          'pylab.setp':pylab.setp,
          'pylab.semilogx':pylab.semilogx,
          'pylab.semilogy':pylab.semilogy,
          'pylab.show':pylab.show,
          'pylab.specgram':pylab.specgram,
          'pylab.stem':pylab.stem,
          'pylab.subplot':pylab.subplot,
          'pylab.subplots_adjust':pylab.subplots_adjust,
          'pylab.table':pylab.table,
          'pylab.text':pylab.text,
          'pylab.title':pylab.title,
          'pylab.xlabel':pylab.xlabel,
          'pylab.ylabel':pylab.ylabel,
          'pylab.pie':pylab.pie,
          'pylab.polar':pylab.polar}
"""