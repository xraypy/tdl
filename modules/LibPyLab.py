############################################################################
# Make pylab plotting functions available to tdl
# T. Trainor
############################################################################
# imports
import sys
import types
from Util import PrintExceptErr

# setup for matplotlib with Tk as default
DEFAULT_BACKEND = "TkAgg"
PLOT_ROOT = None

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

    #################################
    _func_ = {'pylab.axes':pylab.axes,
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