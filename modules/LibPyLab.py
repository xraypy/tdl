############################################################################
# A very simple pylab wrapper for plotting.
# This is an interm solution for wx plotting till LibPlotter has
# a wx wrapper
#
#
############################################################################
# imports
# setup for matplotlib with Tk as default
import sys
import matplotlib
from Num import Num

# default to Tk, should work from system shell
DEFAULT_BACKEND = "TkAgg"
ISINIT  = False
PLOT_ROOT = None

PLOTTER = "_plot.plotter"

####################################################################
def init_matplotlib(backend=None):
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
    
    txt = "INIT MATPLOTLIB, backend = %s\n" % backend
    sys.__stdout__.write(txt)
    import matplotlib
    matplotlib.use(backend)
    matplotlib.interactive(True)
    import matplotlib.pylab
    matplotlib.pylab.show._needmain=False
    if backend=="TkAgg":
        try:
            matplotlib.pylab.plot_root = PLOT_ROOT                
        except:
            print_except_err("Error in plotter")
            return  
    ISINIT = True
    return 

#######################################################
# get the plot obj 
#
def get_plot(tdl):
    if not tdl:
        print "No TDL"
        return(None)
       
    if tdl.symbolTable.hasSymbol(PLOTTER):
        try:
            p = tdl.symbolTable.getSymbolValue(PLOTTER)
            if p == None:
                new_plot(tdl)
                return tdl.symbolTable.getSymbolValue(PLOTTER)
            else:
                return p
        except:
            new_plot(tdl)
            return tdl.symbolTable.getSymbolValue(PLOTTER) 
    else:
        new_plot(tdl)
        return tdl.symbolTable.getSymbolValue(PLOTTER)

def new_plot(tdl):
    backend = tdl.symbolTable.getSymbolValue('_builtin.GUI')
    if ISINIT == False:
        init_matplotlib(backend=backend)
    p = plotter(backend)
    tdl.symbolTable.addVariable(PLOTTER,p,constant=True)
    return 

#######################################################
#
def close_plot(tdl):
    if tdl.symbolTable.hasSymbol(PLOTTER):
        p = tdl.symbolTable.getSymbol(PLOTTER)
        p.close()
        tdl.symbolTable.deleteSymbol(PLOTTER,override=True)
    else:
        return 

#######################################################
# helper class for plots 
#
class plotter:
    def __init__(self,backend):
        self.plt = None
        self.cfig = 1
        self.x_cur = None
        self.y_cur = None
        self.backend = backend
        self.root = PLOT_ROOT

    def plot(self,x,y,fmt='',ylog=False,xlog=False,over=False,figure=1):
        matplotlib.pylab.figure(figure)
        if over == False:
            matplotlib.pylab.clf()
            self.plt = matplotlib.pylab.subplot(1,1,1)
        else:
            matplotlib.pylab.figure(figure)
            self.plt = matplotlib.pylab.subplot(1,1,1)

        if ylog and xlog:
            self.plt.loglog(x,y,fmt)
        elif ylog:
            self.plt.semilogy(x,y,fmt)
        elif xlog:
            self.plt.semilogx(x,y,fmt)
        else:
            self.plt.plot(x,y,fmt)
        matplotlib.pylab.show()
        return
    
    ################################
    # get cursor values from plot
    # get the x and y coords, flip y from top to bottom
    def on_click(self,event):
        x, y = event.x, event.y
        #if event.button==1:
        if event.button==3:
            if event.inaxes is not None:
                #print 'data coords', event.xdata, event.ydata
                self.x_cur = event.xdata
                self.y_cur = event.ydata
                self.clicked = True
    ##############################
    def is_clicked(self):
        if self.clicked == True:
            self.clicked = False
            return(True)
        return(False)
    ##############################
    def get_cur(self):
        matplotlib.pylab.connect('button_press_event', self.on_click)
        self.clicked = False
        while self.is_clicked() == False:
            if self.backend == "TkAgg":
                self.root.update()
            elif self.backend == "WXAgg":
                time.sleep(.05)
                wx.YieldIfNeeded()
                #wx.UpdateUIEvent()
        return (self.x_cur, self.y_cur)
    ##############################
    def close(self):
        return matplotlib.pylab.close()


#######################################################
# plot function
#
def plot(x, y=None, fmt = '', ylog=False,xlog=False,over = True, norm_dat = None, norm_max = 0, tdl = None):
    """Plot data
       plot x,y,[fmt=fmt,norm=norm,max=max,over=True]
       x, y = variables or eval type expressions
       norm = variable for y normalization
       max  = normalize the max of y to specified value
       fmt  = matlab style format string"""

    p = get_plot(tdl)
    if p == None:
        print "No Plotter"
        return 

    if y == None: y = x

    if norm_dat:
        for j in range(len(y)):
            y[j] = y[j]/norm_dat[j]

    if norm_max != 0:
        m1 = max(y)
        c = norm_max/m1
        for j in range(len(y)):
            y[j] = y[j] * c

    p.plot(x,y,ylog=ylog,xlog=xlog,over=over,fmt=fmt)
    return

def newplot(x, y=None, fmt = '', ylog=False,xlog=False,norm_dat = None, norm_max = 0, tdl = None):
    plot(x=x, y=y, fmt = fmt,ylog=ylog,xlog=xlog, over = False, norm_dat = norm_dat,
         norm_max = norm_max, tdl = tdl)
    
############################################################    

"""
def cur_fun(ds_arg, arg_str, **kw):

    # parse the arg string based on comma delimiters
    opt = {'name':'cursor'}
    args = dsi.ParseArgStr(arg_str,nreq=0,opt=opt)
    
    # get the arguments
    var_name = opt['name']

    p = DSPlot.get_plot(ds_arg)
    if p == None:
        print "no plot"
        return(FAILURE)

    print "Select point"
    (x,y) = p.get_cur() 
    d = dsi.AddNode(ds_arg,{'x':x,'y':y},var_name)
        
    return(SUCCESS)
    
"""
## add to list
#t = """Cursor"""
#d = """Get coordinates from plot"""
#u = """cursor [name=cursor]"""
#######################################################################


######################################################################
# plot function
#
"""
def close_fun(ds_arg, arg_str, **kw):
    return(DSPlot.close_plot(ds_arg))

"""    

#######################################################################


_groups_ = [('_plot',True)]
_var_    = {'_plot.plotter':None}
_func_ = {'_plot.plot':(plot,None),
          '_plot.newplot':(newplot,None)}


