#!/usr/bin/env python
##
## A simple plotting frame, wrapping Tk and matplotlib

import TkPlotter

class Plotter:
    def __init__(self, root=None, exit_callback=None, **kwds):
        self.plotter= None
        self.root = root
        self.kwds = kwds
        self.exit_callback = exit_callback

    def createPlotter(self,**kwds):
        if self.plotter is None:
            self.plotter= TkPlotter.PlotFrame(self.root,
                                              exit_callback=self.onExit,
                                              **kwds)
            self.root = self.plotter.root
            
    def onExit(self,event=None):
        self.plotter.destroy()
        self.plotter = self.root = None
        if self.exit_callback is not None: self.exit_callback()

        
    def plot(self,x,y=None,**kw):
        """plot after cleaing current plot """        
        self.plotter.plot(x,y=y,**kw)
        
    def oplot(self,x,y=None,**kw):
        """generic plotting method, overplotting any existing plot """
        self.plotter.oplot(x,y=y,**kw)

    def show_map(self,z, **kw):
        """generic plotting method, overplotting any existing plot """
        self.plotter.show_map(z,**kw)

    def update_line(self,t,x,y,**kw):
        """overwrite data for trace t """
        self.plotter.update_line(t,x,y,**kw)

    def set_xylims(self,xylims,**kw):
        """overwrite data for trace t """
        self.plotter.set_xylims(xylims,**kw)

    def get_xylims(self):
        """overwrite data for trace t """
        return self.plotter.get_xylims()

    def clear(self):
        """clear plot """
        self.plotter.clear()

    def unzoom_all(self,event=None):
        """zoom out full data range """
        self.plotter.unzoom_all(event=event)

    def unzoom(self,event=None):
        """zoom out 1 level, or to full data range """
        self.plotter.unzoom(event=event)
        
    def set_title(self,s):
        "set plot title"
        self.plotter.set_title(s)
        
    def set_xlabel(self,s):
        "set plot xlabel"        
        self.plotter.set_xlabel(s)

    def set_ylabel(self,s):
        "set plot xlabel"
        self.plotter.set_ylabel(s)        

    def save_figure(self,event=None):
        """ save figure image to file"""
        self.plotter.save_figure(event=event)

    def configure(self,event=None):
        self.plotter.configure(event=event)

    def setStatusText(self,text,**kw):
        self.plotter.setStatusText(text,**kw)



plotter_symbol = '_plot.plotter'
def _initPlot(tdl=None):
    if tdl is None: return None
    def _onClosePlot():
        tdl.symbolTable.deleteSymbol(plotter_symbol,override=True)

    p = Plotter(exit_callback=_onClosePlot)
    p.createPlotter()
    tdl.symbolTable.addVariable(plotter_symbol,value=p,constant=True)

def _getPlot(tdl=None):
    if tdl is None: return None
    p = tdl.symbolTable.getSymbol(plotter_symbol)
    if p is None or p.value is None:
        _initPlot(tdl=tdl)
        p = tdl.symbolTable.getSymbol(plotter_symbol)
        if p is None: return None

    return p.value

def tdl_plot(x,y=None,tdl=None,new=False,**kw):
    p = _getPlot(tdl=tdl)
    if p is None:
        print 'cannot plot?'
        return None
    plotfcn = p.plot
    if not new: plotfcn = p.oplot
    plotfcn(x,y=y,**kw)

def tdl_newplot(x,y=None,tdl=None,**kw):
    tdl_plot(x,y=y,tdl=None,new=True,**kw)

# functions to add to namespace
HelpPlot = """
  Plotting in TDL (simple plots)
"""

_help_ = {'plotting': HelpPlot}
_func_ = {'_builtin.plot':(tdl_plot, None),
          '_builtin.newplot':(tdl_newplot, None)}
