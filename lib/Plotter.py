#!/usr/bin/env python
##
## A simple plotting frame, wrapping Tk and matplotlib


import tdl.TkPlotter as TkPlotter
import Tkinter as Tk
import Pmw
import types

plot_group = '_plot'
plot_obj   = 'plotter'
_plotter   = "%s.%s" % (plot_group,plot_obj)
class Plotter:
    plot_attribs = {'title':'set_title',
                    'xlabel':'set_xlabel',
                    'ylabel':'set_ylabel'}
    def __init__(self, tdl, root=None,
                 exit_callback=None,
                 cursor_callback=None,
                 **kwds):
        self.plotter= None
        
        self.tdl  = tdl
        self.root = root
        self.kwds = kwds
        self.exit_callback = exit_callback
        self.cursor_callback = cursor_callback
        self.symtable = self.tdl.symbolTable
        self.symtable.addVariable(_plotter,value=self,constant=True)
        
        self.plotter= TkPlotter.PlotFrame(self.root,
                                          exit_callback=self.onExit,
                                          cursor_callback=self.onCursor,
                                          **kwds)
        self.root = self.plotter.root

            
    def onCursor(self,x=None,y=None,event=None):
        self.symtable.setVariable('%s.cursor_x'%plot_group,x)
        self.symtable.setVariable('%s.cursor_y'%plot_group,y)
        
    def onExit(self,event=None):
        self.plotter.destroy()
        self.plotter = self.root = None
        self.symtable.deleteSymbol(_plotter,override=True)


    def setPlotOptions(self):
        " sets plotting options based on tdl variables (all in the _plot group)"
        if self.tdl is None: return
        getval = self.symtable.getSymbolValue
        for name,method in self.plot_attribs.items():
            val = getval("%s.%s" % (plot_group,name),default=None)
            if val is not None:
                val = r"%s" % val
                getattr(self.plotter,method)(val)
    
    def plot(self,x,y=None,**kw):
        """plot after cleaing current plot """        
        self.setPlotOptions()
        self.plotter.plot(x,y=y,**kw)
        
    def oplot(self,x,y=None,**kw):
        """generic plotting method, overplotting any existing plot """
        self.setPlotOptions()
        self.plotter.oplot(x,y=y,**kw)

    def show_map(self,z, **kw):
        """generic plotting method, overplotting any existing plot """
        self.setPlotOptions()
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
        self.plotter.draw()

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



def _getPlot(tdl=None):
    if tdl is None: return None
    p = tdl.symbolTable.getSymbol(_plotter)
    if p is None or p.value is None:
        px = Plotter(tdl)
        p  = tdl.symbolTable.getSymbol(_plotter)
        if p is None: return None
    return p.value

def tdl_plot(x,y=None,new=False,tdl=None,**kw):
    p = _getPlot(tdl=tdl)
    if p is not None:
        pfunc = p.oplot
        if new: pfunc = p.plot
        # eval expression arguments here
        if type(x) is types.StringType: x = tdl.eval(x)
        if type(y) is types.StringType: y = tdl.eval(y)
        if x is not None:
            pfunc(x,y=y,**kw)
        else:
            print 'clear: ' 
            p.clear()
        
    else:       print 'cannot plot?'

def tdl_newplot(x,y=None,tdl=None,**kw):
    tdl_plot(x,y=y,tdl=tdl,new=True,**kw)

##########################
# todo:
#  -  control of legend, titles, etc  (tdl functions??)
#  -  draw box on zooming
#  -  GUI config

title = 'Tk Plotting library'

HelpPlot = """
  Simple Plotting in tdl:  there are two basic plotting commands:
     plot(x,y, ....)
  to add a trace to the existing plot, or the slight variation
     newplot(x,y, ....)
  which will clear the existing plot before plotting x and y

  There are several options to set line color, size, marker type, etc,
  that all come in keyword=value pairs  (say, color='red')
      option      meaning           example arguments
      ---------------------------------------------------------
      color       line color       'red', 'blue', '#EE00DD'
      style       line style       'solid','dashed', 'dotted'
      linewidth   width of line    1,2,3,4, 0 = 'no line'
      marker      plot symbol      '+','.','o','square',....
      markersize  size of marker   1,2,3,4, 0 = 'no marker'
      label       text label       'trace 1'
      yaxis       side for yaxis   'left', 'right'
      
  examples:
     x = arange(100)
     y = sin(x/5)
     z = cos(x/4)

     newplot(x,y,color='red')
     plot(x,z,color='black')
     plot(x,z,color='blue',style='dotted',linewidth=5)
     plot(x,z,color='green',linewidth=0,marker='+')           

   
"""

_help_ = {'plotting': HelpPlot}
_func_ = {'_builtin.plot':(tdl_plot, None),
          '_builtin.newplot':(tdl_newplot, None)}
