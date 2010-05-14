import MPlot
import sys
import wx
import time

class Plotter(MPlot.PlotFrame):
    def __init__(self, id=0, larch=None, **kws):
        MPlot.PlotFrame.__init__(self, exit_callback=self.onExit, **kws)
        self.Show()
        self.Raise()
        self.cursor_pos = []
        print dir(self)
        self.panel.cursor_callback = self.onCursor
        self.id = id
        self.larch = larch
        if larch is not None:
            symtable = larch.symtable
            if not symtable.has_group('_plotter'):
                symtable.newgroup('_plotter')
            symtable.set_symbol('_plotter.win%i' % id, self)
        
    def onExit(self,o, **kw):
        symtable = self.larch.symtable
        if symtable.has_group('_plotter'):
            try:
                symtable.del_symbol('_plotter.win%i' % self.id)
            except:
                pass
        self.Destroy()
        
    def onCursor(self,x=None, y=None,**kw):
        symtable = self.larch.symtable
        if not symtable.has_group('_plotter'):
            symtable.newgroup('_plotter')
        symtable.set_symbol('_plotter.win%i_x' % self.id, x)
        symtable.set_symbol('_plotter.win%i_y' % self.id, y)
        

def _getPlotter(win=0, larch=None):
    """make a plotter"""
    if larch is None:
        return
    plotter = larch.symtable.get_symbol('_plotter.win%i' %win, create=True)
    if plotter is None:
        plotter = Plotter(id=win, larch=larch)
        t0 = time.time()
        larch.symtable.set_symbol('_plotter.win%i' %win, plotter)
        try:
            plotter = larch.symtable.get_symbol('_plotter.win%i' %win)
        except:
            time.sleep(0.5)
            try:
                plotter = larch.symtable.get_symbol('_plotter.win%i' %win)
            except:
                pass
    return plotter

def _plot(x,y, win=0, larch=None, **kws):
    """plot doc"""
    plotter = _getPlotter(win=win, larch=larch)
    if plotter is not None:
        plotter.plot(x, y, **kws)    
    
def _oplot(x,y, win=0, larch=None, **kws):
    """oplot doc"""
    plotter = _getPlotter(win=win, larch=larch)
    if plotter is not None:
        plotter.oplot(x, y, **kws)
    
def register():
    return '_plotter', {'plot':_plot,
                        'make_plotter':_getPlotter,
                        'oplot': _oplot}

        
