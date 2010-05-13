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
        # self.plotframe.reportLeftDown = self.customLeftDown
        self.id = id
        self.larch = larch
        if larch is not None:
            print 'setting Pllotter!'
            if not '_plotter' in larch.symtable:
                larch.symtable.newgroup('_plotter')
            larch.symtable.set_symbol('_plotter.win%i' % id, self)
        
    def onExit(self,o, **kw):
        if '_plotter' in self.larch.symtable:
            self.larch.symtable.del_symbol('_plotter.win%i' % self.id)
            
    def customLeftDown(self,event=None,**kw):
        print 'Plotter left down'
        if event == None: return                
        if self.larch is not None:
            if not '_plotter' in larch.symbtable:
                larch.symtable.newgroup('_plotter')
                larch.symtable.set_symbol('_plotter.win%i_x' % self.id, event.xdata)
                larch.symtable.set_symbol('_plotter.win%i_y' % self.id, event.ydata)
        
        self.plotframe.write_message("%f, %f" % (event.xdata,event.ydata), panel=1)


def _getPlotter(win=0, larch=None):
    if larch is None:
        return
    plotter = larch.symtable.get_symbol('_plotter.win%i' %win, create=True)
    print 'get1 ', plotter
    if plotter is None:
        Plotter(id=win, larch=larch)
        t0 = time.time()
        while time.time()-t0 < 10:
            try:
                plotter = larch.symtable.get_symbol('_plotter.win%i' %win)
            except:
                pass
            if plotter is not None:
                break
            time.sleep(0.5)
            print 'get... ', plotter
    print '_getPlotter - ', plotter
            
    return plotter

def _plot(x,y, win=0, larch=None, **kws):
    print 'Plot!! ', win, larch
    plotter = _getPlotter(win=win, larch=larch)
    print 'plotter ', plotter
    if plotter is not None:
        plotter.plot(x, y, **kws)    
    
def _oplot(x,y, win=0, larch=None, **kws):
    plotter = _getPlotter(win=win, larch=larch)
    if plotter is not None:
        plotter.oplot(x, y, **kws)
    
def register():
    return '_plotter', {'plot':_plot, 'oplot': _oplot}

        
