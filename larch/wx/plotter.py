import MPlot
import sys
import wx
import time

class PlotDisplay(MPlot.PlotFrame):
    def __init__(self, id=0, larch=None, **kws):
        MPlot.PlotFrame.__init__(self, exit_callback=self.onExit, **kws)
        self.Show()
        self.Raise()
        self.cursor_pos = []
        self.panel.cursor_callback = self.onCursor
        self.id = id
        self.larch = larch
        if larch is not None:
            symtable = larch.symtable
            if not symtable.has_group('_plotter'):
                symtable.newgroup('_plotter')
            symtable.set_symbol('_plotter.plot%i' % id, self)
        
    def onExit(self,o, **kw):
        symtable = self.larch.symtable
        if symtable.has_group('_plotter'):
            try:
                symtable.del_symbol('_plotter.plot%i' % self.id)
            except:
                pass
        self.Destroy()
        
    def onCursor(self,x=None, y=None,**kw):
        symtable = self.larch.symtable
        if not symtable.has_group('_plotter'):
            symtable.newgroup('_plotter')
        symtable.set_symbol('_plotter.plot%i_x' % self.id, x)
        symtable.set_symbol('_plotter.plot%i_y' % self.id, y)
        

class ImageDisplay(MPlot.ImageFrame):
    def __init__(self, id=0, larch=None, **kws):
        MPlot.ImageFrame.__init__(self, exit_callback=self.onExit, **kws)
        self.Show()
        self.Raise()
        self.cursor_pos = []
        self.panel.cursor_callback = self.onCursor
        self.id = id
        self.larch = larch
        if larch is not None:
            symtable = larch.symtable
            if not symtable.has_group('_plotter'):
                symtable.newgroup('_plotter')
            symtable.set_symbol('_plotter.img%i' % id, self)
        
    def onExit(self,o, **kw):
        symtable = self.larch.symtable
        if symtable.has_group('_plotter'):
            try:
                symtable.del_symbol('_plotter.img%i' % self.id)
            except:
                pass
        self.Destroy()
        
    def onCursor(self,x=None, y=None,**kw):
        symtable = self.larch.symtable
        if not symtable.has_group('_plotter'):
            symtable.newgroup('_plotter')
        symtable.set_symbol('_plotter.img%i_x' % self.id, x)
        symtable.set_symbol('_plotter.img%i_y' % self.id, y)


def _getDisplay(win=0, larch=None, image=False):
    """make a plotter"""
    if larch is None:
        return
    symname = '_plotter.plot%i' %win
    title   = 'Larch Plot Display Window %i' % win
    creator = PlotDisplay
    if image:
        creator = ImageDisplay
    display = larch.symtable.get_symbol(symname, create=True)
    if display is None:
        display = creator(id=win, larch=larch)
        larch.symtable.set_symbol(symname, display)
        t0 = time.time()
        while display is None:
            display = larch.symtable.get_symbol(symname)
            time.sleep(0.05)
            if t0 - time.time() > 5.0:
                break
    if display is not None:
        display.set_title(title)
    return display

def _plot(x,y, win=0, larch=None, **kws):
    """plot doc"""
    plotter = _getDisplay(win=win, larch=larch)
    if plotter is not None:
        plotter.plot(x, y, **kws)    
    
def _oplot(x,y, win=0, larch=None, **kws):
    """oplot doc"""
    plotter = _getDisplay(win=win, larch=larch)
    if plotter is not None:
        plotter.oplot(x, y, **kws)

def _imshow(map, win=0, larch=None, **kws):
    """imshow doc"""
    img = _getDisplay(win=win, larch=larch, image=True)
    if img is not None:
        print img
        img.display(map, **kws)
    
def registerPlugin():
    return ('_plotter', True,
            {'plot':_plot,
             'oplot': _oplot,
             'imshow':_imshow})

        
