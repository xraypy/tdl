"""
This holds the items for menu->window
"""
#######################################################################

from PythonCard import model
import wx, string, sys, os, time
from wx import stc

from ..shellutil import mod_import

from .wxPlotSelection_rsrc import data as r_wxPlotSelection
from .wxXrf_rsrc           import data as r_wxXrf
from .wxXrrBuilder_rsrc    import data as r_wxXrrBuilder
from .wxXrrIntModel_rsrc   import data as r_wxXrrIntModel
from .wxSpecData_rsrc      import data as r_wxSpecData
from .wxCtrData_rsrc       import data as r_wxCtrData
from .wxScanSelect_rsrc    import data as r_wxScanSelect

from . import wxXrf, wxSpecData
from . import wxPlotSelection, wxCtrData
from . import wxFilter, wxScanSelect, wxIntegrator
from . import wxXrrBuilder, wxXrrIntModel

#######################################################################
def show_win(self, cls, rsrc):
    win = model.childWindow(self, cls, rsrc=rsrc)
    win.Position = (150, 50) # win.CenterOnScreen()
    win.visible = True
    return win

class menuApps:

    def on_menuAppsPlotSelection_select(self, event):
        cls = mod_import(wxPlotSelection).wxPlotSelection
        self.PlotSelectionWindow = show_win(self, cls, r_wxPlotSelection)

    def on_menuAppsXRF_select(self, event):
        cls = mod_import(wxXrf).wxXrf
        self.wxXrfWindow = show_win(self, cls, rsrc=r_wxXrf)

    def on_menuAppsXRRBuild_select(self, event):
        cls = mod_import(wxXrrBuilder).wxXrrBuilder
        self.wxXrrBuilder = show_win(self, cls, r_wxXrrBuilder)

    def on_menuAppsXRRModel_select(self, event):
        cls = mod_import(wxXrrIntModel).wxXrrIntModel
        self.wxXrrModl = show_win(self, cls, r_wxXrrIntModel)

    def on_menuAppsSpecData_select(self, event):
        cls = mod_import(wxSpecData).wxSpecData
        self.wxSpecData = show_win(self, cls, r_wxSpecData)

    def on_menuAppsCtrData_select(self, event):
        p = mod_import(wxCtrData).wxCtrData
        self.wxCtrData = show_win(self, p, r_wxCtrData)

    def on_menuAppsScanSelect_select(self, event):
        p = mod_import(wxScanSelect).wxScanSelect
        self.wxScanSelect = show_win(self, p, r_wxScanSelect)

    def on_menuAppsFilter_select(self, event):
        from . import wxFilter
        Filter = mod_import(wxFilter).filterGUI
        self.wxFilter = Filter(self)
        self.wxFilter.CenterOnScreen()

    def on_menuAppsIntegrator_select(self, event):
        from . import wxIntegrator
        Integrator = mod_import(wxIntegrator, True).Integrator
        self.wxIntegrator = Integrator(self)
        self.wxIntegrator.CenterOnScreen()
