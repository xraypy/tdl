"""
This holds the items for menu->window
"""
#######################################################################

from PythonCard import model
import wx, string, sys, os, time
from wx import stc

from pds.lib.shellutil import mod_import

#######################################################################
class menuApps:
    
    def on_menuAppsPlotSelection_select(self, event):
        from wxPlotSelection import wxPlotSelection
        filename = os.path.join(self.rsrc_path,'wxPlotSelection.rsrc.py')
        self.PlotSelectionWindow = model.childWindow(self, wxPlotSelection,
                                                     filename=filename)
        self.PlotSelectionWindow.position = (200, 5)
        self.PlotSelectionWindow.visible = True

    """
    def on_menuAppsLTeQ_select(self, event):
        from LTEQ.wxLTEQ import vLTEQ
        filename = os.path.join(self.rsrc_path,'LTEQ','wxLTEQ.rsrc.py')
        self.PlotSelectionWindow = model.childWindow(self, vLTEQ,
                                                     filename=filename)
        self.PlotSelectionWindow.position = (200, 5)
        self.PlotSelectionWindow.visible = True                
    """
    
    def on_menuAppsXRF_select(self, event):
        import wxXrf
        wxXrf = mod_import(wxXrf)
        wxXrf = wxXrf.wxXrf
        filename = os.path.join(self.rsrc_path,'wxXrf.rsrc.py')
        self.XrfWindow = model.childWindow(self, wxXrf,filename=filename)
        self.XrfWindow.position = (200, 5)
        self.XrfWindow.visible = True

    def on_menuAppsXRRBuild_select(self, event):
        #__import__('wxXrr')
        #from wxXrr import wxXrr
        import wxXrrBuilder
        wxXrrBuilder = mod_import(wxXrrBuilder)
        wxXrrBuilder = wxXrrBuilder.wxXrrBuilder
        filename = os.path.join(self.rsrc_path,'wxXrrBuilder.rsrc.py')
        self.XrrBuildWindow = model.childWindow(self, wxXrrBuilder,
                                                filename=filename)
        self.XrrBuildWindow.position = (200, 5)
        self.XrrBuildWindow.visible = True

    def on_menuAppsXRRModel_select(self, event):
        #__import__('wxXrr')
        #from wxXrr import wxXrr
        import wxXrrIntModel
        wxXrrIntModel = mod_import(wxXrrIntModel)
        wxXrrIntModel = wxXrrIntModel.wxXrrIntModel
        filename = os.path.join(self.rsrc_path,'wxXrrIntModel.rsrc.py')
        self.XrrIntModelWindow = model.childWindow(self, wxXrrIntModel,
                                                   filename=filename)
        self.XrrIntModelWindow.position = (200, 5)
        self.XrrIntModelWindow.visible = True

    def on_menuAppsSpecData_select(self, event):
        import wxSpecData
        wxSpecData = mod_import(wxSpecData)
        wxSpecData = wxSpecData.wxSpecData
        filename = os.path.join(self.rsrc_path,'wxSpecData.rsrc.py')
        self.wxSpecData = model.childWindow(self, wxSpecData,
                                            filename=filename)
        self.wxSpecData.position = (200, 5)
        self.wxSpecData.visible = True

    def on_menuAppsCtrData_select(self, event):
        import wxCtrData
        wxCtrData = mod_import(wxCtrData)
        wxCtrData = wxCtrData.wxCtrData
        filename = os.path.join(self.rsrc_path,'wxCtrData.rsrc.py')
        self.wxCtrData = model.childWindow(self, wxCtrData,
                                           filename=filename)
        self.wxCtrData.position = (200, 5)
        self.wxCtrData.visible = True
