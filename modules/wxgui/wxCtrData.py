"""
CTR Gui

Authors/Modifications:
----------------------
Tom Trainor (tptrainor@alaska.edu)


"""
###############################################################################

from PythonCard import model, dialog
import wx
import os
import glob
import types
import copy
import time
import numpy as num
from matplotlib import pyplot

from wxUtil import wxUtil
import scandata
import ctr_data

###############################################################################

CTR_HEADER = """Number of points       = %s
New/selected points    = %s
Current selected point = %s
H=%3.2f, K=%3.2f, L=%6.5f
Scan Type = '%s'
Labels: I='%s', Inorm='%s', Ierr='%s', Ibgr='%s'
Geom='%s'
Beam slits = %s
Det slits = %s
Sample  = %s
Scale = %s
I=%6.5g, Ierr=%6.5g, Ibgr=%6.5g, ctot=%6.5f
F=%6.5g, Ferr=%6.5g
"""

###############################################################################

class wxCtrData(model.Background, wxUtil):

    ###########################################################
    # Init and util methods
    ###########################################################
    
    def on_initialize(self, event):
        """
        Initialization

        including sizer setup, do it here
        self.setupSizers()
        """
        self.startup  = True
        self.dir      = '.'
        self.point = 0
        self.set = []
        
        # set up shell
        self.shell = None
        self.init_shell()

        # Make sure Scan Data is loaded
        self.exec_line("import scandata")

        # init all the menus
        self.init_shell_items()
        self.init_gui()

    def init_shell_items(self,):
        self.init_CtrDataVar()
        self.init_AppendScanVar()

    def on_Update_mouseClick(self,event):
        self.init_shell_items()
    
    ###########################################################
    #   Menu
    ###########################################################

    def on_menuFileExit_select(self,event): 
        self.close()

    def on_menuHelpParams_select(self,event):
        pass

    ###########################################################
    #  Ctr Data
    ###########################################################

    def init_CtrDataVar(self):
        """
        Initialize the list of variables. 
        """
        ctr_var = self.components.CtrDataVar.stringSelection
        tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        #tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False,instance=scandata.CtrData)
        tmp = tmp['var'] + tmp['ins']
        tmp.sort()
        self.components.CtrDataVar.items = tmp
        if ctr_var in tmp:
            self.components.CtrDataVar.stringSelection = ctr_var
        else:
            self.components.CtrDataVar.text = ''
        return

    def on_CtrDataVar_select(self,event):
        "select a Ctr data name and check it out"
        self.init_gui()
        if self.components.AutoPlotCtr.checked==True:
            self._plot_ctr()
        return

    def on_CtrDataVar_keyDown(self,event):
        "select a variable name and check it out"
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init_gui()
            if self.components.AutoPlotCtr.checked==True:
                self._plot_ctr()
        else:
            event.skip()
        return

    def get_ctr_name(self,):
        if len(self.components.CtrDataVar.stringSelection) > 0:
            self.components.CtrDataVar.text = self.components.CtrDataVar.stringSelection
        ctr = self.components.CtrDataVar.text
        if len(ctr.strip()) == 0: return None
        name = "%s" % ctr.strip()
        return name

    def get_ctr(self):
        name = self.get_ctr_name()
        if name == None:
            print "Please provide a Ctr Data instance name"
            return None
        ctr = self.get_data(name)
        return ctr

    def set_ctr(self,ctr):
        name = self.get_ctr_name()
        return self.set_data(name,ctr)

    def check_ctr(self,):
        try:
            name = self.get_ctr_name()
            ctr  = self.get_data(name)
            #if isinstance(ctr,scandata.CtrData):
            if hasattr(ctr,'L'):
                self.post_message("Valid ctr object")
                return True
            else:
                self.post_message("Invalid ctr object")
                return False
        except:
            self.post_message("Invalid ctr object")
            return False

    ###########################################################
    #  Scan Data
    ###########################################################

    def init_AppendScanVar(self):
        """
        Initialize the list of variables. 
        """
        scan_var = self.components.AppendScanVar.stringSelection
        tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        #tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False,
        #                                                  instance=scandata.ScanData)
        tmp = tmp['var'] + tmp['ins']
        tmp.sort()
        tmp.insert(0,'')
        self.components.AppendScanVar.items = tmp
        if scan_var in tmp:
            self.components.AppendScanVar.stringSelection = scan_var
        else:
            self.components.AppendScanVar.stringSelection = ''
        return

    def on_AppendScanVar_select(self,event):
        """
        select a scan data name and check it out
        """
        s = self.get_scan()
        return

    def get_scan(self):
        scan = self.components.AppendScanVar.stringSelection
        if len(scan.strip()) == 0: return None
        name = "%s" % scan.strip()
        try:
            scan = self.get_data(name)
            #if isinstance(ctr,scandata.ScanData):
            if type(scan) != types.ListType:
                scan = [scan]
            if hasattr(scan[0],'get_positioner'):
                self.post_message("Valid scan object")
                return scan
            else:
                self.post_message("Invalid scan object")
                self.components.AppendScanVar.stringSelection = ''
                return None
        except:
            self.post_message("Invalid scan object")
            self.components.AppendScanVar.stringSelection = ''
            return None

    def on_AppendScan_mouseClick(self,event):
        """
        Select a scan data name and check it out
        """
        ctr = self.get_ctr()
        if ctr == None: return
        npts = len(ctr.L)
        scans = self.get_scan()
        if scans == None: return
        ctr.append_scans(scans)
        set  = range(npts,len(ctr.L))
        self.set = set
        self.update_gui_from_ctr()
        self.process_set()
        self.components.PointNum.stringSelection = str(npts)
        #
        self.components.AppendScanVar.stringSelection = ''
        return

    ###########################################################
    #  Point Nums
    ###########################################################
    
    def init_PointNum(self,pnts=None):
        """
        Initialize the point nums. 
        """
        if pnts == None:
            self.components.PointNum.items = ['0']
            self.components.PointNum.stringSelection = '0'
            self.point = 0
            self.set = []
        else:
            if len(self.components.PointNum.stringSelection) > 0:
                self.components.PointNum.text = self.components.PointNum.stringSelection
            try:
                point  = int(self.components.PointNum.text)
            except:
                print "Enter valid point"
                return
            if point not in pnts: point = 0
            self.components.PointNum.items = map(str,pnts)
            self.components.PointNum.stringSelection = str(point)
            self.components.PointNum.text = str(point)
            
    def init_AnchorPointNum(self,pnts=None):
        """
        Initialize the point nums. 
        """
        if pnts == None:
            self.components.AnchorPointNum.items = ['','0']
            self.components.AnchorPointNum.stringSelection = ''
        else:
            point= self.components.AnchorPointNum.stringSelection
            pnts = map(str,pnts)
            pnts.insert(0,'')
            if point not in pnts: point = ''
            self.components.AnchorPointNum.items = pnts
            self.components.AnchorPointNum.stringSelection = point

    def on_PointNum_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN or keyCode == 372:
            self.update_point()
        else:
            event.skip()
        return
    
    def on_PointNum_select(self,event):
        self.update_point()
        return

    def on_NextPoint_mouseClick(self,event):
        try:
            point  = int(self.components.PointNum.stringSelection)
        except:
            print "Enter valid point"
            return
        point = str(point + 1)
        if point in self.components.PointNum.items:
            self.components.PointNum.stringSelection = point
            self.update_point()
        return

    def on_PrevPoint_mouseClick(self,event):
        try:
            point  = int(self.components.PointNum.stringSelection)
        except:
            print "Enter valid point"
            return
        point = str(point - 1)
        if point in self.components.PointNum.items:
            self.components.PointNum.stringSelection = point
            self.update_point()
        return

    def on_PointNum_select(self,event):
        self.update_point()
        return

    def on_SetPointClick_mouseClick(self,event):
        ctr = self.get_ctr()
        if ctr == None: return
        point = ctr.get_idx()
        try:
            point = str(point[0])
            if point in self.components.PointNum.items:
                self.components.PointNum.stringSelection = point
                self.update_point()
        except:
            pass
        return

    def on_SetPointsClick_mouseClick(self,event):
        ctr = self.get_ctr()
        if ctr == None: return
        points = ctr.get_points()
        if points!=None:
            self.set = points
            self.update_gui_from_ctr()
            self.process_set()
        return

    def on_PointUpdate_mouseClick(self,event):
        self.update_point()

    def on_SetUpdate_mouseClick(self,event):
        self.process_set()

    def update_point(self,update_gui=True):
        ctr = self.get_ctr()
        if ctr == None: return
        apnt = self.components.AnchorPointNum.stringSelection
        if len(apnt) > 0:
            apnt = int(apnt)
            (intpar,corrpar) = ctr_data.get_params(ctr,apnt)
            if self.components.SetCorrParams.checked == False:
                corrpar = {}
            if self.components.SetIntParams.checked == False:
                intpar = {}
            point = int(self.components.PointNum.stringSelection)
            ctr_data.set_params(ctr,point,intpar=intpar,corrpar=corrpar)
        if self.components.AutoIntegrate.checked == True:
            self.integrate_point(update_plots=update_gui)
        if update_gui:
            self.update_gui_from_ctr()
            if self.components.AutoPlotScan.checked == True:
                self._plot_scan()

    def process_set(self):
        for point in self.set:
            self.components.PointNum.stringSelection = str(point)
            print "Process point ", point
            self.update_point(update_gui=False)
        self.components.PointNum.stringSelection = str(self.set[0])
        self.update_gui_from_ctr()
        if self.components.AutoPlotScan.checked == True:
            self._plot_scan()
    
    def integrate_point(self,update_plots=True):
        ctr = self.get_ctr()
        if ctr== None: return
        self.point = int(self.components.PointNum.text)
        point = self.point
        if update_plots == True:
            if self.components.AutoPlotIntegration.checked == True:
                plot = True
                fig = 3
            else:
                plot = False
                fig = None
            ctr.integrate_point(point,plot=plot,fig=fig)
            if self.components.AutoPlotCtr.checked==True:
                self._plot_ctr()
        else:
            ctr.integrate_point(point)
        
    ###########################################################
    #  Parameters
    ###########################################################
    
    def init_IntParamList(self,params=None):
        self.components.IntParam.text = ''
        self.components.IntParamDescr.text = ''
        if params == None:
            self.components.IntParamList.items = []
        else:
            items = []
            for (key,val) in params.items():
                items.append([key,str(val)])
            items.sort()
            self.components.IntParamList.items = items
    
    def init_CorrParamList(self,params=None):
        self.components.CorrParam.text = ''
        self.components.CorrParamDescr.text = ''
        if params == None:
            self.components.CorrParamList.items = []
            self.components.CorrParam.text = ''
        else:
            items = []
            for (key,val) in params.items():
                items.append([key,str(val)])
            items.sort()
            self.components.CorrParamList.items = items
        
    def on_IntParamList_select(self,event):
        selected = self.components.IntParamList.getStringSelection()
        (name,val) = selected[0]
        self.components.IntParam.text = val
        return

    def on_CorrParamList_select(self,event):
        selected = self.components.CorrParamList.getStringSelection()
        (name,val) = selected[0]
        self.components.CorrParam.text = val
        return

    ######################################################
    # update stuff from ctr etc
    ######################################################
    
    def init_gui(self):
        if self.startup:
            self.startup = False
            cmap = copy.copy(pyplot.cm.cmapnames)
            cmap.insert(0,'')
            self.components.ColorMap.items = cmap 
            self.components.ColorMap.stringSelection = ''
        check = self.check_ctr()
        if check == False:
            self.init_PointNum()
            self.init_AnchorPointNum()
            self.init_IntParamList()
            self.init_CorrParamList()
            # others...
            return
        else:
            self.update_gui_from_ctr()

    def update_gui_from_ctr(self,):
        """
        update gui from a ctr data instance
        """
        ctr = self.get_ctr()
        if ctr == None: return
        npts  = len(ctr.L)
        pnts  = range(npts)
        self.init_PointNum(pnts=pnts)
        self.init_AnchorPointNum(pnts=pnts)
        self.point = int(self.components.PointNum.text)
        #
        point = self.point
        set   = self.set
        header   = CTR_HEADER % (str(len(ctr.L)),str(set),str(point),
                                 ctr.H[point],ctr.K[point],ctr.L[point],
                                 str(ctr.scan_type[point]),
                                 str(ctr.labels['I'][point]),
                                 str(ctr.labels['Inorm'][point]),
                                 str(ctr.labels['Ierr'][point]),
                                 str(ctr.labels['Ibgr'][point]),
                                 str(ctr.corr_params[point].get('geom')),
                                 str(ctr.corr_params[point].get('beam_slits')),
                                 str(ctr.corr_params[point].get('det_slits')),
                                 str(ctr.corr_params[point].get('sample')),
                                 str(ctr.corr_params[point].get('scale')),
                                 ctr.I[point],ctr.Ierr[point],ctr.Ibgr[point],
                                 ctr.ctot[point],ctr.F[point],ctr.Ferr[point])
        self.components.PointData.text = header
        #
        (intpar,corrpar) = ctr_data.get_params(ctr,point)
        self.init_IntParamList(params=intpar)
        self.init_CorrParamList(params=corrpar)

    ###########################################################
    #   Plot
    ###########################################################

    ######################################################
    def on_PlotCtr_mouseClick(self,event):
        self._plot_ctr()

    def _plot_ctr(self):
        ctr = self.get_ctr()
        if ctr == None: return
        if self.components.CtrPlotInt.checked == True:
            ctr.plot_I(fig=0)
        else:
            ctr.plot(fig=0)
            
    ######################################################
    def on_PlotScanData_mouseClick(self,event):
        self._plot_scan()
        
    def _plot_scan(self):
        ctr = self.get_ctr()
        if ctr == None: return
        pnt = int(self.components.PointNum.stringSelection)
        if ctr.scan_type[pnt] == 'image':
            if self.components.ColorMap.stringSelection.strip()!='':
                cmap = self.components.ColorMap.stringSelection.strip()
                cmap = str(cmap)
            else:
                cmap=None
            ctr.plot_point(idx=pnt,fig=1,cmap=cmap)
        else:
            ctr.plot_point(idx=pnt,fig=1)

##################################################
if __name__ == '__main__':
    app = model.Application(wxCtrData)
    app.MainLoop()
