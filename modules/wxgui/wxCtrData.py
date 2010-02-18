"""
CTR Gui

Authors/Modifications:
----------------------
Tom Trainor (tptrainor@alaska.edu)


"""
############################################################################

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

########################################################################
CTR_HEADER = """Number of points       = %s
New/selected points    = %s
Current selected point = %s
Scan Type = '%s'
Labels: I='%s', Inorm='%s', Ierr='%s',Ibgr='%s'
Geom='%s'
Beam slits = %s
Det slits = %s
Sample  = %s
Scale = %s
H=%3.2f, K=%3.2f, L=%6.5f
I=%6.5g, Ierr=%6.5g, ctot=%6.5f
F=%6.5g, Ferr=%6.5g
"""

############################################################################

class wxCtrData(model.Background, wxUtil):

    ###########################################################
    # Init and util methods
    ###########################################################
    
    def on_initialize(self, event):
        # Initialization
        # including sizer setup, do it here
        # self.setupSizers()
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
        self.init_ShellItems()
        self.init_GUI()

    def init_ShellItems(self,):
        self.init_CtrDataVar()
        self.init_AppendScanVar()

    def on_Update_mouseClick(self,event):
        self.init_ShellItems()
    
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
        self.init_GUI()
        return

    def on_CtrDataVar_keyDown(self,event):
        "select a variable name and check it out"
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init_GUI()
        else:
            event.skip()
        return

    def get_CtrName(self,):
        if len(self.components.CtrDataVar.stringSelection) > 0:
            self.components.CtrDataVar.text = self.components.CtrDataVar.stringSelection
        ctr = self.components.CtrDataVar.text
        if len(ctr.strip()) == 0: return None
        name = "%s" % ctr.strip()
        return name

    def get_Ctr(self):
        name = self.get_CtrName()
        if name == None:
            print "Please provide a Ctr Data instance name"
            return None
        ctr = self.get_data(name)
        return ctr

    def set_Ctr(self,ctr):
        name = self.get_CtrName()
        return self.set_data(name,ctr)

    def CheckCtr(self,):
        try:
            name = self.get_CtrName()
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
        self.components.AppendScanVar.items = tmp
        if scan_var in tmp:
            self.components.AppendScanVar.stringSelection = scan_var
        else:
            self.components.AppendScanVar.text = ''
        return

    def get_Scan(self):
        if len(self.components.AppendScanVar.stringSelection) > 0:
            self.components.AppendScanVar.text = self.components.AppendScanVar.stringSelection
        scan = self.components.AppendScanVar.text
        if len(scan.strip()) == 0: return None
        name = "%s" % scan.strip()
        scan = self.get_data(name)
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
                return None
        except:
            self.post_message("Invalid scan object")
            return None

    ###########################################################
    #  Point Nums
    ###########################################################
    def init_PointNum(self):
        """
        Initialize the point nums. 
        """
        self.components.PointNum.items = ['0']
        self.components.PointNum.stringSelection = '0'
        self.point = 0
        self.set = []

    def init_AnchorPointNum(self):
        """
        Initialize the point nums. 
        """
        self.components.AnchorPointNum.items = ['0']
        self.components.AnchorPointNum.stringSelection = '0'
        
    ######################################################
    # here update stuff from ctr
    def init_GUI(self):
        if self.startup:
            self.startup = False
            cmap = copy.copy(pyplot.cm.cmapnames)
            cmap.insert(0,'')
            self.components.ColorMap.items = cmap 
            self.components.ColorMap.stringSelection = ''
        check = self.CheckCtr()
        if check == False:
            self.init_PointNum()
            self.components.PointNum.StringSelection = '0'
            return
        else:
            self.UpdateGuiFromCtr()

    def UpdateGuiFromCtr(self,):
        """
        update gui from a ctr data instance
        """
        ctr = self.get_Ctr()
        if ctr == None: return
        self.point = int(self.components.PointNum.StringSelection)
        point  = self.point
        set = self.set
        header   = CTR_HEADER % (str(len(ctr.L)),str(set),str(point),
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
                                 ctr.H[point],ctr.K[point],ctr.L[point],
                                 ctr.I[point],ctr.Ierr[point],ctr.ctot[point],
                                 ctr.F[point],ctr.Ferr[point])
        self.components.PointData.text = header
        #
        print self.get_Scan()
    
    ###########################################################
    #  Data Point
    ###########################################################

    def on_PointNum_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            pass
        else:
            event.skip()
        return

    ###########################################################
    #   Plot
    ###########################################################

    ######################################################
    def on_PlotImg_mouseClick(self,event):
        var_name = self.components.CtrDataVar.text
        self._plot_img(var_name)

    ######################################################
    def _plot_img(self,var_name):
        data = self.get_data(var_name)
        pnt = int(self.components.PointNum.stringSelection)
        """
        s = "%s.image_plot(%s.image.image[%s]" % (var_name,str(pnt))
        if self.components.ColorMap.stringSelection.strip()!='':
            cmap = self.components.ColorMap.stringSelection.strip()
            s = "%s,cmap='%s'" % (s,cmap)
        s = "%s)" % s
        #print s
        self.exec_line(s)
        """
        
##################################################
if __name__ == '__main__':
    app = model.Application(wxCtrData)
    app.MainLoop()
