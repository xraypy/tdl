"""
Scan Selector GUI
Author: Craig Biwer (cbiwer@uchicago.edu)
Last modified: 7.14.2011
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

from pds.pcgui.wxUtil import wxUtil
from pds.shellutil import mod_import
from tdl.modules import ana as scandata

rsrc_path = os.path.dirname(__file__)
#rsrc_path = os.path.split(os.path.abspath(__file__))[0]

############################################################################

class wxScanSelect(model.Background, wxUtil):

    ###########################################################
    # Init and util methods
    ###########################################################
    def on_initialize(self, event):
        # Initialization
        # including sizer setup, do it here
        # self.setupSizers()
        self.startup  = True
        self.dir      = '.'
        self.rsrc_path = rsrc_path

        # set up shell
        self.shell = None
        self.init_shell()

        # Make sure Scan Data is loaded
        self.exec_line("import ana")

        # init all the menus
        #self.init_ShellItems()
        self.init_GUI()
        self.on_Update_mouseClick(None)
        self.current = True

   ###########################################################

    def init_ShellItems(self,):
        #self.init_Grp()
        #self.init_Reader()
        #self.init_ScanVar()
        #self.init_BadMcas()
        #self.init_McaTaus()
        #self.init_XrfLines()
        pass

    def on_Update_mouseClick(self,event):
        self.init_ShellItems()
        self.init_SpecFile()
    
    ###########################################################
    #   Menu
    ###########################################################

    def on_menuFileExit_select(self,event): 
        self.close()
    
    def on_menuHelpDocumentation_select(self,event):
        self.exec_line("web 'http://cars9.uchicago.edu/ifeffit/tdl/Pds/SpecGui'")

    ######################################################

    ###########################################################
    #   Group / Reader
    ###########################################################

    ######################################################

    ######################################################
    # here update stuff from reader
    def init_GUI(self):
        if self.startup:
            self.startup = False
            cmap = copy.copy(pyplot.cm.cmapnames)
            cmap.insert(0,'')
            self.components.SpecPath.text = os.getcwd()

    ######################################################

    ###########################################################
    #   Path/files
    ###########################################################

    ######################################################
    def init_SpecPath(self):
        self.components.SpecPath.text = ''
        
    def on_SpecPathSel_mouseClick(self,event):        
        cdir = self.eval_line("pwd()")
        result = dialog.directoryDialog(self, 'Open', cdir)
        if result.accepted:
            dir = result.path
            #dir = dir.replace("\\","\\\\")
            self.components.SpecPath.text = dir
            self.init_SpecFile()
    
    def on_SpecPath_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN or keyCode == 372:
            if os.path.isdir(self.components.SpecPath.text.strip()):
                dir = self.components.SpecPath.text.strip()
                self.init_SpecFile()
        else:
            event.skip()
        return

    def init_SpecFile(self):
        dir = self.components.SpecPath.text.replace("\\","\\\\")
        dir = os.path.join(dir,'*.spc')
        files = glob.glob(dir)
        list = []
        for f in files:
            list.append(os.path.basename(f))
        sfile = self.components.SpecFile.stringSelection
        self.components.SpecFile.items = list
        if sfile in list:
            self.components.SpecFile.stringSelection=sfile
            
    def on_SpecFileSel_mouseClick(self,event):
        sdir = self.components.SpecPath.text
        cdir = self.eval_line("pwd()")
        if len(sdir.strip()) == 0:
            dir = cdir
        else:
            dir = sdir
        #
        result = dialog.fileDialog(self, 'Open', dir, '',"*")
        if result.accepted == False:
            self.post_message("File selection cancelled.")
            return
        #
        files = []
        #print result.paths
        for j in range(len(result.paths)):
            path        = result.paths[j]
            dir,fname   = os.path.split(path)
            if j == 0:
                if os.path.abspath(dir) != os.path.abspath(sdir):
                    self.components.SpecPath.text = dir
                    #self.init_SpecFile()
                    self.components.SpecFile.items = []
                fname0 = fname
            if fname not in files:
                files.append(fname)
        items = self.components.SpecFile.items
        for f in files:
            if f not in items: items.append(f)
        items.sort()
        self.components.SpecFile.items = items
        #
        self.components.SpecFile.stringSelection = fname0
        self.SpecFileSelect()
        #idx = items.index(fname0)
        #self.components.SpecFile.SetSelection(idx)
 
    ######################################################
    def on_SpecFile_select(self,event):
        self.SpecFileSelect()
        
    def SpecFileSelect(self):
        spath = str(self.components.SpecPath.text)
        sfile = str(self.components.SpecFile.stringSelection)
        #reader = self.get_Reader()
        #if reader == None: return None
        #if reader.spec_path != spath:
        #    reader.spec_path = os.path.abspath(spath)
        #reader.read_spec(sfile)
        #
        #self.init_ScanVar()
        #
        #for s in reader.spec_files:
        #    if s.fname == sfile:
        #        s.read()
        #        min = s.min_scan
        #        max = s.max_scan
        #        idx = num.arange(min,max+1,dtype=int)
        #        idx = map(str,idx)
        #        self.components.ScanNum.items = idx
        #        self.components.ScanNum.stringSelection=idx[-1]
        #        self.ScanNumSelect()

    ######################################################

    ###########################################################
    #   ScanVar/Num
    ###########################################################

    ######################################################

    #def on_ScanNum_select(self,event):
    #    self.ScanNumSelect()

    #def on_ScanNum_keyDown(self,event):
    #    keyCode = event.keyCode
    #    if keyCode == wx.WXK_RETURN:
    #        self.ScanNumSelect()
    #        self.ReadScan()
    #    else:
    #        event.skip()
    #    return

    #def ScanNumSelect(self,):
    #    sfile = str(self.components.SpecFile.stringSelection)
    #    snum = str(self.components.ScanNum.stringSelection)
    #    if self.components.AutoCalcVarName.checked:
    #        if self.components.LongName.checked:
    #            scan_var = u"%s.s%s" % (sfile,snum)
    #            scan_var = scan_var.replace('.','_')
    #        else:
    #            scan_var = u"s%s" % (snum)
    #        #
    #        items = self.components.ScanVar.items
    #        if scan_var not in items:
    #            items.append(scan_var)
    #            items.sort()
    #            self.components.ScanVar.items = items
    #        self.components.ScanVar.stringSelection = scan_var
    #    else:
    #        scan_var = self.components.ScanVar.stringSelection
    #        if  len(scan_var) == 0:
    #            self.components.ScanVar.text = 'tmp'
    #    if len(self.components.ScanVar.stringSelection) > 0:
    #        self.components.ScanVar.text = self.components.ScanVar.stringSelection
            
    #def init_ScanVar(self,):
    #    var = self.components.ScanVar.stringSelection
    #    tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False)
    #    self.components.ScanVar.items = tmp['var'] + tmp['ins']
    #    self.components.ScanVar.stringSelection = var
    #    return

    #def on_ScanVar_select(self,event):
    #    tmp = self.components.ScanVar.text
    #    if len(self.components.ScanVar.stringSelection) > 0:
    #        self.components.ScanVar.text = self.components.ScanVar.stringSelection
    #    var_name = self.components.ScanVar.text
    #    if len(var_name)>0:
    #        data = self.get_data(var_name)
    #        if (data != None) and (hasattr(data,'primary_axis')==True):
    #            self.UpdateGuiParmsFromData(data)
    #            if self.components.AutoUpdateCheck.checked==True:
    #                self.AutoPlot(var_name=var_name)

    #def on_Read_mouseClick(self,event):
    #    self.ReadScan()
        
    def on_Filter_mouseClick(self,event):
        if self.components.SpecFile.stringSelection=="" or self.components.SpecPath.text.strip()==" ":
            event.skip()
            return
        if not os.path.isfile(os.path.join(str(self.components.SpecPath.text), str(self.components.SpecFile.stringSelection))):
            event.skip()
            return
        from pds.pcgui import wxScanFilter
        wxScanFilter = mod_import(wxScanFilter)
        wxScanFilter = wxScanFilter.wxScanFilter
        filename = os.path.join(self.rsrc_path,'wxScanFilter.rsrc.py')
        self.wxScanFilter = model.childWindow(self, wxScanFilter,
                                            filename=filename)
        self.wxScanFilter.position = (320, 160)
        self.wxScanFilter.visible = True
        self.wxScanFilter.filename = os.path.join(str(self.components.SpecPath.text), str(self.components.SpecFile.stringSelection))

    ######################################################    
    #def ReadScan(self):
    #    ####
    #    self.UpdateReaderMedImgPar()

    #    ####
    #    fname = self.components.SpecFile.stringSelection
    #    if len(fname.strip())==0: return
    #    #
    #    scan_num = self.components.ScanNum.stringSelection
    #    if len(scan_num.strip())==0: return
    #    scan_num = int(scan_num)
    #    #
    #    var_name = self.components.ScanVar.text
    #    if len(var_name.strip())==0:
    #        #var_name = '__tmp__'
    #        print "Please enter a variable name"
    #        return 
    #    ####
    #    reader = self.get_Reader()
    #    if reader == None: return
    #    ###
    #    data = reader.spec_scan(scan_num,file=fname)
    #    self.set_data(var_name,data)
    #    ###
    #    
    #    self.UpdateGuiParmsFromData(data)        
    #    ####
    #    if self.components.AutoUpdateCheck.checked==True:
    #        self.AutoPlot(var_name=var_name)


##################################################
#if __name__ == '__main__':
#    app = model.Application(wxXrayData)
#    app.MainLoop()
