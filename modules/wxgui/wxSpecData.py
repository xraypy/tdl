############################################################################
"""
Tom Trainor (fftpt@uaf.edu)
Data / Spec GUI

Modifications:
--------------

"""
############################################################################

from PythonCard import model, dialog
import wx
import os
import glob
import types
import copy
import time
import numpy as Num

from   wxUtil import wxUtil
import scandata

############################################################################

class wxSpecData(model.Background, wxUtil):

    ###########################################################
    # Init and util methods
    ###########################################################
    def on_initialize(self, event):
        # Initialization
        # including sizer setup, do it here
        # self.setupSizers()
        self.startup  = True
        self.dir      = '.'

        # set up shell
        self.shell = None
        self.init_shell()

        # Make sure Scan Data is loaded
        self.exec_line("import scandata")

        # init all the menus
        self.init_ShellItems()
        self.init_GUI()
        self.current = True

   ###########################################################

    def init_ShellItems(self,):
        self.init_Grp()
        self.init_Reader()
        self.init_ScanVar()
        self.init_BadMcas()
        self.init_McaTaus()
        self.init_XrfLines()

    def on_Update_mouseClick(self,event):
        self.init_ShellItems()
        self.init_SpecFile()
    
    ###########################################################
    #   Menu
    ###########################################################

    def on_menuFileExit_select(self,event): 
        self.close()

    def on_menuHelpParams_select(self,event): 
        import wxXrayDataHelp
        wxXrayDataHelp = mod_import(wxXrayDataHelp)
        dir       = os.path.dirname(wxXrayDataHelp.__file__)
        filename  = os.path.join(dir,'wxXrayDataHelp.rsrc.py')
        wxXrayDataHelp = wxXrayDataHelp.wxXrayDataHelp
        self.wxXrayDataHelp = model.childWindow(self,wxXrayDataHelp,
                                                filename=filename)
        self.wxXrayDataHelp.position = (200, 5)
        self.wxXrayDataHelp.visible = True

    ######################################################

    ###########################################################
    #   Group / Reader
    ###########################################################

    ######################################################

    def init_Grp(self):
        " Initialize the menu    "
        grp = self.components.Grp.text
        tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        self.components.Grp.items = tmp['ins']
        self.components.Grp.text = grp
        return

    def on_Grp_select(self, event):
        "Re-init reader list given the new grp name"
        grp = self.components.Grp.stringSelection
        self.components.Grp.text = grp
        self.init_Model()
        return

    def on_Grp_keyDown(self,event):
        "select a variable name and check it out"
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init_Model()
        else:
            event.skip()
        return

    def init_Reader(self):
        """
        Initialize the menu. Use the group thats
        selected in the group menu
        """
        grp = self.components.Grp.text  
        if len(grp) == 0: grp = None
        reader = self.components.Reader.stringSelection
        tmp = self.shell.interp.symbol_table.list_symbols(symbol=grp,tunnel=False)
        tmp = tmp['var'] + tmp['ins']
        tmp.sort()
        self.components.Reader.items = tmp
        if reader in tmp:
            self.components.Reader.stringSelection = reader
        else:
            #self.components.Reader.stringSelection = 'reader'
            self.components.Reader.text = 'reader'
        return

    def on_Reader_select(self,event):
        "select a Reader name and check it out"
        self.init_GUI()
        return

    def on_Reader_keyDown(self,event):
        "select a variable name and check it out"
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init_GUI()
        else:
            event.skip()
        return

    def get_ReaderName(self,):
        if len(self.components.Reader.stringSelection) > 0:
            self.components.Reader.text = self.components.Reader.stringSelection
        reader = self.components.Reader.text
        if len(reader.strip()) == 0: return None
        name = "%s" % reader.strip()
        return name

    def get_Reader(self):
        name = self.get_ReaderName()
        if name == None:
            print "Please provide a reader name"
            return None
        reader = self.get_data(name)
        if reader == None:
            reader = scandata.Reader()
            self.set_Reader(reader)
        return reader

    def set_Reader(self,reader):
        name = self.get_ReaderName()
        return self.set_data(name,reader)

    def CheckReader(self,):
        try:
            name = self.get_ReaderName()
            r    = self.get_data(name)
            if type(r) == types.InstanceType:
                if hasattr(r,'read_spec'):
                    self.post_message("Valid reader object: %s" % name)
                    return True
                else:
                    self.post_message("Invalid reader object: %s" % name)
                    return False
            else:
                self.post_message("Invalid reader object")
                return False
        except:
            self.post_message("Invalid reader object")
            return False

    ######################################################
    # here update stuff from reader
    def init_GUI(self):
        if self.startup:
            self.startup = False
        check = self.CheckReader()
        if check == False:
            return
        else:
            self.UpdateGuiFromReader()

    def UpdateGuiFromReader(self,):
        """
        update gui from a reader
        """
        reader = self.get_Reader()
        if reader == None: return
        
        self.components.SpecPath.text = reader.spec_path
        #
        fname0 = ''
        items = []
        for j in range(len(reader.spec_files)):
            s = reader.spec_files[j].fname
            if j == 0: fname0 = s
            items.append(s)
        self.components.SpecFile.items = items
        if fname0:
            self.components.SpecFile.stringSelection = fname0
        #
        self.UpdateGuiMedImgPar(reader)
        #
        self.SpecFileSelect()

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
            dir = dir.replace("\\","\\\\")
            self.components.SpecPath.text = dir
            self.init_SpecFile()

    def init_SpecFile(self):
        dir = self.components.SpecPath.text
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
        reader = self.get_Reader()
        if reader == None: return None
        if reader.spec_path != spath:
            reader.spec_path = os.path.abspath(spath)
        reader.read_spec(sfile)
        #
        self.init_ScanVar()
        #
        for s in reader.spec_files:
            if s.fname == sfile:
                min = s.min_scan
                max = s.max_scan
                idx = Num.arange(min,max+1,dtype=int)
                idx = map(str,idx)
                self.components.ScanNum.items = idx
                self.components.ScanNum.stringSelection=idx[-1]
                self.ScanNumSelect()

    ######################################################

    ###########################################################
    #   ScanVar/Num
    ###########################################################

    ######################################################

    def on_ScanNum_select(self,event):
        self.ScanNumSelect()

    def on_ScanNum_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.ScanNumSelect()
            self.ReadScan()
        else:
            event.skip()
        return

    def ScanNumSelect(self,):
        sfile = str(self.components.SpecFile.stringSelection)
        snum = str(self.components.ScanNum.stringSelection)
        if self.components.AutoCalcVarName.checked:
            if self.components.LongName.checked:
                scan_var = u"%s.S%s" % (sfile,snum)
                scan_var = scan_var.replace('.','_')
            else:
                scan_var = u"S%s" % (snum)
            #
            items = self.components.ScanVar.items
            if scan_var not in items:
                items.append(scan_var)
                items.sort()
                self.components.ScanVar.items = items
            self.components.ScanVar.stringSelection = scan_var
        else:
            scan_var = self.components.ScanVar.stringSelection
            if  len(scan_var) == 0:
                self.components.ScanVar.text = 'tmp'
        if len(self.components.ScanVar.stringSelection) > 0:
            self.components.ScanVar.text = self.components.ScanVar.stringSelection
            
    def init_ScanVar(self,):
        var = self.components.ScanVar.stringSelection
        tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        self.components.ScanVar.items = tmp['var'] + tmp['ins']
        self.components.ScanVar.stringSelection = var
        return

    def on_ScanVar_select(self,event):
        tmp = self.components.ScanVar.text
        if len(self.components.ScanVar.stringSelection) > 0:
            self.components.ScanVar.text = self.components.ScanVar.stringSelection
        var_name = self.components.ScanVar.text
        if len(var_name)>0:
            data = self.get_data(var_name)
            if (data != None) and (hasattr(data,'primary_axis')==True):
                self.UpdateGuiParmsFromData(data)
                if self.components.AutoUpdateCheck.checked==True:
                    self.AutoPlot(var_name=var_name)

    def on_Read_mouseClick(self,event):
        self.ReadScan()

    ######################################################    
    def ReadScan(self):
        ####
        self.UpdateReaderMedImgPar()

        ####
        fname = self.components.SpecFile.stringSelection
        if len(fname.strip())==0: return
        #
        scan_num = self.components.ScanNum.stringSelection
        if len(scan_num.strip())==0: return
        scan_num = int(scan_num)
        #
        var_name = self.components.ScanVar.text
        if len(var_name.strip())==0:
            #var_name = '__tmp__'
            print "Please enter a variable name"
            return 
        ####
        reader = self.get_Reader()
        if reader == None: return
        ###
        data = reader.spec_scan(scan_num,file=fname)
        self.set_data(var_name,data)
        ###
        self.UpdateGuiParmsFromData(data)        
        ####
        if self.components.AutoUpdateCheck.checked==True:
            self.AutoPlot(var_name=var_name)

    def UpdateGuiParmsFromData(self,data):
        """
        update fields
        """
        ####
        x = data.positioners.keys()
        x.sort()
        x.insert(0,'')
        scalers = data.scalers.keys()
        scalers.sort()
        scalers.insert(0,'')
        xrf_lines = data.xrf_lines
        if xrf_lines == None: xrf_lines = []
        ####
        # scaler default
        default = self.components.DefaultScalerAxes.checked
        ####
        # Plot X
        tmp = self.components.XPlot.stringSelection
        self.components.XPlot.items = x + scalers
        if (len(tmp)>0) and (tmp in x) and (default==False):
            self.components.XPlot.stringSelection = tmp
        else:
            self.components.XPlot.stringSelection = data.primary_axis[0]
        ####
        # Plot Y
        tmp = self.components.YPlot.stringSelection
        self.components.YPlot.items = scalers + xrf_lines
        if (len(tmp)> 0) and (tmp in scalers) and (default==False):
            self.components.YPlot.stringSelection = tmp
        else:
            self.components.YPlot.stringSelection = data.primary_det
        ####
        # Plot norm
        tmp = self.components.NormPlot.stringSelection
        self.components.NormPlot.items = scalers
        if len(tmp)>0 and (tmp in scalers) and (default==False):
            self.components.NormPlot.stringSelection = tmp
        else:
            self.components.NormPlot.stringSelection = ''
        ####
        # Med Deadtime
        self.components.DTX.items = scalers + x 
        self.components.DTX.stringSelection = 'io'
        self.components.DTNorm.items = scalers
        self.components.DTNorm.stringSelection = 'Seconds'
        ###
        # Med
        npts = data.dims[0]
        pts = Num.arange(npts,dtype=int)
        pts = map(str,pts)
        self.components.ScanPnt.items = pts
        self.components.ScanPnt.stringSelection = '0'
        # check bad_idx
        tmp = str(self.components.BadMcas.text).strip()
        #print tmp
        if (len(tmp) == 0) and (len(data.med)>0):
            if (len(data.med[0].bad_mca_idx) > 0):
                self.components.BadMcas.text = repr(data.med[0].bad_mca_idx)
     
    ######################################################

    ###########################################################
    #   Med and Image Params
    ###########################################################

    ######################################################
    def init_XrfLines(self,):
        " Initialize the menu    "
        var = self.components.XrfLines.stringSelection
        tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        self.components.XrfLines.items = tmp['var']+ tmp['ins']
        self.components.XrfLines.stringSelection = var
        return

    def init_BadMcas(self,):
        " Initialize the menu    "
        var = self.components.BadMcas.stringSelection
        tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        self.components.BadMcas.items = tmp['var']+ tmp['ins']
        self.components.BadMcas.stringSelection = var
        return

    def init_McaTaus(self,):
        " Initialize the menu    "
        var = self.components.McaTaus.stringSelection
        tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        self.components.McaTaus.items = tmp['var'] + tmp['ins']
        self.components.McaTaus.stringSelection = var
        return

    def init_ImagePath(self):
        self.components.ImagePath.text = ''
 
    def on_MedPathSel_mouseClick(self,event):        
        cdir = self.eval_line("pwd()")
        result = dialog.directoryDialog(self, 'Open', cdir)
        if result.accepted:
            dir = result.path
            dir = dir.replace("\\","\\\\")
            self.components.MedPath.text = dir
       
    def on_ImgPathSel_mouseClick(self,event):        
        cdir = self.eval_line("pwd()")
        result = dialog.directoryDialog(self, 'Open', cdir)
        if result.accepted:
            dir = result.path
            dir = dir.replace("\\","\\\\")
            self.components.ImagePath.text = dir
            
    def UpdateGuiMedImgPar(self,reader):
        """
        update the GUI from reader
        """
        tmp = reader.med_path
        if tmp == None:
            self.components.MedPath.text = ''
        else:
            self.components.MedPath.text = str(tmp)
        self.components.ReadMed.checked = reader.med
        self.components.ReadXrf.checked = reader.xrf
        self.components.XrfLines.stringSelection = repr(reader.xrf_lines)
        self.components.BadMcas.stringSelection = repr(reader.med_params['bad_mca_idx'])
        self.components.McaTaus.stringSelection = repr(reader.med_params['tau'])
        self.components.Emin.text = repr(reader.med_params['emin'])
        self.components.Emax.text = repr(reader.med_params['emax'])
        self.components.Total.checked = reader.med_params['total']
        self.components.Align.checked = reader.med_params['align'] 
        self.components.CorrectData.checked = reader.med_params['correct']
        # missing fields for det_idx and nfmt       
        #
        tmp = reader.image_path
        if tmp == None:
            self.components.ImagePath.text = ''
        else:
            self.components.ImagePath.text = str(tmp)
        self.components.ReadImg.checked = reader.img
        
    def UpdateReaderMedImgPar(self):
        """
        update reader from GUI
        """
        reader = self.get_Reader()
        if reader == None: return
        #
        med_path = str(self.components.MedPath.text).strip()
        if len(med_path) > 0:
            reader.med_path=med_path
        else:
            reader.med_path=None
        #
        reader.med = self.components.ReadMed.checked
        reader.xrf = self.components.ReadXrf.checked
        #
        xrf_lines = str(self.components.XrfLines.stringSelection).strip()
        if len(xrf_lines) > 0:
            reader.xrf_lines = self.eval_line(xrf_lines)
        else:
            reader.xrf_lines = None
        #
        #bad_mcas  = str(self.components.BadMcas.stringSelection).strip()
        bad_mcas  = str(self.components.BadMcas.text).strip()
        if len(bad_mcas) > 0:
            reader.med_params['bad_mca_idx'] = self.eval_line(bad_mcas)
        else:
            reader.med_params['bad_mca_idx'] = []
        #
        mca_taus  = str(self.components.McaTaus.stringSelection).strip()
        if len(mca_taus) > 0:
            reader.med_params['tau'] = self.eval_line(mca_taus)
        else:
            reader.med_params['tau'] = None
        #
        emin = str(self.components.Emin.text).strip()
        if len(emin) > 0:
            reader.med_params['emin'] = self.eval_line(emin)
        else:
            reader.med_params['emin'] = -1.
        emax = str(self.components.Emax.text).strip()
        if len(emax) > 0:
            reader.med_params['emax'] = self.eval_line(emax)
        else:
            reader.med_params['emax'] = -1.
        #
        reader.med_params['total'] = self.components.Total.checked
        reader.med_params['align'] = self.components.Align.checked
        reader.med_params['correct'] = self.components.CorrectData.checked
        # missing fields for det_idx and nfmt       
        #
        image_path = str(self.components.ImagePath.text).strip()
        if len(image_path) > 0:
            reader.image_path=image_path
        else:
            reader.image_path=None
        reader.img = self.components.ReadImg.checked

    def on_FitDeadtime_mouseClick(self,event):
        #reader_name = self.get_ReaderName()
        tau_name = str(self.components.McaTaus.text).strip()
        if len(tau_name) == 0:
            print "Please give a 'Tau' variable name"
        var_name = self.components.ScanVar.text
        x = str(self.components.DTX.stringSelection)
        norm = str(self.components.DTNorm.stringSelection)
        #s = "%s.med_params['tau'] = scandata.fit_deadtime(%s,"
        #s = s % (reader_name, var_name)
        #s = s + "x='io',y='Med',norm='Seconds')"
        s = "%s = scandata.fit_deadtime(%s,"  % (tau_name, var_name)
        s = s + "x='%s',y='Med',norm='%s')" % (x, norm)
        #print s
        self.exec_line(s)

    ######################################################

    ###########################################################
    #   Plot
    ###########################################################

    ######################################################
    def on_PlotScaler_mouseClick(self,event):
        var_name = self.components.ScanVar.text
        self._plot_scaler(var_name)
        
    def on_PlotMed_mouseClick(self,event):
        var_name = self.components.ScanVar.text
        self._plot_med(var_name)
        
    def on_PlotImg_mouseClick(self,event):
        var_name = self.components.ScanVar.text
        self._plot_img(var_name)

    ######################################################
    def AutoPlot(self,var_name=None):
        if var_name == None:
            var_name = self.components.ScanVar.text
        if len(var_name.strip()) == 0:
            print "Please provide a scan variable name"
            return
        if self.components.ScalerPlot.checked:
            self._plot_scaler(var_name)
        if self.components.ReadMed.checked:
            self._plot_med(var_name)
        if self.components.ReadImg.checked:
            self._plot_img(var_name)

    ######################################################
    def _plot_scaler(self,var_name):
        import pylab
        pylab.figure(1)
        hold = str(self.components.HoldCheck.checked)
        xlog = str(self.components.XlogCheck.checked)
        ylog = str(self.components.YlogCheck.checked)
        #
        x = self.components.XPlot.stringSelection
        y = self.components.YPlot.stringSelection
        norm = self.components.NormPlot.stringSelection
        #
        if len(norm.strip()) > 0:
            s = "plot(%s['%s'],%s['%s']/%s['%s']" % (var_name,x,
                                                     var_name,y,
                                                     var_name,norm)
        else:
            s = "plot(%s['%s'],%s['%s']" % (var_name,x,
                                         var_name,y)
        #
        s = s + ",xlog=%s,ylog=%s,hold=%s)" % (xlog,ylog,hold)
        #print s
        self.exec_line(s)

    ######################################################
    def _plot_med(self,var_name):
        data = self.get_data(var_name)
        if data == None: return
        if len(data.med) == 0: return
        import pylab
        pylab.figure(2)
        hold = str(self.components.MedHold.checked)
        ylog = str(self.components.MedYlog.checked)
        pnt = int(self.components.ScanPnt.stringSelection)
        s = "scandata.med_plot(%s,scan_pnt=%s,hold=%s,ylog=%s)" % (var_name,
                                                                   str(pnt),
                                                                   hold,
                                                                   ylog)
        #print s
        self.exec_line(s)

    ######################################################
    def _plot_img(self,var_name):
        data = self.get_data(var_name)
        if data == None: return
        if len(data.image) == 0: return
        import pylab
        pylab.figure(3)
        pylab.clf()
        pnt = int(self.components.ScanPnt.stringSelection)
        s = "scandata.image_plot(%s.image[%s])" % (var_name,str(pnt))
        #print s
        self.exec_line(s)
        
##################################################
if __name__ == '__main__':
    app = model.Application(wxXrayData)
    app.MainLoop()
