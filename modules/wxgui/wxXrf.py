########################################################################
"""
Tom Trainor (fftpt@uaf.edu)
This should be run as a child window from wxGUI

Modifications:
--------------

"""
########################################################################

from PythonCard import model, dialog
import wx
import os
import types
import numpy as num

from   wxUtil import wxUtil
import xrf_lookup
import xrf_data

#########################################################################

class wxXrf(model.Background, wxUtil):

    ###########################################################
    # Init and util methods
    ###########################################################
    def on_initialize(self, event):
        # Initialization
        # including sizer setup, do it here
        # self.setupSizers()
        self.components.PkParams._autoresize = 0
        self.is_scan        = False
        self.init_xrf_lines = True
        self.startup        = True
        self.dir            = '.'

        # set up shell
        self.shell = None
        self.init_shell()

        # Make sure Xrf is loaded
        self.exec_line("import xrf")

        # init all the menus
        self.init_shell_list_items()
        self.init_xrf()

    def init_shell_list_items(self):
        # the below lists are just updated
        # to include tdl variables ie list items
        # (should not change selections)
        self.init_GrpItems()
        self.init_NodeItems()
        self.init_SaveParVarItems()
        self.init_ExtractResultsVarItems()

    def init_xrf(self):
        #print "init"
        if self.startup:
            self.init_ScanParams()
            self.startup = False

        self.init_PlotParams()
        self.init_BgrParams()
        self.init_XrfLineItems()
        self.init_PkParams()
        self.init_FitParams()

        check = self.check_xrf_var()
        if check == False: return

        if self.components.SetParFromSave.checked:
            self.update_gui_from_save()
            self.update_xrf_from_gui()
        else:
            self.update_gui_from_xrf()

        # probably easiest to do this here if xrf var check ok
        if self.components.AutoUpdateCheck.checked:
            self.plot_cmd()
        return

    ###########################################################
    #             Menus                                       #
    ###########################################################
    def on_menuFileExit_select(self,event): 
        self.close()

    def on_menuFileReadXrf_select(self,event):
        
        result = dialog.fileDialog(self, 'Open...', self.dir, '',"*")
        if result.accepted:
            path      = result.paths[0]
            self.dir  = os.path.dirname(path)
        else:
            self.post_message("File selection cancelled.")
            return
        
        #print path, self.dir
        grp = self.components.Grp.text
        if len(grp)>0:
            var_name = "%s.xrf_data" % grp
        else:
            var_name = "xrf_data"

        line = "%s = []" % var_name  
        self.exec_line(line)
        for path in result.paths:
            path = path.replace("\\","\\\\")
            line = '__temp__ = xrf.read("%s")' % (path)
            self.exec_line(line)
            line = "%s.append(__temp__)" % (var_name)
            self.exec_line(line)
        self.exec_line("del __temp__")

        # set select list to default grp etc..
        # self.components.Grp.text  = 'xrf.data'
        self.components.Node.text = var_name
        self.init_shell_list_items()

        return

    def on_menuHelpParams_select(self,event): 
        import wxXrfHelp
        wxXrfHelp = mod_import(wxXrfHelp)
        dir       = os.path.dirname(wxXrfHelp.__file__)
        filename  = os.path.join(dir,'wxXrfHelp.rsrc.py')
        #print filename
        wxXrfHelp = wxXrfHelp.wxXrfHelp
        self.XrfHelpWindow = model.childWindow(self,wxXrfHelp,filename=filename)
        self.XrfHelpWindow.position = (200, 5)
        self.XrfHelpWindow.visible = True
    
    ###########################################################
    #       Update Variables              
    ###########################################################
    def on_UpdateShell_mouseClick(self, event):
        "use this to force update incase the data list has changed"
        self.init_shell_list_items()
        return

    ###########################################################
    #   Group and Node
    ###########################################################
    def get_xrf_var_name(self,ignore_idx=False):
        node = self.components.Node.text
        #idx = self.components.Node.selection
        #if idx > -1:
        #    if node != self.components.Node.items[idx]:
        #        node = self.components.Node.items[idx]
        #        self.components.Node.text = node
        #print "node:",node
        if len(node.strip()) == 0: return None

        idx  = self.components.ScanIdxSelect.stringSelection
        if len(idx.strip()) == 0:
            idx = None
        else:
            idx = int(idx)
        
        name = "%s" % node
        if (idx != None) and (self.is_scan == True) and (ignore_idx == False):
            name = "%s[%d]" % (name,idx)
        #print "xrf name", name, self.is_scan

        return name

    def get_xrf(self):
        name = self.get_xrf_var_name()
        if name == None: return None
        return self.get_data(name)

    def set_xrf(self,xrf):
        name = self.get_xrf_var_name()
        return self.set_data(name,xrf)

    def check_xrf_var(self):
        try:
            name = self.get_xrf_var_name(ignore_idx=True)
            m    = self.get_data(name)
            if hasattr(m,'xrf'):
                node = self.components.Node.text + '.xrf'
                self.components.Node.text = node
                m = m.xrf
            if type(m) in (types.ListType, num.ndarray):
                self.is_scan = True
                idx  = self.get_scan_idx()

                # set list items for xrf index
                tmp  = self.components.ScanIdxSelect.stringSelection
                self.components.ScanIdxSelect.items = idx
                if tmp in idx:
                    self.components.ScanIdxSelect.stringSelection = tmp
                else:
                    self.components.ScanIdxSelect.stringSelection = idx[0]

                # set list items for fit scan
                tmp = self.components.SaveParVar.text
                if len(tmp) > 0:
                    idx = [tmp] + [u'-1'] + idx
                else:
                    idx = [u'-1'] + idx
                self.components.ScanInitIdxSelect.items = idx
                self.components.ScanInitIdxSelect.stringSelection = u'-1'

                # get given scan
                name = self.get_xrf_var_name(ignore_idx=False)
                m    = self.get_data(name)
            else:
                self.is_scan = False
            self.ScanParamsToggle()

            if type(m) == types.InstanceType:
                if hasattr(m,'get_energy'):
                    self.post_message("Valid XRF object: %s" % name)
                    return True
                else:
                    self.post_message("Invalid XRF object: %s" % name)
                    return False
            else:
                self.post_message("Invalid XRF object")
                return False
            
        except:
            self.post_message("Invalid XRF object")
            return False

    def get_scan_idx(self,):
        scan_idx = []
        name = self.get_xrf_var_name(ignore_idx=True)
        m    = self.get_data(name)
        if type(m) in (types.ListType, num.ndarray):
            for k in range(len(m)):
                scan_idx.append(str(k))
            return scan_idx
        else:
            return [u'0']

    ######################################################
    def init_GrpItems(self):
        " Initialize the menu    "
        grp = self.components.Grp.text
        tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        self.components.Grp.items = tmp['ins']
        self.components.Grp.text = grp
        return

    def on_Grp_select(self, event):
        "Re-init Node list given the new grp name"
        grp = self.components.Grp.stringSelection
        self.components.Grp.text = grp
        self.init_NodeItems()
        return

    def on_Grp_keyDown(self,event):
        """
        select a variable name and check it out
        """
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init_NodeItems()
        else:
            event.skip()
        return

    def init_NodeItems(self):
        """Initialize the menu. Use the group thats
        selected in the group menu"""
        grp = self.components.Grp.text  
        if len(grp) == 0: grp = None
        node = self.components.Node.text
        tmp = self.shell.interp.symbol_table.list_symbols(symbol=grp,tunnel=False)
        tmp = tmp['var'] + tmp['ins']
        tmp.sort()
        self.components.Node.items = tmp
        #if node in tmp:
        self.components.Node.text = node
        #else:
        #    self.components.Node.text = ''
        return

    def on_Node_select(self,event):
        "select a variable name and check it out"
        self.init_xrf()
        return

    def on_Node_keyDown(self,event):
        "select a variable name and check it out"
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init_xrf()
        else:
            event.skip()
        return

    #def on_Node_textUpdate(self,event):
    #    "select a variable name and check it out"
    #    #print "update Node text"
    #    self.init_xrf()
    #    return
    
    ###########################################################
    # Scan Parameters
    ###########################################################
    def init_ScanParams(self):
        "Initialize the scan parameter data (except grp and node)"
        self.components.ScanIdxSelect.enabled = False
        self.components.ScanIdxSelect.text    = u'0'
        self.components.ScanIdxSelect.items   = [u'0'] 
        self.ScanParamsToggle()
        return
        
    def ScanParamsToggle(self):
        if self.is_scan == True:
            self.components.ScanIdxSelect.backgroundColor  = (255, 255, 255)
            self.components.ScanIdxSelect.enabled          = True
            #self.components.ScanIdxSelect.stringSelection  = u'0'
            self.components.FitScan.enabled                = True
        else:
            self.components.ScanIdxSelect.backgroundColor = (192, 192, 192)
            self.components.ScanIdxSelect.enabled         = False
            self.components.ScanIdxSelect.stringSelection = u'0'
            self.components.ScanIdxSelect.items           = [u'0'] 
            self.components.FitScan.enabled               = False
            self.components.ScanInitIdxSelect.items       = [u'-1'] 
            self.components.ScanInitIdxSelect.stringSelection = u'-1' 
    
    def on_ScanIdxSelect_textUpdate(self,event):
        self.init_xrf()

    def on_ScanIdxSelect_select(self,event):
        self.init_xrf()
    
    ###########################################################
    # Copy/Del/Update Parameters Buttons
    ###########################################################

    def on_ClearParamsAll_mouseClick(self,event):
        #self.init_PlotParams()
        self.init_BgrParams()
        self.init_PkParams()
        # keep energy offset and slope
        en_off = self.components.FitEnOffset.text
        en_slope = self.components.FitEnSlope.text 
        self.init_FitParams()
        self.components.FitEnOffset.text = en_off
        self.components.FitEnSlope.text = en_slope
        self.update_xrf_from_gui()
        pass
    
    ###########################################################
    #       Bgr Parameters                                    
    ###########################################################
    def init_BgrParams(self):
        self.components.BgrBtmWdth.text            = '4.0'
        self.components.BgrTopWdth.text            = '0.0'
        self.components.BgrTangent.stringSlection  = '0'
        self.components.BgrExp.text                = '2'
        self.components.BgrCompress.text           = '4'
        self.components.BgrBtmWdthFlag.stringSelection = '1'
        self.components.BgrTopWdthFlag.stringSelection = '1'
        self.on_BgrCheck_mouseClick(None)
        return
    
    def on_BgrCheck_mouseClick(self,event):
        if self.components.BgrCheck.checked == True:
            self.components.BgrBtmWdthFlag.editable = True
            self.components.BgrBtmWdthFlag.backgroundColor = (255, 255, 255)
            self.components.BgrTopWdthFlag.editable = True
            self.components.BgrTopWdthFlag.backgroundColor = (255, 255, 255)
        else:
            self.components.BgrBtmWdthFlag.editable = False
            self.components.BgrBtmWdthFlag.backgroundColor = (192, 192, 192)
            self.components.BgrTopWdthFlag.editable = False
            self.components.BgrTopWdthFlag.backgroundColor = (192, 192, 192)
        return

    def on_BgrDefault_mouseClick(self,event):
        self.components.BgrBtmWdth.text            = '4.0'
        self.components.BgrTopWdth.text            = '0.0'
        self.components.BgrTangent.stringSlection  = '0'
        self.components.BgrExp.text                = '2'
        self.components.BgrCompress.text           = '4'
        self.components.BgrBtmWdthFlag.stringSelection = '1'
        self.components.BgrTopWdthFlag.stringSelection = '1'
        return

    def get_BgrParams(self):
        """ return list of all Bgr parameters """
        return self.get_BgrPar_fields()

    def get_BgrPar_fields(self):
        "ret a list of all the peak parameter info in the entry fields"
        bgr_params = {}
        bgr_params['exponent']     = self.components.BgrExp.text.strip()
        bgr_params['top_width']    = self.components.BgrTopWdth.text.strip()
        bgr_params['bottom_width'] = self.components.BgrBtmWdth.text.strip()
        bgr_params['tangent']      = self.components.BgrTangent.stringSelection.strip()
        bgr_params['compress']     = self.components.BgrCompress.text.strip()
        bgr_params['bottom_width_flag'] = self.components.BgrBtmWdthFlag.stringSelection.strip()
        bgr_params['top_width_flag']    = self.components.BgrTopWdthFlag.stringSelection.strip()
        return bgr_params
    
    def put_BgrPar_fields(self,bgr_params):
        "reverse above"
        if bgr_params == None: return
        if len(bgr_params) < 5: return
        self.components.BgrExp.text       = str(bgr_params['exponent'])  
        self.components.BgrTopWdth.text   = str(bgr_params['top_width']) 
        self.components.BgrBtmWdth.text   = str(bgr_params['bottom_width']) 
        self.components.BgrTangent.text   = str(bgr_params['tangent']) 
        self.components.BgrCompress.text  = str(bgr_params['compress'])
        self.components.BgrBtmWdthFlag.stringSelection = str(bgr_params['bottom_width_flag'])
        self.components.BgrTopWdthFlag.stringSelection = str(bgr_params['top_width_flag'])
        return

    ###########################################################
    # XRF lines
    ###########################################################

    def init_XrfLineItems(self):
        " Initialize the menu    "
        if self.init_xrf_lines == True:
            if self.components.Emin:
                Emin = float(self.components.Emin.text)
            else:
                Emin = 0.01
            if self.components.Emax:
                Emax = float(self.components.Emax.text)
            else:
                Emax = 40.0
            #print Emin, Emax
            self.components.XrfLine.items = xrf_lookup.list_lines(Emin,Emax,short=True)
            self.init_xrf_lines = False
            #print "set to false"
    
    def on_XrfLine_select(self, event):
        "select an xrf line"
        if self.init_xrf_lines == True:
            self.init_XrfLineItems()
            #print "do update"
            return
        else:
            line = self.components.XrfLine.stringSelection
            en = xrf_lookup.lookup_xrf_line(line)
            if en:
                en = str(en)
            else:
                en = ''
            self.components.PkLbl.text                 = line
            self.components.PkEn.text                  = en
            self.components.PkAmp.text                 = '0.0'
            self.components.PkFWHM.text                = '0.0'
            self.components.PkEnFlag.stringSelection   = '1'
            self.components.PkFwhmFlag.stringSelection = '1'
            self.components.PkAmpFactor.text           = '0.0'
            self.components.PkIgnore.checked           = False
            return

    ###########################################################
    #       Peak Parameters                                   
    ###########################################################
    ## Note we also need the parameters in PkParams list to
    ## get updated on an enter event in any of the edit fields
    ## or mouseclick on sliders ...
    ## at the moment need to click add to update...
    ## see the event thing for energy fields

    def init_PkParams(self):
        self.components.PkParams.items = []
        self.components.PkLbl.text                 = ''
        self.components.PkEn.text                  = '0.0'
        self.components.PkAmp.text                 = '0.0'
        self.components.PkFWHM.text                = '0.0'
        self.components.PkEnFlag.stringSelection   = '1'
        self.components.PkFwhmFlag.stringSelection = '1'
        self.components.PkAmpFactor.text           = '0.0'
        self.components.PkIgnore.checked           = False
        return

    def _pk_par_dict(self,params):
        pk_par = {}
        pk_par['label']       = str(params[0]) 
        pk_par['energy']      = float(params[1]) 
        pk_par['ampl']        = float(params[2])   
        pk_par['fwhm']        = float(params[3])  
        pk_par['energy_flag'] = int(params[4])  
        pk_par['fwhm_flag']   = int(params[5])  
        pk_par['ampl_factor'] = float(params[6]) 
        pk_par['ignore']      = eval(str(params[7])) 
        pk_par['area']        = float(params[8])
        return pk_par

    def _pk_par_list(self,pk_par):
        params = ['']*9
        params[0] = str(pk_par['label'])   
        params[1] = str(pk_par['energy']) 
        params[2] = str(pk_par['ampl'])      
        params[3] = str(pk_par['fwhm'])      
        params[4] = str(pk_par['energy_flag'])  
        params[5] = str(pk_par['fwhm_flag'])    
        params[6] = str(pk_par['ampl_factor']) 
        params[7] = str(pk_par['ignore']) 
        params[8] = str(pk_par['area'])
        return params
    
    def on_PkParams_select(self,event):
        # print self.components.PkParams.items
        selected = self.components.PkParams.getStringSelection()
        params   = self._pk_par_dict(selected[0])
        self.put_PkPar_fields(params)
        return

    def get_PeakParams(self):
        """ return list of all Peak Parameters from GUI """
        params = self.components.PkParams.items
        pk_par = []
        for p in params:
            pk_par.append(self._pk_par_dict(p))
        return pk_par
    
    def PkParams_update(self,pk_params):        
        "update the pk params list"
        list   = self.components.PkParams.items
        params = self._pk_par_list(pk_params)
        
        if len(list) == 0:
            list.append(params)
            sel = 0
        else:
            found = False
            for j in range(len(list)):
                if list[j][0] == params[0]: 
                    list[j] = params
                    found = True
                    sel = j
                    break
            if found == False:
                list.append(params)
                sel = last = len(list)-1
        self.components.PkParams.items = list
        
        # this should invoke the on_PkParams_select event
        #self.components.PkParams.SetSelection(sel)

        return

    def get_PkPar_fields(self):
        "ret a list of all the peak parameter info in the entry fields"
        params = ['']*9
        # gather all the info
        lbl = self.components.PkLbl.text.strip()
        if len(lbl) < 1: return None
        params[0] = lbl
        params[1] = self.components.PkEn.text.strip()
        params[2] = self.components.PkAmp.text.strip()  
        params[3] = self.components.PkFWHM.text.strip() 
        params[4] = self.components.PkEnFlag.stringSelection.strip() 
        params[5] = self.components.PkFwhmFlag.stringSelection.strip() 
        params[6] = self.components.PkAmpFactor.text.strip()
        params[7] = str(self.components.PkIgnore.checked)
        params[8] = '0.0'
        pk_params = self._pk_par_dict(params) 
        return pk_params

    def put_PkPar_fields(self,pk_params):
        "reverse of above"
        if pk_params == None: return
        if len(pk_params) < 8: return
        params = self._pk_par_list(pk_params)
        #print params
        self.components.PkLbl.text                 = params[0]
        self.components.PkEn.text                  = params[1]
        self.components.PkAmp.text                 = params[2]
        self.components.PkFWHM.text                = params[3]
        self.components.PkEnFlag.stringSelection   = params[4]
        self.components.PkFwhmFlag.stringSelection = params[5]
        self.components.PkAmpFactor.text           = params[6]
        self.components.PkIgnore.checked           = eval(params[7])
        return

    ###########################################################
    def on_LineAdd_mouseClick(self,event):
        "add an xrf line"
        # get peak params from fields
        pk_params = self.get_PkPar_fields()
        if pk_params == None: return
        self.PkParams_update(pk_params)
        return

    def on_LineDel_mouseClick(self,event):
        selected =  self.components.PkParams.getStringSelection()
        if selected == None: return
        list = self.components.PkParams.items
        for sel in selected:
            j = 0
            while j < len(list):
                if sel[0] == list[j][0] and sel[1] == list[j][1]:
                    list.pop(j)
                    break
                j = j + 1
        self.components.PkParams.items = list

        # this should invoke the on_PkParams_select event
        self.components.PkParams.SetSelection(0)

        return
    
    ###########################################################
    # XRF Object
    ###########################################################

    def on_UpdateXRF_mouseClick(self,event):
        "update the xrf data variable given input data"
        self.update_xrf_from_gui()
        self.update_gui_from_xrf()
        if self.components.AutoUpdateCheck.checked:
            self.plot_cmd()
        return

    def compare_fit_params(self):
        """ return True if the xrf fit params are the same as
        the parameters listed in the GUI"""
        return False

    def get_params_from_gui(self):
        bgr_params = self.get_BgrParams()
        pk_params  = self.get_PeakParams()
        fit_params = self.get_FitParams()
        params = {'fit':fit_params,'bgr':bgr_params,'pk':pk_params}
        return params

    def update_xrf_from_gui(self):
        " update the xrf object from data in GUI"
        xrf = self.get_xrf()
        if xrf == None: return
        params = self.get_params_from_gui()
        guess  = self.components.PeakGuess.checked
        xrf.init(params=params,guess=guess,calc=True)
        # post updated xrf
        # do we need to do this?
        self.set_xrf(xrf)
        return

    def update_gui_from_xrf(self):
        xrf = self.get_xrf()
        if xrf == None:
            self.post_message("No XRF variable")
            return 
        try:
            # get parameters 
            p   = xrf.get_params()
            fit = p['fit']
            bgr = p['bgr']
            pk  = p['pk']

            # Fit params
            self.put_FitParams(fit)
            
            # BgrParams
            self.init_BgrParams()
            try:
                self.put_BgrPar_fields(bgr)
            except:
                pass

            #PeakParams
            self.init_PkParams()
            try:
                for j in range(len(pk)):
                    self.PkParams_update(pk[j])
            except:
                pass
            
            # Lines for extract results
            lst    = []
            for peak in pk:
                lst.append(peak['label'])
            self.components.ExtractLineSelect.items = lst

        except:
            self.post_message("Failed to read XRF object ")

        return

    ###########################################################
    # Save/Restore
    ###########################################################

    def init_SaveParVarItems(self):
        lst = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        lst = lst['var']
        t   = self.components.SaveParVar.text
        self.components.SaveParVar.items = lst
        self.components.SaveParVar.text  = t
        return

    def on_SaveSavePar_mouseClick(self,event):
        self.save_parameters()
        self.init_shell_list_items()
        return

    def on_RestoreSavePar_mouseClick(self,event):
        self.update_gui_from_save()
        return    

    def save_parameters(self):
        var_name = self.components.SaveParVar.text
        if len(var_name.strip()) == 0: return
        bgr_par  = self.get_BgrParams()
        peak_par = self.get_PeakParams()
        fit_par  = self.get_FitParams()

        xrf_params =  {'fit':fit_par,'bgr':bgr_par,'pk':peak_par}
        self.set_data(var_name,xrf_params)
        return

    def update_gui_from_save(self):
        #try:
        #    xrf = self.get_xrf() 
        #except:
        #    return

        # get save parameters
        var_name = self.components.SaveParVar.text
        if len(var_name.strip()) == 0: return
        xrf_params = self.get_data(var_name)

        #try:
        bgr_par = xrf_params['bgr']
        pk_par  = xrf_params['pk']
        fit_par = xrf_params['fit']
        
        # put stuff
        self.put_BgrPar_fields(bgr_par)
        for pk in pk_par:
            self.PkParams_update(pk)
        self.put_FitParams(fit_par)
        #except:
        #    self.post_message("failed to read save parameters")
        return
    
    ###########################################################
    # Calc/Fit
    ###########################################################

    def on_Calc_mouseClick(self,event):
        self.calc()
        return

    def calc(self):
        self.update_xrf_from_gui()
        xrf = self.get_xrf()
        xrf.calc()
        self.set_xrf(xrf)
        self.update_gui_from_xrf()
        if self.components.AutoUpdateCheck.checked == True:
            self.plot_cmd()
        return

    def init_FitParams(self):
        self.components.FitEnOffset.text            = '0.0'
        self.components.FitEnSlope.text             = '1.0'
        self.components.FitEnFlag.stringSelection   = '0'
        self.components.FitFwhmOffset.text          = '0.15'
        self.components.FitFwhmSlope.text           = '0.0'
        self.components.FitFwhmFlag.stringSelection = '0'
        self.components.FitChiExp.stringSelection   = '0.0'
        self.components.BgrCheck.checked            = False 

    def get_FitParams(self):
        fit_args = {}
        fit_args['energy_offset'] = self.components.FitEnOffset.text
        fit_args['energy_slope']  = self.components.FitEnSlope.text
        fit_args['energy_flag']   = self.components.FitEnFlag.stringSelection 
        fit_args['fwhm_offset']   = self.components.FitFwhmOffset.text
        fit_args['fwhm_slope']    = self.components.FitFwhmSlope.text
        fit_args['fwhm_flag']     = self.components.FitFwhmFlag.stringSelection
        fit_args['chi_exp']       = self.components.FitChiExp.stringSelection
        return fit_args

    def put_FitParams(self,fit_args):
        self.components.FitEnOffset.text            = str(fit_args['energy_offset']) 
        self.components.FitEnSlope.text             = str(fit_args['energy_slope'])
        self.components.FitEnFlag.stringSelection   = str(fit_args['energy_flag']) 
        self.components.FitFwhmOffset.text          = str(fit_args['fwhm_offset'])
        self.components.FitFwhmSlope.text           = str(fit_args['fwhm_slope'])
        self.components.FitFwhmFlag.stringSelection = str(fit_args['fwhm_flag']) 
        self.components.FitChiExp.stringSelection   = str(fit_args['chi_exp']) 
        return

    def on_Fit_mouseClick(self,event):
        self.fit()
        return

    def on_FitScan_mouseClick(self,event):
        self.fit_scan()
        return

    def fit(self):
        self.update_xrf_from_gui()
        xrf     = self.get_xrf()
        guess   = self.components.PeakGuess.checked
        opt_bgr = self.components.BgrCheck.checked 
        xrf.fit(guess = guess, opt_bgr=opt_bgr)
        
        # turn off guess after a fit!
        # other wise itll reguess if calc again
        # (e.g. plotting!)
        self.components.PeakGuess.checked = False

        #xrf.fit()
        self.set_xrf(xrf)
        self.update_gui_from_xrf()
        if self.components.AutoUpdateCheck.checked == True:
            self.plot_cmd()
        return

    def fit_scan(self):
        # must be a scan
        fit_init = self.components.ScanInitIdxSelect.stringSelection
        try:
            fit_init = int(fit_init)
            params   = {}
        except:
            params = self.get_data(fit_init)
            fit_init = -1

        use_prev = self.components.UsePrevFit.checked
        guess    = self.components.PeakGuess.checked
        name     = self.get_xrf_var_name(ignore_idx=True)
        m        = self.get_data(name)
        if type(m) != types.ListType:
            raise "Variable must be a list"
        xrf_data.fit(m,xrf_params=params,use_prev_fit=use_prev,
                    fit_init=fit_init,guess=guess)

        # turn off guess after a fit!
        # other wise itll reguess if calc again
        # (e.g. plotting!)
        self.components.PeakGuess.checked = False

        # update the data group if this is a scan
        self.update_scan_data_group()

        return

    def update_scan_data_group(self):
        try:
            grp = self.components.Grp.text
            if len(grp.strip())==0: return
            data = self.get_data(grp)
            if type(data) == types.InstanceType:
                if hasattr(data,'xrf'):
                    data.update_xrf_peaks()
        except:
            return

    ###########################################################
    # Extract Results
    ###########################################################
    def init_ExtractResultsVarItems(self):
        " Initialize the items from tdl    "
        # variables
        lst = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        lst = lst['var']
        t = self.components.ExtractResultsVar.text
        self.components.ExtractResultsVar.items = lst
        self.components.ExtractResultsVar.text  = t
        return

    def on_ExtractResults_mouseClick(self,event):

        res_name = self.components.ExtractResultsVar.text
        if len(res_name.strip()) == 0: return

        line = self.components.ExtractLineSelect.stringSelection
        if len(line.strip()) == 0: return

        xrf_name = self.get_xrf_var_name(ignore_idx=True)
        xrf = self.get_data(xrf_name)
        
        results = xrf_data.peak_areas(xrf,line)

        self.set_data(res_name,results)

        self.init_shell_list_items()

        return

    ###########################################################
    # Plotting
    ###########################################################

    def on_Plot_mouseClick(self,event):
        "make a plot"
        # update data and plot params
        self.update_xrf_from_gui()
        self.plot_cmd()
        self.update_gui_from_xrf()
        return

    def init_PlotParams(self):
        pass

    def plot_cmd(self):
        "run xrf.plot cmd"
        # update data and plot params
        #self.update_xrf_from_gui()
        if self.check_xrf_var() == False: return

        hold = self.components.HoldCheck.checked
        ylog = self.components.YlogCheck.checked
        xlog = self.components.XlogCheck.checked
        yerr = self.components.YerrCheck.checked

        data_check = self.components.DataPlotCheck.checked
        if data_check:
            data_str = self.components.DataFmt.stringSelection
        else:
            data_str = None

        fit_check = self.components.FitPlotCheck.checked
        if fit_check:
            fit_str  = self.components.FitFmt.stringSelection
        else:
            fit_str = None

        show_pk = self.components.ShowPeak.checked
        show_pk_bgr = self.components.ShowPeakBgr.checked
        if show_pk:
            if show_pk_bgr:
                pk_str = 'Peaks+Bgr'
            else:
                pk_str = 'Peaks'
        else:
            pk_str = None
            
        # build cmd
        xrf_name = self.get_xrf_var_name()
        xrf = self.get_data(xrf_name)
        xrf_data.xrf_plot(xrf,d=data_str,f=fit_str,p=pk_str,
                         ylog=ylog,xlog=xlog,hold=hold) 

        return


##################################################
if __name__ == '__main__':
    app = model.Application(wxXRF)
    app.MainLoop()
