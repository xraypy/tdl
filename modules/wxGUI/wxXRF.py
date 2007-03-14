############################################################################
# Tom Trainor (fftpt@uaf.edu)
# This should be run as a child window from wxGUI
#
# --------------
# Modifications:
# --------------
#  Next to do...
# - why if run fit when only a single det is chose do we get back data
#   for all the dets???? (when total is false)
# - add output of results
# - add fit xrf scan
#
############################################################################

from PythonCard import model, dialog
import wx
import os
import xrf_lookup
from wxTdlUtil import wxTdlUtil

"""
XRF GUI
"""
############################################################################
#
#ScanPars = {'grp':'None','node':'None','scan':'False','sc_start':'0',
#           'sc_stop':'0','sc_inc':'1','sc_step':'0'}
#
#DataPars = {'total':'True','align':'True','correct':'True','bad_mcas':'[]',
#          'mca_taus':'[]','emin':'0.0','emax':'0.0'}
#
#PlotPars = {'plt_data':True,'plt_fit':False,'plt_auto_update':False,
#          'plt_hold':False,'plt_components':False,
#          'plt_xlog':False,'plt_ylog':False,'plt_yerr':False}
#
#BgrPars =[[DetSelect,BgrExp,BgrTopWdth,BgrBtmWdth,
#           BgrTangent,BgrCompress,BgrCheck] ... ]
#
#PeakPars = [[Det,Label,PkEn,PkAmp,PkFWHM,PkEnFlag,PkFwhmFlag,
#             PkAmpFactor,PkIgnore,PkArea] ... ]
#
#FitPars = [FitFWHMFlag,FitEnFlag,FitChiExp,InitScanIdx,UsePrevFit,FitBgr]
#
#############################################################################

class wxXRF(model.Background, wxTdlUtil):

    ###########################################################
    # Init and util methods
    ###########################################################
    def on_initialize(self, event):
        # Initialization
        # including sizer setup, do it here
        # self.setupSizers()
        self.components.PkParams._autoresize = 0
        self.components.BgrParams._autoresize = 0
        self.dir = '.'
        self.init_xrf_lines = True
        self.startup = True

        # set up tdl
        self.shell = None
        self.tdl = None
        self.init_tdl()

        # init all the menus
        self.init_tdl_list_items()
        self.init()

    def init_tdl_list_items(self):
        # the below lists are just updated
        # to include tdl variables ie list items
        # (should not change selections)
        self.init_GrpItems()
        self.init_NodePfxItems()
        self.init_SaveParVarItems()
        self.init_DetAndTauItems()
        self.init_ExtractResultsVarItems()

    def init(self):
        # parameters and data
        #self.post_message("Initialize")
        #print "init"

        if self.startup:
            self.init_ScanParams()
            self.startup = False

        self.init_DataParams()
        self.init_EnRange()
        self.init_XrfLineItems()

        self.init_PkParams()
        self.init_BgrParams()
        self.init_FitParams()
        self.init_PlotParams()

        # below updates depend on how we want it to work 
        if self.components.SetParFromSave.checked:
            self.update_gui_from_save()
            self.update_xrf_from_gui()
        elif self.check_xrf_var():
            self.update_gui_from_xrf()

        # probably easiest to do this here if xrf var check ok
        if self.components.AutoUpdateCheck.checked:
            self.plot_cmd()
        return

    def get_xrf_var_name(self,grp=None,node=None):
        if grp == None: grp = self.components.Grp.text  
        if node == None: node = self.components.NodePfx.text
        if len(grp.strip()) == 0:
            var_name = node
        else:
            var_name = "%s.%s" % (grp,node)
        return var_name

    def get_xrf(self):
        var_name = self.get_xrf_var_name()
        return self.getValue(var_name)

    def set_xrf(self,xrf):
        var_name = self.get_xrf_var_name()
        return self.setValue(var_name,xrf)

    def check_xrf_var(self):
        try:
            var_name = self.get_xrf_var_name()
            m = self.getValue(var_name)
            num_mca = m.med.n_detectors
            if num_mca > 0:
                #self.post_message("Valid XRF object")
                return True
        except:
            self.post_message("Invalid XRF object")
            return False

    ###########################################################
    #             Menus                                       #
    ###########################################################
    def on_menuFileExit_select(self,event): 
        self.close()
        
    def on_menuFileReadXrf_select(self,event):
        dir = self.dir
        result = dialog.fileDialog(self, 'Open...', dir, '',"*")
        if result.accepted:
            path     = result.paths[0]
            self.dir = os.path.dirname(path)
        else:
            self.post_message("File selection cancelled.")
            return 
        for path in result.paths:
            line = "xrf.read '%s'" % path
            self.execLine(line)
        
        # set select list to default grp etc..
        self.components.Grp.text = 'xrf_data'
        self.components.NodePfx.text = ''
        self.init_tdl_list_items()
        return

    ###########################################################
    #       Update Tdl Variables              
    ###########################################################

    def on_UpdateTdl_mouseClick(self, event):
        "use this to force update incase the data list has changed"
        self.init_tdl_list_items()
        return

    ###########################################################
    #   Group and Node
    ###########################################################

    def init_GrpItems(self):
        " Initialize the menu    "
        grp = self.components.Grp.text
        if len(grp) == 0:
            grp = self.getValue('_builtin.data_group')
        self.components.Grp.items = self.listGroups()
        self.components.Grp.text = grp
        return

    def on_Grp_select(self, event):
        "Re-init Node list given the new grp name"
        self.init_NodePfxItems()
        return

    def init_NodePfxItems(self):
        """Initialize the menu. Use the group thats
        selected in the group menu"""
        grp = self.components.Grp.text  
        if len(grp) == 0: grp = None
        node = self.components.NodePfx.text
        self.components.NodePfx.items = self.listDataGroup(grp)
        self.components.NodePfx.text = node
        return

    def on_NodePfx_select(self,event):
        "select a variable name and check it out"
        self.init()
        return

    #def on_NodePfx_textUpdate(self,event):
    #    "select a variable name and check it out"
    #    #print "update Node text"
    #    self.init()
    #    return

    def on_NodePfx_keyDown(self,event):
        "select a variable name and check it out"
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init()
        else:
            event.skip()
        return
    
    ###########################################################
    # Scan Parameters
    ###########################################################

    def init_ScanParams(self):
        " Initialize the scan parameter data (except grp and node)"
        self.components.ScanCheck.checked = False
        self.components.ScanNodePfx.text = ''
        self.components.ScanStart.text = '0' 
        self.components.ScanEnd.text = '0'
        self.components.ScanFmt.text = '3'
        self.components.ScanInc.text = '0'
        self.ScanParamsToggle()
        return

    def on_ScanCheck_mouseClick(self,event):
        self.ScanParamsToggle()
        
    def ScanParamsToggle(self):
        if self.components.ScanCheck.checked == True:
            self.components.ScanNodePfx.editable = True
            self.components.ScanNodePfx.backgroundColor = (255, 255, 255)
            self.components.ScanStart.editable = True 
            self.components.ScanStart.backgroundColor = (255, 255, 255)
            self.components.ScanEnd.editable = True
            self.components.ScanEnd.backgroundColor = (255, 255, 255)
            self.components.ScanFmt.editable = True
            self.components.ScanFmt.backgroundColor = (255, 255, 255)
            self.components.ScanInc.editable = True
            self.components.ScanInc.backgroundColor = (255, 255, 255)
            self.components.ScanDown.enabled = True
            self.components.ScanUp.enabled = True
            self.components.FitScan.enabled = True
        else:
            self.components.ScanNodePfx.editable = False
            self.components.ScanNodePfx.backgroundColor = (192, 192, 192)
            self.components.ScanStart.editable = False 
            self.components.ScanStart.backgroundColor = (192, 192, 192)
            self.components.ScanEnd.editable = False
            self.components.ScanEnd.backgroundColor = (192, 192, 192)
            self.components.ScanFmt.editable = False
            self.components.ScanFmt.backgroundColor = (192, 192, 192)
            self.components.ScanInc.editable = False
            self.components.ScanInc.backgroundColor = (192, 192, 192)
            self.components.ScanDown.enabled = False
            self.components.ScanUp.enabled = False
            self.components.FitScan.enabled = False
            
    def get_ScanPar(self):
        """ Return ScanPar from the info in GUI components
        Note all the data is stored as strings"""
        scan_par = {}
        scan = self.components.ScanCheck.checked
        scan_par['scan'] = str(scan)
        node_pfx = self.components.ScanNodePfx.text
        if len(node_pfx) == 0: node_pfx = ''
        scan_par['node_pfx'] = node_pfx
        start = self.components.ScanStart.text
        scan_par['start'] = str(int(start))
        end = self.components.ScanEnd.text
        scan_par['end'] = str(int(end))
        fmt = self.components.ScanFmt.text
        scan_par['fmt'] = str(int(fmt))
        inc = self.components.ScanInc.text
        scan_par['inc'] = str(int(inc))
        return scan_par

    def put_ScanPar(self,scan_par):
        self.components.ScanCheck.checked = eval(scan_par['scan'])
        self.components.ScanNodePfx.text = scan_par['node_pfx']
        self.components.ScanStart.text = scan_par['start'] 
        self.components.ScanEnd.text = scan_par['end']
        self.components.ScanFmt.text = scan_par['fmt']
        self.components.ScanInc.text = scan_par['inc']

    def on_ScanStart_textUpdate(self,event):
        try:
            st = int(self.components.ScanStart.text)
        except:
            print "need integer"

    def on_AutoIncStop_textUpdate(self,event):
        try:
            en = int(self.components.ScanEnd.text)
        except:
            print "need integer"

    def on_ScanDown_mouseClick(self,event):
        try:
            st = int(self.components.ScanStart.text)
            en = int(self.components.ScanEnd.text)
            cur = int(self.components.ScanInc.text)
            if cur - 1 < st:
                self.components.ScanInc.text = str(en)
            else:
                self.components.ScanInc.text = str(cur-1)
            self.ScanExec()
        except:
            print "error in scan params"

    def on_ScanUp_mouseClick(self,event):
        try:
            st = int(self.components.ScanStart.text)
            en = int(self.components.ScanEnd.text)
            cur = int(self.components.ScanInc.text)
            if cur + 1 > en:
                self.components.ScanInc.text = str(st)
            else:
                self.components.ScanInc.text = str(cur+1)
            self.ScanExec()
        except:
            print "error in scan params"

    def on_ScanInc_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.ScanExec()
        else:
            event.skip()
        return

    def ScanExec(self):
        if not self.components.ScanCheck.checked:
            return
        inc = int(self.components.ScanInc.text)
        node_name = self.BuildScanName(inc)
        self.components.NodePfx.text = node_name
        # execute
        self.init()
        #if self.components.AutoUpdateCheck.checked:
        #    self.plot_cmd()

    def BuildScanName(self,inc):
        try:
            grp = self.components.Grp.text
            node_pfx = self.components.ScanNodePfx.text
            nc = self.components.ScanFmt.text
            fmt = "%s" + "%" + nc + "." + nc + "d"
            node_name = fmt % (node_pfx,inc)
            self.post_message(node_name)
            return node_name
        except:
            self.post_message("error building var name")

    ###########################################################
    # Data Parameters
    ###########################################################

    def init_DataParams(self):
        " Initialize the data param fields (en range and mca/tau done seperate)    "
        self.components.Total.checked = True
        self.components.Align.checked = True
        self.components.CorrectData.checked = True
        self.components.NumMcas.text = "NumMcas = 0"
        self.components.BadMcas.text = '[]'
        self.components.McaTaus.text = ''
        return

    def init_DetAndTauItems(self):
        " Initialize the items from tdl    "
        lst = self.listAllData()
        t1 = self.components.BadMcas.text
        t2 = self.components.McaTaus.text
        self.components.BadMcas.items = lst
        self.components.McaTaus.items = lst
        self.components.BadMcas.text = t1
        self.components.McaTaus.text = t2
        return
    
    def init_EnRange(self):
        if self.components.Emin.text != '2.0':
            self.components.Emin.text = '2.0'
            self.init_xrf_lines == True
        if self.components.Emax.text != '20.0':
            self.components.Emax.text = '20.0'
            self.init_xrf_lines == True

    def on_EnRange_mouseClick(self,event):
        """ Select energy Range from plot"""
        self.plot_cmd()
        # --> select Emin/Emax
        self.set_xrf_data()
        self.init_xrf_lines = True
        return

    def on_Emin_textUpdate(self,event):
        """ new emin """
        self.init_xrf_lines = True

    def on_Emax_textUpdate(self,event):
        """ new emax """
        self.init_xrf_lines = True

    #def on_BadMcas_textUpdate(self,event):
    #    pass

    #def on_McaTaus_textUpdate(self,event):
    #    pass

    def get_DataPar(self):
        """ Return DataPar from the info in GUI components
        Note all the data is stored as strings"""
        data_par = {}
        total = self.components.Total.checked
        data_par['total'] = str(total)
        align = self.components.Align.checked
        data_par['align'] = str(align)
        correct = self.components.CorrectData.checked
        data_par['correct'] = str(correct)
        bad_mcas = self.components.BadMcas.text
        data_par['bad_mcas'] = bad_mcas
        mca_taus = self.components.McaTaus.text
        data_par['mca_taus'] = mca_taus
        Emin = self.components.Emin.text
        data_par['emin'] = Emin
        Emax = self.components.Emax.text
        data_par['emax'] = Emax
        return data_par

    def put_DataPar(self,data_par):
        self.components.Total.checked = eval(data_par['total'])
        self.components.Align.checked = eval(data_par['align'])
        self.components.CorrectData.checked = eval(data_par['correct'])
        self.components.BadMcas.text = data_par['bad_mcas']
        self.components.McaTaus.text = data_par['mca_taus'] 
        self.components.Emin.text = data_par['emin'] 
        self.components.Emax.text = data_par['emax']
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
            self.components.PkLbl.text       = line
            self.components.PkEn.text        = en
            self.components.PkAmp.text       = '0.0'
            self.components.PkFWHM.text      = '0.0'
            self.components.PkEnFlag.text    = '0'
            self.components.PkFwhmFlag.text  = '0'
            self.components.PkAmpFactor.text = '0.0'
            self.components.PkIgnore.checked = False
            return

    def on_XrfLineAdd_mouseClick(self,event):
        "add an xrf line"
        # get peak params from fields
        pk_params = self.get_PkPar_fields()
        if pk_params == None: return
        # update the list
        self.PkParams_update(pk_params)
        return

    def on_XrfLineDel_mouseClick(self,event):
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
    #       Peak Parameters                                   
    ###########################################################
    ## Note we also need the parameters in PkParams list to
    ## get updated on an enter event in any of the edit fields
    ## or mouseclick on sliders ...
    ## at the moment need to click add to update...
    ## see the event thing for energy fields

    def init_PkParams(self):
        #self.components.DetSelect.items = ['0']
        #self.components.DetSelect.stringSelection = '0'
        self.components.PkParams.items = []
        self.components.PkLbl.text                = ''
        self.components.PkEn.text                 = '0.0'
        self.components.PkAmp.text                = '0.0'
        self.components.PkFWHM.text               = '0.0'
        self.components.PkEnFlag.text             = '0'
        self.components.PkFwhmFlag.text           = '0'
        self.components.PkAmpFactor.text          = '0.0'
        self.components.PkIgnore.checked          = False
        return
    
    def on_PkParams_select(self,event):
        #print self.components.PkParams.items
        selected =  self.components.PkParams.getStringSelection()
        pk_params = selected[0]
        self.put_PkPar_fields(pk_params)
        return

    def get_PeakParams(self):
        """ return list of all Peak Parameters from GUI """
        return self.components.PkParams.items
    
    def PkParams_update(self,pk_params):        
        "update the pk params list"
        list = self.components.PkParams.items

        if len(list) == 0:
            list.append(pk_params)
            sel = 0
        else:
            found = False
            for j in range(len(list)):
                if list[j][0] == pk_params[0] and list[j][1] == pk_params[1]:
                    list[j] = pk_params
                    found = True
                    sel = j
                    break
            if found == False:
                list.append(pk_params)
                sel = last = len(list)-1
        self.components.PkParams.items = list
        
        # this should invoke the on_PkParams_select event
        #self.components.PkParams.SetSelection(sel)

        return

    def get_PkPar_fields(self):
        "ret a list of all the peak parameter info in the entry fields"
        pk_params = ['']*10
        # gather all the info
        det = self.components.DetSelect.stringSelection
        lbl = self.components.PkLbl.text.strip()
        if len(lbl) < 1: return None
        pk_params[0] = det
        pk_params[1] = lbl
        pk_params[2] = self.components.PkEn.text.strip()
        pk_params[3] = self.components.PkAmp.text.strip()  
        pk_params[4] = self.components.PkFWHM.text.strip() 
        pk_params[5] = self.components.PkEnFlag.text.strip() 
        pk_params[6] = self.components.PkFwhmFlag.text.strip() 
        pk_params[7] = self.components.PkAmpFactor.text.strip()
        pk_params[8] = str(self.components.PkIgnore.checked)
        pk_params[9] = '0.0'
        return pk_params

    def put_PkPar_fields(self,pk_params):
        "reverse of above"
        if pk_params == None: return
        if len(pk_params) < 9: return
        self.components.DetSelect.stringSelection = pk_params[0]
        self.components.PkLbl.text                = pk_params[1]
        self.components.PkEn.text                 = pk_params[2]
        self.components.PkAmp.text                = pk_params[3]
        self.components.PkFWHM.text               = pk_params[4]
        self.components.PkEnFlag.text             = pk_params[5]
        self.components.PkFwhmFlag.text           = pk_params[6]
        self.components.PkAmpFactor.text          = pk_params[7]
        self.components.PkIgnore.checked          = eval(pk_params[8])
        return
    
    ###########################################################
    #       Bgr Parameters                                    
    ###########################################################

    def init_BgrParams(self):
        self.components.BgrParams.items = []
        self.components.BgrCheck.checked == True
        self.components.BgrExp.text = '2'
        self.components.BgrTopWdth.text = '0.0'
        self.components.BgrBtmWdth.text = '4.0'
        self.components.BgrTangent.text = '0'
        self.components.BgrCompress.text = '4'
        return
    
    def on_BgrCheck_mouseClick(self,event):
        if self.components.BgrCheck.checked == True:
            self.components.BgrExp.editable = True
            self.components.BgrExp.backgroundColor = (255, 255, 255)
            self.components.BgrTopWdth.editable = True
            self.components.BgrTopWdth.backgroundColor = (255, 255, 255)
            self.components.BgrBtmWdth.editable = True
            self.components.BgrBtmWdth.backgroundColor = (255, 255, 255)
            self.components.BgrTangent.editable = True
            self.components.BgrTangent.backgroundColor = (255, 255, 255)
            self.components.BgrCompress.editable = True
            self.components.BgrCompress.backgroundColor = (255, 255, 255)
            self.components.BgrDefault.enabled = True
        else:
            self.components.BgrExp.editable = False
            self.components.BgrExp.backgroundColor = (192, 192, 192)
            self.components.BgrTopWdth.editable = False
            self.components.BgrTopWdth.backgroundColor = (192, 192, 192)
            self.components.BgrBtmWdth.editable = False
            self.components.BgrBtmWdth.backgroundColor = (192, 192, 192)
            self.components.BgrTangent.editable = False
            self.components.BgrTangent.backgroundColor = (192, 192, 192)
            self.components.BgrCompress.editable = False
            self.components.BgrCompress.backgroundColor = (192, 192, 192)
            self.components.BgrDefault.enabled = False
        return

    def on_BgrDefault_mouseClick(self,event):
        self.components.BgrExp.text = '2'
        self.components.BgrTopWdth.text = '0.0'
        self.components.BgrBtmWdth.text = '4.0'
        self.components.BgrTangent.text = '0'
        self.components.BgrCompress.text = '4'
        return

    def on_AddBgr_mouseClick(self,event):
        "add bgr params"
        # get bgr params from fields
        bgr_params = self.get_BgrPar_fields()
        if bgr_params == None: return
        # update bgr params list
        self.BgrParams_update(bgr_params)
        return

    def on_DelBgr_mouseClick(self,event):
        selected =  self.components.BgrParams.getStringSelection()
        if selected == None: return
        list = self.components.BgrParams.items
        for sel in selected:
            j = 0
            while j < len(list):
                if sel[0] == list[j][0] and sel[1] == list[j][1]:
                    list.pop(j)
                    break
                j = j + 1
        self.components.BgrParams.items = list

        # this should invoke the on_BgrParams_select event
        self.components.BgrParams.SetSelection(0)

        return

    def on_BgrParams_select(self,event):
        #print self.components.PkParams.items
        selected =  self.components.BgrParams.getStringSelection()
        bgr_params = selected[0]
        self.put_BgrPar_fields(bgr_params)
        return

    def get_BgrParams(self):
        """ return list of all Bgr parameters """
        return self.components.BgrParams.items

    def BgrParams_update(self,bgr_params):        
        "given list of bgr_pars update list"
        list = self.components.BgrParams.items

        if len(list) == 0:
            list.append(bgr_params)
            sel = 0
        else:
            found = False
            for j in range(len(list)):
                if list[j][0] == bgr_params[0] :
                    list[j] = bgr_params
                    found = True
                    sel = j
                    break
            if found == False:
                list.append(bgr_params)
                sel = last = len(list)-1
        self.components.BgrParams.items = list
        
        # this should invoke the on_PkParams_select event
        #self.components.BgrParams.SetSelection(sel)
        return

    def get_BgrPar_fields(self):
        "ret a list of all the peak parameter info in the entry fields"
        #bgr_params = BgrParams.copy()
        bgr_params = ['']*6
        bgr_params[0] = self.components.DetSelect.stringSelection
        bgr_params[1] = self.components.BgrExp.text.strip() 
        bgr_params[2] = self.components.BgrTopWdth.text.strip()
        bgr_params[3] = self.components.BgrBtmWdth.text.strip()
        bgr_params[4] = self.components.BgrTangent.text.strip()
        bgr_params[5] = self.components.BgrCompress.text.strip()
        #bgr_params[6] = str(self.components.BgrCheck.checked)  
        return bgr_params
    
    def put_BgrPar_fields(self,bgr_params):
        "reverse above"
        if bgr_params == None: return
        if len(bgr_params) < 6: return

        self.components.DetSelect.stringSelection = bgr_params[0] 
        self.components.BgrExp.text       = bgr_params[1]  
        self.components.BgrTopWdth.text   = bgr_params[2] 
        self.components.BgrBtmWdth.text   = bgr_params[3] 
        self.components.BgrTangent.text   = bgr_params[4] 
        self.components.BgrCompress.text  = bgr_params[5] 
        return

    ###########################################################
    # Copy/Del/Update Parameters Buttons
    ###########################################################
    def on_CopyParams_mouseClick(self,event):
        det      = self.components.DetSelect.stringSelection
        det_lst  = self.components.DetSelect.items

        # copy Bgr
        list = self.get_BgrParams()
        if len(list) > 0:
            # find the bgr params
            found = False
            for j in range(len(list)):
                if list[j][0] == det:
                    bgr_params = list[j]
                    found = True
                    break
            if found == True:
                for d in det_lst:
                    bgr_params[0] = d
                    self.BgrParams_update(bgr_params)
        
        # copy peak
        list = self.get_PeakParams()
        if len(list) > 0:
            # find the peaks params
            pk_params = []
            found = False
            for j in range(len(list)):
                if list[j][0] == det:
                    pk_params.append(list[j])
                    found = True
            if found == True:
                for d in det_lst:
                    for pk_par in pk_params:
                        pk_par[0] = d
                        self.PkParams_update(pk_par)
        return
    
    def on_ClearParamsDet_mouseClick(self,event):
        det = self.components.DetSelect.stringSelection

        # Bgr
        list = self.get_BgrParams()
        if len(list) > 0:
            # find the bgr params
            for j in range(len(list)):
                if list[j][0] == det:
                    list.pop(j)
                    break
            self.components.BgrParams.items = list

        # peak
        list = self.get_PeakParams()
        if len(list) > 0:
            # find the pk params
            while True:
                n = len(list)
                for j in range(n):
                    if list[j][0] == det:
                        list.pop(j)
                        break
                if j >= n-1: break
            self.components.PkParams.items = list
        
        return 
    
    def on_ClearParamsAll_mouseClick(self,event):
        self.init_BgrParams()
        self.init_PkParams()
        pass


    ###########################################################
    # XRF Object
    ###########################################################

    def on_UpdateXRF_mouseClick(self,event):
        "update the xrf data variable given input data"
        self.update_xrf_from_gui()
        self.update_gui_from_xrf()
        return

    def compare_data_params(self):
        """ return True is the xrf data params are the same as
        the parameters listed in the GUI, ie does the xrf object
        need an update"""
        xrf = self.get_xrf()
        if xrf == None: return False

        data_par = self.get_DataPar()
        if eval(data_par['total'])   != xrf.total: return False 
        if eval(data_par['align'])   != xrf.align: return False
        if eval(data_par['correct']) != xrf.correct: return False
        #if data_par['emin']          != xrf.emin: return False
        #if data_par['emax']          != xrf.emax: return False

        bad_mca = self.str_to_list_var(data_par['bad_mcas'],conv=int)
        if bad_mca != xrf.bad_mca_idx: return False

        # note should fix xrf to make this simpler, ie no None's        
        mca_taus = self.str_to_list_var(data_par['mca_taus'],conv=float)
        if mca_taus != xrf.med.tau:
            if mca_taus == [] and xrf.med.tau == None:
                pass
            else:
                return False

        return True

    def compare_fit_params(self):
        """ return True is the xrf fit params are the same as
        the parameters listed in the GUI"""
        return False

    def update_xrf_from_gui(self):
        " update the xrf object from data in GUI"
        xrf = self.get_xrf()
        if xrf == None: return

        # data parameters
        if self.compare_data_params() == False:
            self.post_message("reset data parameters")
            data_par  = self.get_DataPar()
            bad_mca   = self.str_to_list_var(data_par['bad_mcas'],conv=int)
            tau       = self.str_to_list_var(data_par['mca_taus'],conv=float)
            total     = eval(data_par['total'])
            align     = eval(data_par['align'])
            correct   = eval(data_par['correct'])
            #emin      = float(data_par['emin'])
            #emax      = float(data_par['emax'])
            if tau == []:
                #print 'tau is none'
                tau = None
            xrf.init_data(bad_mca_idx=bad_mca,total=total,align=align,
                          correct=correct,tau=tau,init_params=False)

        # fit parameters
        if self.compare_fit_params() == False:
            xrf.init_params()
            
            # bgr parameters
            bgr_params  = self.components.BgrParams.items
            for bgr in bgr_params:
                DetSelect    = bgr[0]
                det_idx      = int(DetSelect)
                BgrExp       = int(bgr[1])  
                BgrTopWdth   = float(bgr[2]) 
                BgrBtmWdth   = float(bgr[3]) 
                BgrTangent   = float(bgr[4]) 
                BgrCompress  = int(bgr[5]) 
                xrf.set_bgr(slope=None,exponent=BgrExp,top_width=BgrTopWdth,
                            bottom_width=BgrBtmWdth, tangent=BgrTangent,
                            compress=BgrCompress,det_idx=det_idx)

            # peak parameters
            peak_params = self.components.PkParams.items
            for pk in peak_params:
                #print pk
                DetSelect   = pk[0]
                det_idx     = int(DetSelect)
                PkLbl       = pk[1]
                PkEn        = float(pk[2])
                PkAmp       = float(pk[3])
                PkFWHM      = float(pk[4])
                PkEnFlag    = int(pk[5])
                PkFwhmFlag  = int(pk[6])
                PkAmpFactor = float(pk[7])
                PkIgnore    = eval(pk[8])
                pk_idx = xrf.init_peak_en(label=PkLbl,energy=PkEn,det_idx=det_idx)

                xrf.set_peak(pk_idx=pk_idx,energy=PkEn,ampl=PkAmp,fwhm=PkFWHM,
                             energy_flag=PkEnFlag,fwhm_flag=PkFwhmFlag,
                             ampl_factor=PkAmpFactor,ignore=PkIgnore,det_idx=det_idx)

        # post updated xrf
        self.set_xrf(xrf)
        return

    def update_gui_from_xrf(self):
        xrf = self.get_xrf()
        if xrf == None:
            self.post_message("No XRF variable")
            return 
        try:
            #data parameters
            num_mca = xrf.med.n_detectors
            self.components.NumMcas.text = "NumMcas = %i" % num_mca

            self.components.Total.checked = xrf.total
            self.components.Align.checked = xrf.align
            self.components.CorrectData.checked = xrf.correct
            #emin      = str(xrf.emin)
            #emax      = str(xrf.emax)
            self.components.BadMcas.text = str(xrf.bad_mca_idx)
            self.components.McaTaus.text = str(xrf.med.tau)

            # update det select for fitting
            # should this be a seperate function
            # should the DetSelect list skip 'bads'?
            det_idx_str = []
            for idx in range(xrf.ndet):
                det_idx_str.append(str(idx))
            det_sel = self.components.DetSelect.stringSelection
            self.components.DetSelect.items = det_idx_str
            if det_sel in det_idx_str:
                self.components.DetSelect.stringSelection = det_sel
            else:
                self.components.DetSelect.stringSelection = '0'
            
            # get parameters 
            (fit,bgr,pk) = xrf.get_params()

            # BgrParams
            self.init_BgrParams()
            try:
                for j in range(len(bgr)):
                    bgr_params = ['']*6
                    bgr_params[0] = str(j)
                    bgr_params[1] = str(bgr[j].exponent)
                    bgr_params[2] = str(bgr[j].top_width)
                    bgr_params[3] = str(bgr[j].bottom_width)
                    bgr_params[4] = str(bgr[j].tangent)
                    bgr_params[5] = str(bgr[j].compress)
                    self.BgrParams_update(bgr_params)
            except:
                pass

            #PeakParams
            self.init_PkParams()
            try:
                for j in range(len(pk)):
                    for k in range(len(pk[j])):
                        pk_params = ['']*10
                        pk_params[0] = str(j)
                        pk_params[1] = pk[j][k].label
                        pk_params[2] = str(pk[j][k].energy)
                        pk_params[3] = str(pk[j][k].ampl)
                        pk_params[4] = str(pk[j][k].fwhm)
                        pk_params[5] = str(pk[j][k].energy_flag)
                        pk_params[6] = str(pk[j][k].fwhm_flag)
                        pk_params[7] = str(pk[j][k].ampl_factor)
                        pk_params[8] = str(pk[j][k].ignore)
                        pk_params[9] = str(pk[j][k].area)
                        self.PkParams_update(pk_params)
            except:
                pass

            #done
            return    
        except:
            self.post_message("Failed to read XRF object ")
            return

    ###########################################################
    # Save/Restore
    ###########################################################

    def init_SaveParVarItems(self):
        " Initialize the items from tdl    "
        lst = self.listAllData()
        t = self.components.SaveParVar.text
        self.components.SaveParVar.items = lst
        self.components.SaveParVar.text = t
        return

    def on_SaveSavePar_mouseClick(self,event):
        self.save_parameters()
        self.init_tdl_list_items()
        return

    def on_RestoreSavePar_mouseClick(self,event):
        self.update_gui_from_save()
        return    

    def save_parameters(self):
        var_name = self.components.SaveParVar.text
        if len(var_name.strip()) == 0: return
        data_par = self.get_DataPar()
        bgr_par = self.get_BgrParams()
        pk_par  = self.get_PeakParams()
        fit_par = self.get_FitArgs_fields()

        vars = {'data_par':data_par,'bgr_par':bgr_par,
                'pk_par':pk_par,'fit_par':fit_par}
        self.setValue(var_name,vars)
        return

    def update_gui_from_save(self):
        try:
            xrf = self.get_xrf() 
            num_mca = xrf.med.n_detectors
            self.components.NumMcas.text = "NumMcas = %i" % num_mca
        except:
            self.components.NumMcas.text = "NumMcas = 0"

        # get save parameters
        var_name = self.components.SaveParVar.text
        if len(var_name.strip()) == 0: return
        save_var = self.getValue(var_name)
        try:
            data_par = save_var['data_par']
            bgr_par = save_var['bgr_par']
            pk_par = save_var['pk_par']
            fit_par = save_var['fit_par']
            
            # put stuff
            self.put_DataPar(data_par)
            self.components.BgrParams.items = bgr_par
            self.components.PkParams.items  = pk_par
            self.put_FitArgs_fields(fit_par)
        except:
            self.post_message("failed to read save parameters")
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
        xrf.calc_peaks()
        self.set_xrf(xrf)
        self.update_gui_from_xrf()
        if self.components.AutoUpdateCheck.checked == True:
            self.plot_cmd()
        return

    def init_FitParams(self):
        self.components.FitFWHMFlag.stringSelection = '1' 
        self.components.FitEnFlag.stringSelection = '1'
        self.components.FitChiExp.stringSelection = '0.0'
        self.components.BgrCheck.checked = True 

    def get_FitArgs_fields(self):
        fit_args = ['']*4
        fit_args[0] = self.components.FitFWHMFlag.stringSelection
        fit_args[1] = self.components.FitEnFlag.stringSelection 
        fit_args[2] = self.components.FitChiExp.stringSelection
        fit_args[3] = str(self.components.BgrCheck.checked) 
        #print fit_args
        return fit_args

    def put_FitArgs_fields(self,fit_args):
        self.components.FitFWHMFlag.stringSelection = fit_args[0] 
        self.components.FitEnFlag.stringSelection = fit_args[1]
        self.components.FitChiExp.stringSelection = fit_args[2]
        self.components.BgrCheck.checked = eval(fit_args[3]) 
        return

    def on_Fit_mouseClick(self,event):
        self.fit()
        return

    def on_FitScan_mouseClick(self,event):
        self.fit_scan()
        return

    def fit(self):
        self.update_xrf_from_gui()
        xrf         = self.get_xrf()
        fit_args    = self.get_FitArgs_fields()
        FitFWHMFlag = int(fit_args[0]) 
        FitEnFlag   = int(fit_args[1])  
        FitChiExp   = float(fit_args[2])
        fit_bgr     = eval(fit_args[3])
        if self.components.Total.checked == True:
            xrf.fit(fwhm_flag=FitFWHMFlag,energy_flag=FitEnFlag,
                    chi_exp=FitChiExp,fit_bgr=fit_bgr)
        else:
            det_idx = int(self.components.DetSelect.stringSelection)
            if fit_bgr: xrf._fit_bgr(det_idx)
            xrf._fit_peaks(det_idx,fwhm_flag=FitFWHMFlag,
                           energy_flag=FitEnFlag,chi_exp=FitChiExp)
        self.set_xrf(xrf)
        self.update_gui_from_xrf()
        if self.components.AutoUpdateCheck.checked == True:
            self.plot_cmd()
        return

    def fit_scan(self):
        # must be a scan
        if not self.components.ScanCheck.checked:
            return
        # must have save params set 
        if not self.components.SetParFromSave.checked:
            self.components.SetParFromSave.checked = True
        if len(self.components.SaveParVar.text.strip()) == 0:
            self.post_message("You must chose save parameters")
            return 
        # shut off autoupdate
        if self.components.AutoUpdateCheck.checked:
            reset_auto = True
            self.components.AutoUpdateCheck.checked = False
        else:
            reset_auto = False
            
        # fit all the spectra
        st = int(self.components.ScanStart.text)
        en = int(self.components.ScanEnd.text)
        for j in range(st,en+1):
            node_name = self.BuildScanName(j)
            self.components.NodePfx.text = node_name
            # not sure need below?
            self.update_xrf_from_gui()
            xrf = self.get_xrf()
            fit_args = self.get_FitArgs_fields()
            FitFWHMFlag = int(fit_args[0]) 
            FitEnFlag   = int(fit_args[1])  
            FitChiExp   = float(fit_args[2]) 
            fit_bgr     = eval(fit_args[3])
            xrf.fit()
            self.set_xrf(xrf)
            
        # tun auto update back on if needed
        if reset_auto:
            self.components.AutoUpdateCheck.checked = True
            
        # shut off set parameters from save
        self.components.SetParFromSave.checked = False
        
        return

    ###########################################################
    # Extract Results
    ###########################################################
    def init_ExtractResultsVarItems(self):
        " Initialize the items from tdl    "
        lst = self.listAllData()
        t = self.components.ExtractResultsVar.text
        self.components.ExtractResultsVar.items = lst
        self.components.ExtractResultsVar.text = t
        return

    def on_ExtractResults_mouseClick(self,event):
        self.extract_results()
        self.init_tdl_list_items()
        return

    def extract_results(self):
        res_name = self.components.ExtractResultsVar.text
        if len(res_name.strip()) == 0: return
        results = []
        if self.components.ScanCheck.checked == False:
            xrf = self.get_xrf()
            results.append(xrf.get_peaks())
            self.setValue(res_name,results)
        else:
            st = int(self.components.ScanStart.text)
            en = int(self.components.ScanEnd.text)
            for j in range(st,en+1):
                node_name = self.BuildScanName(j)
                var_name = self.get_xrf_var_name(node=node_name)
                self.post_message(var_name)
                xrf = self.getValue(var_name)
                print xrf
                results.append(xrf.get_peaks())
            self.setValue(res_name,results)
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
        if not self.check_xrf_var():
            return

        # build cmd
        xrf   = self.get_xrf_var_name()
        if self.components.FitPlotCheck.checked:
            if self.components.Total.checked == True:
                cmd_str = "xrf.plot_fit(%s" % xrf
            else:
                det_idx = int(self.components.DetSelect.stringSelection)
                cmd_str = "xrf.plot_fit_det(%s, det=%d" % (xrf,det_idx)
        else:
            if self.components.Total.checked == True:
                cmd_str = "xrf.plot(%s" % xrf
            else:
                det_idx = int(self.components.DetSelect.stringSelection)
                cmd_str = "xrf.plot_det(%s, det=%d" % (xrf, det_idx)

        if self.components.YlogCheck.checked:
            cmd_str = cmd_str + ',logplt=True)'
        else:
            cmd_str = cmd_str + ')'

        self.post_message(cmd_str)
        self.execLine(cmd_str)
        #self.update_gui_from_xrf()
        return

    def plot_par_update(self):
        """Update self.plot_par from components"""
        pass
    
    def plot_par_display(self):
        """Update components from self.plot_par"""
        pass


##################################################
if __name__ == '__main__':
    app = model.Application(wxXRF)
    app.MainLoop()
