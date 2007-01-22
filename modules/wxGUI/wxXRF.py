############################################################################
# Tom Trainor (fftpt@uaf.edu)
# This should be run as a child window from wxGUI
#
# --------------
# Modifications:
# --------------
#
# --------------
# Todo
# --------------
# - work on plotting....(need tdl to plot!)
# - note need some error checking so set_data and plot commands dont
#   get run if there is an invalid variable name
#   - e.g. have data_par_update return True/False if var is defined ??
#
############################################################################

from PythonCard import model, dialog
import wx
import xrf_lookup

"""
XRF GUI
"""

DataPar= {'grp':'None','node':'None','total':'True','align':'True',
          'scan':'False','sc_start':'0','sc_stop':'0','sc_inc':'1','sc_step':'0',
          'mcas':'[]','mca_taus':'[]','emin':'0.0','emax':'0.0'}

PlotPar= {'plt_data':True,'plt_fit':False,'plt_auto_update':False,
          'plt_hold':False,'plt_components':False,
          'plt_xlog':False,'plt_ylog':False,'plt_yerr':False}
#
#BgrParams=[[DetSelect,BgrExp,BgrTopWdth,BgrBtmWdth,
#           BgrTangent,BgrCompress,BgrCheck]]
#
#PeakParams= [[Det,Label,PkEn,PkAmp,PkFWHM,PkEnFlag,PkFwhmFlag,
#             PkAmpFactor,PkIgnore,PkArea]]
#
#FitArgs= [FitFWHMFlag,FitEnFlag,FitChiExp,InitScanIdx,UsePrevFit]
#

class wxXRF(model.Background):

    ###########################################################
    # Init and util methods
    ###########################################################
    
    def on_initialize(self, event):
        # if you have any initialization
        # including sizer setup, do it here
        # are there potential race conditions?
        # window management???

        # including sizer setup, do it here
        #self.setupSizers()

        self.shell = None
        self.tdl = None
        self.init_tdl()

        # Flags etc
        self.NewData = True
        self.fit_det_select = '0'

        # parameters and data
        self.data_par = DataPar
        self.plot_par = PlotPar
        self.bgr_pars=[]
        self.peak_pars=[]
        self.fit_args = []
        self.components.PkParams._autoresize = 0
        self.components.BgrParams._autoresize = 0

    def init_tdl(self,shell=None):
        # self.shell = shell
        # get the tdl reference from the parent window 
        self.shell = self.getParent().get_shell()
        self.tdl = self.shell.tdl
        self.init_menus()

    def init_menus(self):   
        self.init_GrpItems()
        self.init_NodePfxItems()
        self.init_XrfLineItems()
        self.init_DetAndTauItems()
        self.components.DetSelect.items = ['All']

    def check_xrf_var(self):
        node = self.components.NodePfx.text
        grp = self.components.Grp.text  
        if len(grp) == 0:
            var_name = node
        else:
            var_name = "%s.%s" % (grp,node)
        m = self.getValue(var_name)
        try:
            num_mca = m.med.n_detectors
            self.components.NumMcas.text = "NumMcas = %i" % num_mca
            det = str(m.bad_mca_idx)
            self.components.McaList.text = det
            tau = str(m.med.tau)
            self.components.TauList.text = tau
            # update det select for fitting
            det_idx_str = ['All']
            for idx in range(m.ndet):
                det_idx_str.append(str(idx))
            det_sel = self.components.DetSelect.stringSelection
            self.components.DetSelect.items = det_idx_str
            for det_sel in det_idx_str:
                self.components.DetSelect.stringSelection = det_sel

            # get parameters from the xrf object and update in the gui
            self.PkParams_init()
            self.BgrParams_init()
            self.get_xrf_fit_params()
            
        except:
            self.components.NumMcas.text = "Variable is not an XRF data type"
            self.components.McaList.text = '[]'
            self.components.TauList.text = 'None'
            self.components.DetSelect.items = ['All']
        return

    ###########################################################
    # Tdl Utilities
    ###########################################################

    def setValue(self,var_name,value):
        self.tdl.setVariable(var_name,value)

    def getValue(self,var_name):
        return self.tdl.getVariableValue(var_name)

    def getVariable(self,var_name):
        return self.tdl.getVariable(var_name)

    def listGroups(self):
        return self.tdl.symbolTable.listGroups()

    def listDataGroup(self,grp):
        return self.tdl.symbolTable.listDataGroup(grp)

    def listAllData(self):
        data_dict = self.tdl.symbolTable.listData()
        lst = []
        for grp in data_dict.keys():
            for node in data_dict[grp]:
                name = "%s.%s" % (grp,node)
                lst.append(name)
        lst.sort()
        return lst

    def execLine(self,line):
        #return self.tdl.execute(line)
        return self.tdl.eval(str(line))

    def str_to_list(self,s,conv=float):
        if s == None: return []
        s = s.strip()
        if len(s) == 0: return []
        if s[0] == '=': s = s[1:]
        if s[0] == '[':
            s = s[1:]
            if s[-1] == ']':
                s = s[:-1]
        lst = s.split(',')
        return map(conv, lst)

    def post_message(self,mess):
        """ write a message to the screen"""
        txt = "\n---\n%s\n" % mess 
        #self.components.Messages.text = txt
        self.components.Messages.appendText(txt)
        return

    ###########################################################
    #             Menus                                       #
    ###########################################################

    def on_menuFileExit_select(self,event):        
        self.close()
        
    def on_menuFileReadXrf_select(self,event):
        #self.post_message('hello')
        dir = '.'
        result = dialog.fileDialog(self, 'Open...', dir, '',"*")
        if result.accepted:
            path = result.paths[0]
        else:
            self.status("File insert cancelled.")
            return 0
        for path in result.paths:
            line = "xrf.read '%s'" % path
            self.execLine(line)
        
        # update the variables lists
        self.init_GrpItems()
        # set select list to default grp etc..
        self.components.Grp.text = 'xrf_data'
        self.components.NodePfx.text = ''
        self.init_NodePfxItems()
        
        return


    ###########################################################
    #             Components and Events                       #
    ###########################################################

    #--------------
    # Update - button
    #--------------
    def on_UpdateMenus_mouseClick(self, event):
        "use this to force update incase the data list has changed"
        self.init_menus()
        return

    #--------------
    # GrpItems - choice
    #--------------
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

    #------------------
    # NodePrefix - choice
    #------------------
    def init_NodePfxItems(self):
        """Initialize the menu. Use the group thats selected in the group menu"""
        grp = self.components.Grp.text  
        if len(grp) == 0: grp = None
        node = self.components.NodePfx.text
        self.components.NodePfx.items = self.listDataGroup(grp)
        self.components.NodePfx.text = node
        return

    def on_NodePfx_select(self,event):
        "select a variable name and check it out"
        self.check_xrf_var()
        return

    def on_NodePfx_textUpdate(self,event):
        "select a variable name and check it out"
        self.check_xrf_var()
        return

    #------------------
    # Choice box for det and tau variables
    #------------------
    def init_DetAndTauItems(self):
        " Initialize the menu    "
        lst = self.listAllData()
        t1 = self.components.McaVariable.text
        t2 = self.components.TauVariable.text
        self.components.McaVariable.items = lst
        self.components.TauVariable.items = lst
        self.components.McaVariable.text = t1
        self.components.TauVariable.text = t2
        return
    
    #------------------
    # UpdateData - button
    #------------------
    def on_UpdateData_mouseClick(self,event):
        "update the xrf data variable given input data"
        self.set_xrf_data()
        return 

    #------------------
    # Energy Range - button
    #------------------
    def on_EnRange_mouseClick(self,event):
        """ Select energy Range from plot"""
        self.plot_cmd()
        # select Emin/Emax
        self.set_xrf_data()
        return

    #------------------
    # Plot - button 
    #------------------
    def on_Plot_mouseClick(self,event):
        "make a plot"
        self.plot_cmd()
        return

    #------------------
    # XrfLine - choice
    #------------------
    def init_XrfLineItems(self):
        " Initialize the menu    "
        if self.components.Emin:
            Emin = float(self.components.Emin.text)
        else:
            Emin = 0.01
        if self.components.Emax:
            Emax = float(self.components.Emax.text)
        else:
            Emax = 40.0
        self.components.XrfLine.items = xrf_lookup.list_lines(Emin,Emax)

    def on_XrfLine_select(self, event):
        "select an xrf line"
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

    #------------------
    # XrfLineAdd - button
    #------------------
    def on_XrfLineAdd_mouseClick(self,event):
        "add an xrf line"
        # get peak params from fields
        pk_params = self.get_PkPar_fields()
        if pk_params == None: return
        # update the list
        self.PkParams_update(pk_params)
        return

    #------------------
    # XrfLineDel - button
    #------------------
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

    #------------------
    # Peak params box
    #------------------
    def PkParams_init(self):
        self.components.PkParams.items = []
        return
    
    def on_PkParams_select(self,event):
        #print self.components.PkParams.items
        selected =  self.components.PkParams.getStringSelection()
        pk_params = selected[0]
        self.put_PkPar_fields(pk_params)
        return

    ## Note we also need the parameters in PkParams list to
    ## get updated on an enter event in any of the edit fields
    ## or mouseclick on sliders / ignore box...
    ## at the moment need to click add to update...

    #------------------
    # Bgr params box
    #------------------
    def BgrParams_init(self):
        self.components.BgrParams.items = []
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

    #------------------
    # BgrParamsAdd - button
    #------------------
    def on_AddBgr_mouseClick(self,event):
        "add bgr params"
        # get bgr params from fields
        bgr_params = self.get_BgrPar_fields()
        if bgr_params == None: return
        # update bgr params list
        self.BgrParams_update(bgr_params)
        return

    #------------------
    # BgrDel - button
    #------------------
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

    #------------------
    # Bgr params select
    #------------------
    def on_BgrParams_select(self,event):
        #print self.components.PkParams.items
        selected =  self.components.BgrParams.getStringSelection()
        bgr_params = selected[0]
        self.put_BgrPar_fields(bgr_params)
        return

    #------------------
    # Update params button
    #------------------
    def on_UpdateParams_mouseClick(self,event):
        self.set_xrf_fit_params()

    def on_Calc_mouseClick(self,event):
        self.set_xrf_fit_params()
        self.calc()
        return

    def on_Fit_mouseClick(self,event):
        self.set_xrf_fit_params()
        self.fit()
        return

    def on_FitScan_mouseClick(self,event):
        pass


    ###########################################################
    # Data Handling etc.
    ###########################################################
    #------------------
    # Data pars
    #------------------ 
    def data_par_update(self):
        """ Updata self.data_par from the info in GUI components
        Note all the data is stored as strings"""

        # grp
        grp = self.components.Grp.text  
        if len(grp) == 0: grp = None
        self.data_par['grp'] = grp

        # node
        node = self.components.NodePfx.text
        if len(node) == 0: node = None
        self.data_par['node'] = node

        #total
        total = self.components.Total.checked
        self.data_par['total'] = str(total)

        #align
        align = self.components.Align.checked
        self.data_par['align'] = str(align)

        #scan
        scan = self.components.AutoIncCheck.checked
        self.data_par['scan'] = str(scan)

        #sc_start
        sc_start = self.components.AutoIncStart.text
        self.data_par['sc_start'] = str(int(sc_start))

        #sc_stop
        sc_stop = self.components.AutoIncStop.text
        self.data_par['sc_stop'] = str(int(sc_stop))

        #sc_inc
        sc_inc = self.components.AutoIncInc.text
        self.data_par['sc_inc'] = str(int(sc_inc))

        #sc_step
        sc_step = self.components.NodeInc.value  #??????
        self.data_par['sc_step'] = str(int(sc_step))

        #mcas
        # --> if variable is not '' get vals from tdl
        mcas = self.components.McaList.text
        self.data_par['mcas'] = mcas

        #mca taus
        # --> if variable is not '' get vals from tdl
        mca_taus = self.components.TauList.text
        self.data_par['mca_taus'] = mca_taus

        # Emin/Emax
        Emin = self.components.Emin.text
        Emax = self.components.Emax.text
        OldEmin = self.data_par['emin']
        OldEmax = self.data_par['emax']
        self.data_par['emin'] = Emin
        self.data_par['emax'] = Emax
        
        # this should ret true if E-range changed, ie reset
        # the xrf lines in range
        if (OldEmin != Emin) or (OldEmax != Emax):
            return True
        else:
            return False

    def data_par_display(self):
        """ Update gui components from self.data_par"""
        # update from the xrf obj 
        pass

    #------------------
    # Peak pars
    #------------------
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
        self.components.PkParams.SetSelection(sel)

        return

    #------------------
    # Peak par fields
    #------------------
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
        self.components.PkIgnore.checked = eval(pk_params[8])
        #if pk_params[8] == 'False':
        #    self.components.PkIgnore.checked = False
        #elif pk_params[8] == 'True':
        #    self.components.PkIgnore.checked = True
        return
    
    #------------------
    # Bgr pars
    #------------------
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
        self.components.BgrParams.SetSelection(sel)

        return

    #------------------
    # Bgr par fields
    #------------------
    def get_BgrPar_fields(self):
        "ret a list of all the peak parameter info in the entry fields"
        #bgr_params = BgrParams.copy()
        bgr_params = ['']*7
        bgr_params[0] = self.components.DetSelect.stringSelection
        bgr_params[1] = self.components.BgrExp.text.strip() 
        bgr_params[2] = self.components.BgrTopWdth.text.strip()
        bgr_params[3] = self.components.BgrBtmWdth.text.strip()
        bgr_params[4] = self.components.BgrTangent.text.strip()
        bgr_params[5] = self.components.BgrCompress.text.strip()
        bgr_params[6] = str(self.components.BgrCheck.checked)  
        return bgr_params
    
    def put_BgrPar_fields(self,bgr_params):
        "reverse above"
        if bgr_params == None: return
        if len(bgr_params) < 7: return

        self.components.DetSelect.stringSelection = bgr_params[0] 
        self.components.BgrExp.text       = bgr_params[1]  
        self.components.BgrTopWdth.text   = bgr_params[2] 
        self.components.BgrBtmWdth.text   = bgr_params[3] 
        self.components.BgrTangent.text   = bgr_params[4] 
        self.components.BgrCompress.text  = bgr_params[5] 
        self.components.BgrCheck.checked    = eval(bgr_params[6])  
        return

    #------------------
    # calc/fit args
    #------------------
    def get_FitArgs_fields(self):
        fit_args = ['']*5
        fit_args[0] = self.components.FitFWHMFlag.stringSelection
        fit_args[1] = self.components.FitEnFlag.stringSelection 
        fit_args[2] = self.components.FitChiExp.stringSelection
        fit_args[3] = self.components.InitScanIdx.stringSelection
        fit_args[4] = str(self.components.UsePrevFit.checked)
        print fit_args
        return fit_args

    def put_FitArgs_fields(self):
        pass


    ###########################################################
    # XRF Data parameters, wrappers, etc.
    ###########################################################

    #------------------
    # xrf
    #------------------
    def get_xrf(self):
        group = self.data_par['grp']
        node  = self.data_par['node']
        xrf_name   = "%s.%s" % (group,node)
        return self.getValue(xrf_name)

    def set_xrf(self,xrf):
        group = self.data_par['grp']
        node  = self.data_par['node']
        xrf_name   = "%s.%s" % (group,node)
        return self.setValue(xrf_name,xrf)

    #------------------
    # xrf data
    #------------------
    def set_xrf_data(self):
        " run xrf.set_data command given data in self.data_par"
        # change to work on xrf directly
        
        # update the data parameters
        do_lines = self.data_par_update()

        # build args
        xrf   = "%s.%s" % (self.data_par['grp'],self.data_par['node'])
        det   = self.data_par['mcas']        
        total = self.data_par['total']
        align = self.data_par['align']
        tau   = self.data_par['mca_taus']
        if tau == '' or tau == '[]': tau = 'None'

        # build cmd str
        cmd_str = "xrf.set_data(%s,detectors=%s,total=%s,align=%s,tau=%s)" % \
                  (xrf,det,total,align,tau)
        self.post_message(cmd_str)
        self.execLine(cmd_str)

        # updates
        if do_lines: self.init_XrfLineItems()
        self.check_xrf_var()
        
        return

    #------------------
    # xrf fit parameters
    #------------------
    def set_xrf_fit_params(self):
        # update the xrf obj
        bgr_params = self.components.BgrParams.items
        peak_params = self.components.PkParams.items

        # get xrf
        xrf   = self.get_xrf()
        if xrf == None: return
        
        # blows away everything in xrf
        xrf.init_params()

        # update with new stuff
        for pk in peak_params:
            DetSelect   = pk[0]
            PkLbl       = pk[1]
            PkEn        = float(pk[2])
            PkAmp       = float(pk[3])
            PkFWHM      = float(pk[4])
            PkEnFlag    = int(pk[5])
            PkFwhmFlag  = int(pk[6])
            PkAmpFactor = float(pk[7])
            PkIgnore    = eval(pk[8])
            if DetSelect == "All":
                for idx in range(xrf.array_len):
                    pk_idx = xrf.init_peak_en(idx=idx,label=PkLbl)
                    
                    xrf.set_peak(idx=idx,pk_idx=pk_idx,energy=PkEn,ampl=PkAmp,fwhm=PkFWHM,
                                 energy_flag=PkEnFlag,fwhm_flag=PkFwhmFlag,
                                 ampl_factor=PkAmpFactor,ignore=PkIgnore)
            else:
                idx = int(DetSelect)
                pk_idx = xrf.init_peak_en(idx=idx,label=PkLbl)

                xrf.set_peak(idx=idx,pk_idx=pk_idx,energy=PkEn,ampl=PkAmp,fwhm=PkFWHM,
                             energy_flag=PkEnFlag,fwhm_flag=PkFwhmFlag,
                             ampl_factor=PkAmpFactor,ignore=PkIgnore)

        for bgr in bgr_params:
            DetSelect    = bgr[0] 
            BgrExp       = int(bgr[1])  
            BgrTopWdth   = float(bgr[2]) 
            BgrBtmWdth   = float(bgr[3]) 
            BgrTangent   = float(bgr[4]) 
            BgrCompress  = int(bgr[5]) 
            BgrCheck     = eval(bgr[6])  
            if DetSelect == "All":
                for idx in range(xrf.array_len):
                    xrf.set_bgr(idx=idx,slope=None,exponent=BgrExp,top_width=BgrTopWdth,
                                bottom_width=BgrBtmWdth, tangent=BgrTangent, compress=BgrCompress)
            else:
                xrf.set_bgr(idx=int(DetSelect),slope=None,exponent=BgrExp,top_width=BgrTopWdth,
                            bottom_width=BgrBtmWdth, tangent=BgrTangent, compress=BgrCompress)

        # post updated xrf
        self.set_xrf(xrf)

        return

    ##############################################################
    def get_xrf_fit_params(self):
        "read params from xrf and ret fit_params, bgr_params"
        # get xrf
        xrf   = self.get_xrf()
        if xrf == None: return

        (fit,bgr,pk) = xrf.get_params()
        print fit,bgr,pk  # mcafit, mcabgr,mcapeaks

        #PeakParams= [[Det,Label,PkEn,PkAmp,PkFWHM,PkEnFlag,PkFwhmFlag,
        #      PkAmpFactor,PkIgnore,PkArea]]
        print len(pk)
        for j in range(len(pk)):
            print len(pk[j])
            # need fix here if empty...
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

        #BgrParams=[[DetSelect,BgrExp,BgrTopWdth,BgrBtmWdth,
        #           BgrTangent,BgrCompress,BgrCheck]]
        for j in range(len(bgr)):
            bgr_params = ['']*7
            bgr_params[0] = str(j)
            bgr_params[1] = str(bgr[j].exponent)
            bgr_params[2] = str(bgr[j].top_width)
            bgr_params[3] = str(bgr[j].bottom_width)
            bgr_params[4] = str(bgr[j].tangent)
            bgr_params[5] = str(bgr[j].compress)
            bgr_params[6] = str(True)
            self.BgrParams_update(bgr_params)

        return

    #------------------
    # calc/fit
    #------------------
    def calc(self):
        xrf = self.get_xrf()
        xrf.calc_peaks()
        self.set_xrf(xrf)
        return
    
    def fit(self):
        xrf = self.get_xrf()
        fit_args = self.get_FitArgs_fields()
        FitFWHMFlag = int(fit_args[0]) 
        FitEnFlag   = int(fit_args[1])  
        FitChiExp   = float(fit_args[2]) 
        xrf.fit(fwhm_flag=FitFWHMFlag,energy_flag=FitEnFlag,chi_exp=FitChiExp)
        self.set_xrf(xrf)
        self.get_xrf_fit_params()
        return

    def fit_scan(self):
        xrf = self.get_xrf()
        fit_args = self.get_FitArgs_fields()
        FitFWHMFlag = int(fit_args[0]) 
        FitEnFlag   = int(fit_args[1])  
        FitChiExp   = float(fit_args[2]) 
        InitScanIdx = int(fit_args[3]) 
        UsePrevFit  = eval(fit_args[4]) 
        xrf.fit()
        self.set_xrf(xrf)
        return

    #------------------
    # plotting
    #------------------
    def plot_cmd(self):
        "run xrf.plot cmd"
        # update data and plot params
        #self.set_xrf_data()
        #self.plot_par_update()
        # build cmd
        xrf   = "%s.%s" % (self.data_par['grp'],self.data_par['node'])
        if self.components.FitPlotCheck.checked:
            cmd_str = "xrf.plot_fit(%s" % xrf
        else:
            cmd_str = "xrf.plot(%s" % xrf
        if self.components.YlogCheck.checked:
            cmd_str = cmd_str + ',logplt=True)'
        else:
            cmd_str = cmd_str + ')'
        self.post_message(cmd_str)
        self.execLine(cmd_str)
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