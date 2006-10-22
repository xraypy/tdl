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

DataPar= {'grp':None,'node':None,'total':True,'align':True,
          'scan':False,'sc_start':0,'sc_stop':0,'sc_inc':1,'sc_step':0,
          'mcas':[],'mca_taus':[],'emin':0.0,'emax':0.0}

PlotPar= {'plt_data':True,'plt_fit':False,'plt_auto_update':False,
          'plt_hold':False,'plt_components':False,
          'plt_xlog':False,'plt_ylog':False,'plt_yerr':False}

BgrParams={'do_bgr':True,'exponent':2.0,'topwidth':0.0,
           'bottomwidth':4.0,'tangent':0,'compress':4}

PeakParams= {'label':None,'energy':None,'amplitude':0.0,
             'fwhm':0.0,'en_flag':0,'fwhm_flag':0,
             'amp_factor':0.0,'ignore':False,'area':0.0}

FitArgs= {'fit_fwhm_flag':0,'fit_en_flag':0,'chi_exp':0,
          'init_sc_idx':0,'use_prev':False}

#self.fit_par = {'bgr_pars':[],'peak_pars':[],'fit_args':[]}
#self.params={'data_par':{},'plot_par':{},'fit_pars':{}}

class wxXRF(model.Background):

    ###########################################################
    # Init and util methods
    ###########################################################
    
    def on_initialize(self, event):
        # if you have any initialization
        # including sizer setup, do it here
        # are there potential race conditions?
        # window management???
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
        self.fit_args = FitArgs
        self.xrf=None
        self.components.PkParams._autoresize = 0

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
            det = str(m.detectors)
            self.components.McaList.text = det
            tau = str(m.med.tau)
            self.components.TauList.text = tau
            # update det select for fitting
            det_idx_str = ['All']
            for idx in m.detectors: det_idx_str.append(str(idx))
            det_sel = self.components.DetSelect.stringSelection
            self.components.DetSelect.items = det_idx_str
            if det_sel in det_idx_str:
                self.components.DetSelect.stringSelection = det_sel
            else:
                self.components.DetSelect.stringSelection = 'All'
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
        self.post_message('hello')
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

    #def on_NodePfx_mouseDown(self,event):
    #def on_NodePfx_gainFocus(self,event):
    #def on_NodePfx_textUpdate(self,event):
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
        self.set_data_cmd()
        return 

    #------------------
    # Energy Range - button
    #------------------
    def on_EnRange_mouseClick(self):
        """ Select energy Range from plot"""
        self.plot_cmd()
        # select Emin/Emax
        self.set_data_cmd()
        return

    #------------------
    # Plot - button 
    #------------------
    def on_Plot_mouseClick(self,event):
        "make a plot"
        self.set_data_cmd()
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
        self.components.PkLbl.text = line
        self.components.PkEn.text = en

        self.components.PkAmp.text       = '0'
        self.components.PkFWHM.text      = '0'
        self.components.PkEnFlag.text    = '0'
        self.components.PkFwhmFlag.text  = '0'
        self.components.PkAmpFactor.text = '0'
        self.components.PkIgnore.clicked = False
        return

    #------------------
    # XrfLineAdd - button
    #------------------
    def on_XrfLineAdd_mouseClick(self,event):
        "add an xrf line"
        # gather all the info
        det = self.components.DetSelect.stringSelection
        lbl = self.components.PkLbl.text.strip()
        if len(lbl) < 1: return
        en  = self.components.PkEn.text.strip()
        amp = self.components.PkAmp.text.strip()  
        fwhm = self.components.PkFWHM.text.strip() 
        enfl = self.components.PkEnFlag.text.strip() 
        fwfl = self.components.PkFwhmFlag.text.strip() 
        ampfac = self.components.PkAmpFactor.text.strip()
        if self.components.PkIgnore.clicked == False:
            ignore = '0'
        elif self.components.PkIgnore.clicked == True:
            ignore = '1'

        # grab the PkParam list  
        list = self.components.PkParams.items
        print list
        if len(list) == 0:
            list.append([det,lbl,en,amp,fwhm,enfl,fwfl,ampfac,ignore,'0'])
            sel = 0
        else:
            found = False
            for j in range(len(list)):
                if list[j][0] == det and list[j][1] == lbl:
                    list[j] = [det,lbl,en,amp,fwhm,enfl,fwfl,ampfac,ignore,'0']
                    found = True
                    sel = j
                    break
            if found == False:
                list.append([det,lbl,en,amp,fwhm,enfl,fwfl,ampfac,ignore,'0'])
                sel = last = len(list)-1
        self.components.PkParams.items = list

        # this should invoke the on_PkParams_select event
        self.components.PkParams.SetSelection(sel)

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
    # Peak params select
    #------------------
    def on_PkParams_select(self,event):
        #print self.components.PkParams.items
        selected =  self.components.PkParams.getStringSelection()
        ed = selected[0]
        #print '\nselcted', selected
        #print '\ned', selected[0]
        self.components.DetSelect.stringSelection = ed[0]
        self.components.PkLbl.text                = ed[1]
        self.components.PkEn.text                 = ed[2]
        self.components.PkAmp.text                = ed[3]
        self.components.PkFWHM.text               = ed[4]
        self.components.PkEnFlag.text             = ed[5]
        self.components.PkFwhmFlag.text           = ed[6]
        self.components.PkAmpFactor.text          = ed[7]
        if ed[8] == '0':
            self.components.PkIgnore.clicked = False
        elif ed[8] == '1':
            self.components.PkIgnore.clicked = True
        return

    ## Note we also need the parameters in PkParams list to
    ## get updated on an enter event in any of the edit fields
    ## or mouseclick on sliders / ignore box...
    ## at the moment need to click add to update...

    ## Also note need to make some of above into functions so can re-use code
    
    ###########################################################
    # Data parameters, cmd wrappers, etc.
    ###########################################################

    #------------------
    # data cmd and parameters
    #------------------
    def set_data_cmd(self):
        " run xrf.set_data command given data in self.data_par"
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
        
    def data_par_update(self):
        """ Updata self.data_par from the info in GUI components
        Note all the data is stored as strings to help cmd building"""

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
        
        # test
        #self.post_message(str(self.data_par))

        # this should ret true if E-range changed, ie reset
        # the xrf lines in range
        if (OldEmin != Emin) or (OldEmax != Emax):
            return True
        else:
            return False

    def data_par_display(self):
        """ Update components from self.data_par"""
        pass

    #------------------
    # plotting
    #------------------
    def plot_cmd(self):
        "run xrf.plot cmd"
        # update plot params
        self.plot_par_update()
        # build cmd
        xrf   = "%s.%s" % (self.data_par['grp'],self.data_par['node'])
        cmd_str = "xrf.plot(%s)" % xrf
        self.post_message(cmd_str)
        self.execLine(cmd_str)
        return

        
    def plot_par_update(self):
        """Update self.plot_par from components"""
        pass
    
    def plot_par_display(self):
        """Update components from self.plot_par"""
        pass

    #------------------
    # Fit pars
    #------------------
    def fit_params_update(self):
        pass
    def fit_params_display(self):
        pass

    #------------------
    # Bgr pars
    #------------------
    def bgr_params_update(self):
        pass
    def default_bgr_params(self):
        pass
    
    #------------------
    # peak pars
    #------------------
    def add_pk(self):
        pass
    def pk_params_update(self):
        pass
    def delete_pk(self):
        pass

    #------------------
    # calc/fit args
    #------------------
    def fit_args_update(self):
        pass
    def fit_args_display(self):
        pass

    #------------------
    # general
    #------------------
    def update_params(self):
        pass
    def display_params(self):
        pass
    def write_params_to_tdl(self):
        pass
    def read_params_from_tdl(self):
        pass
    def update_xrf(self):
        pass
    def plot(self):
        pass
    def calc(self):
        pass
    def fit(self):
        pass
    def fit_scan(self):
        pass

##################################################
if __name__ == '__main__':
    app = model.Application(wxXRF)
    app.MainLoop()
