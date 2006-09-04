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
# - look at execute line stuff and working on set data button.
# - next to work on is plotting....
#
############################################################################

from PythonCard import model
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

    ## Utilities
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
        txt = "\n---\n%s" % mess 
        #self.components.Messages.text = txt
        self.components.Messages.appendText(txt)
        return

    #def on_gainFocus(self,event):
    #    self.init_tdl()
    # I think this is exectued any time anything here gains focus?

    ###########################################################
    #             Components and Events                       #
    ###########################################################

    #--------------
    # GrpUpdate - button
    #--------------
    def on_GrpUpdate_mouseClick(self, event):
        "use this to update incase the list has changed"
        self.init_GrpItems()

    #--------------
    # GrpItems - choice
    #--------------
    ### Initialize the menu    
    def init_GrpItems(self):
        grp = self.components.Grp.text
        if len(grp) == 0:
            grp = self.getValue('_builtin.data_group')
        self.components.Grp.items = self.listGroups()
        self.components.Grp.text = grp
        print "text = %s" % (self.components.Grp.text)

    def on_Grp_select(self, event):
        "Re-init Node list"
        self.init_NodePfxItems()

    #------------------
    # NodePrefix - choice
    #------------------
    ### Initialize the menu    
    def init_NodePfxItems(self):
        """use the group thats selected in the group menu
        Note we could get smart here and only chose nodes that
        are actually xrf objects...."""
        grp = self.components.Grp.text  # should this be selected?
        if len(grp) == 0: grp = None
        self.components.NodePfx.items = self.listDataGroup(grp)

    def on_NodePfxmouseDown(self):
        """update"""
        self.init_NodePfxItems()

    #------------------
    # UpdateData - button
    #------------------
    def on_UpdateData_mouseClick(self,event):
        update_xrf_lines = self.data_par_update()
        if update_xrf_lines:
            self.init_XrfLineItems()
        xrf = "%s.%s" % (self.data_par['grp'],self.data_par['node'])
        det = self.data_par['mcas']
        total =str(self.data_par['total'])
        align = str(self.data_par['align'])
        tau = self.data_par['mcas']
        cmd_str = "xrf.set_data(%s,detectors=%s,total=%s,align=%s,tau=%s)" % \
                  (xrf,det,total,align,tau)
        self.post_message(cmd_str)
        self.execLine(cmd_str)

    #------------------
    # Plot - button 
    #------------------
    def on_Plot_mouseClick(self,event):
        do_xrf = self.data_par_update()
        if do_xrf: self.init_XrfLineItems()

    #------------------
    # XrfLine - choice
    #------------------
    ### Initialize the menu    
    def init_XrfLineItems(self):
        """use the group thats selected in the group menu
        Note we could get smart here and only chose nodes that
        are actually xrf objects...."""
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
        "Re-init Node list"
        #print event.target
        #line = self.components.XrfLine.selected
        line = self.components.XrfLine.stringSelection
        en = xrf_lookup.lookup_xrf_line(line)
        if en:
            en = str(en)
        else:
            en = ''
        self.components.PkLbl.text = line
        self.components.PkEn.text = en

    #------------------
    # XrfLineAdd - button
    #------------------
    def on_XrfLineAdd_mouseClick(self,event):
        det = self.components.DetSelect.stringSelection
        lbl=self.components.PkLbl.text.strip()
        #print 'hello', lbl
        if len(lbl) > 0:
            en=self.components.PkEn.text.strip()
            list = self.components.PkParams.items
            list.append([det,lbl,en,'0','0','0','0','0','0','0'])
            self.components.PkParams.items = list
            #print 'hello2',list, self.components.PkParams.items

    ###########################################################
    # Data                       
    ###########################################################

    ## data    
    def en_range_select(self):
        """ Select energy Range from plot"""
        pass
    
    def data_par_update(self):
        """ Updata self.data_par from components"""
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
        self.data_par['total'] = total
        #align
        align = self.components.Align.checked
        self.data_par['align'] = align
        #scan
        scan = self.components.AutoIncCheck.checked
        self.data_par['scan'] = scan
        #sc_start
        sc_start = self.components.AutoIncStart.text
        self.data_par['sc_start'] = int(sc_start)
        #sc_stop
        sc_stop = self.components.AutoIncStop.text
        self.data_par['sc_stop'] = int(sc_stop)
        #sc_inc
        sc_inc = self.components.AutoIncInc.text
        self.data_par['sc_inc'] = int(sc_inc)
        #sc_step
        sc_step = self.components.NodeInc.value  #??????
        self.data_par['sc_step'] = int(sc_step)
        #mcas
        mcas = self.components.McaList.text
        self.data_par['mcas'] = str(self.str_to_list(mcas,conv=int))
        #self.data_par['mcas'] = mcas
        #mca taus
        mca_taus = self.components.TauList.text
        self.data_par['mca_taus'] = str(self.str_to_list(mca_taus,conv=float))
        #self.data_par['mca_taus'] = mca_taus
        # Emin/Emax
        Emin = float(self.components.Emin.text)
        Emax = float(self.components.Emax.text)
        self.data_par['emin'] = Emin
        self.data_par['emax'] = Emax

        self.post_message(str(self.data_par))

        # this should ret true if E-range changed, ie reset
        # the xrf lines in range
        return False
        
        
    def data_par_display(self):
        """ Update components from self.data_par"""
        pass

    ## plot
    def plot_par_update(self):
        """Update self.plot_par from components"""
        pass
    
    def plot_par_display(self):
        """Update components from self.plot_par"""
        pass

    ## fit pars
    def fit_params_update(self):
        pass
    def fit_params_display(self):
        pass

    ## bgr
    def bgr_params_update(self):
        pass
    def default_bgr_params(self):
        pass
    
    ##peaks
    def add_pk(self):
        pass
    def pk_params_update(self):
        pass
    def delete_pk(self):
        pass

    ##calc/fit
    def fit_args_update(self):
        pass
    def fit_args_display(self):
        pass

    ## overall
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
