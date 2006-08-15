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
#
# - go back to rsrc file and simplify component names
# - make a single var menu update function (groups and nodes)
# - maybe add a sync/update button rather than the middle mouse d-click.
#
############################################################################

from PythonCard import model
import xrf_lookup

class wxXRF(model.Background):
    
    def on_initialize(self, event):
        # if you have any initialization
        # including sizer setup, do it here
        # are there potential race conditions?
        # window management???
        self.shell = None
        self.tdl = None
        self.init_tdl()

        # parameters
        self.NewData = True


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
        return self.tdl.execute(line)

    #def on_gainFocus(self,event):
    #    self.init_tdl()
    # I think this is exectued any time anything here gains focus?

    ###########################################################
    #             Components and Events                       #
    ###########################################################

    #--------------
    # GrpUpdate 
    #--------------
    def on_GrpUpdate_mouseClick(self, event):
        "use this to update incase the list has changed"
        self.init_GrpItems()

    #--------------
    # GrpItems 
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
    # NodePrefix 
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
    # UpdateData 
    #------------------
    def on_UpdateData_mouseClick(self,event):
        self.init_XrfLineItems()

    #------------------
    # XrfLine 
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
        if self.components.Emin:
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
        

##################################################
if __name__ == '__main__':
    app = model.Application(wxXRF)
    app.MainLoop()
