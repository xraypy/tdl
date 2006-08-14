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


class wxXRF(model.Background):
    
    def on_initialize(self, event):
        # if you have any initialization
        # including sizer setup, do it here
        self.shell = None
        self.tdl = None
        self.init_tdl()
        # are there potential race conditions?
        # window management???

    def init_tdl(self,shell=None):
        #self.shell = shell
        self.shell = self.getParent().get_shell()
        self.tdl = self.shell.tdl
        
        #### Init Menus
        self.init_XrfGrpItems()
        self.init_XrfNodePrefixItems()

    ## Utilities
    def setValue(self,var_name,value):
        self.tdl.setVariable(var_name,value)

    def getValue(self,var_name):
        return self.tdl.getVariableValue(var_name)

    def getVariable(self,var_name):
        return self.tdl.getVariable(var_name)

    def execLine(self,line):
        return self.tdl.execute(line)

    #def on_gainFocus(self,event):
    #    self.init_tdl()
    # I think this is exectued any time anything here gains focus?

    ###########################################################
    #             Components and Events                       #
    ###########################################################

    #--------------
    # XrfGrpItems 
    #--------------
    
    ### Initialize the menu    
    def init_XrfGrpItems(self): 
        self.components.XrfGrp.items = self.tdl.symbolTable.listGroups()
        data_group = self.tdl.getVariableValue('_builtin.data_group')
        self.components.XrfGrp.text = data_group

    ### Events
    # note text enter just works, ie typing into the field sets
    # the value of .text

    #def on_XrfGrp_mouseDown(self, event):
    #    self.init_XrfGrpItems()
    
    #def on_XrfGrp_mouseDoubleClick(self, event):
    #    self.init_XrfGrpItems()

    def on_XrfGrp_select(self, event):
        "I dont think this is actually needed, should just work"
        #group = self.components.XrfGrp.text
        #print "selected group = %s" % (group)
        #print "text = %s" % (self.components.XrfGrp.text)
        self.init_XrfNodePrefixItems()

    def on_XrfGrp_mouseMiddleDown(self, event):
        "use this to update incase the list has changed"
        print "text = %s" % (self.components.XrfGrp.text)
        txt = self.components.XrfGrp.text
        self.init_XrfGrpItems()
        self.components.XrfGrp.text = txt

    #------------------
    # XrfXrfNodePrefix 
    #------------------
    
    ### Initialize the menu    
    def init_XrfNodePrefixItems(self):
        """use the group thats selected in the group menu
        Note we could get smart here and only chose nodes that
        are actually xrf objects...."""
        group = self.components.XrfGrp.text
        if len(group) == 0: group = None
        #data_group = self.tdl.getVariableValue('_builtin.data_group')
        self.components.XrfNodePrefix.items = self.tdl.symbolTable.listDataGroup(group)

    ### Events
    def on_XrfNodePrefix_mouseMiddleDown(self, event):
        "use this to update incase the list has changed"
        print "text = %s" % (self.components.XrfNodePrefix.text)
        txt = self.components.XrfNodePrefix.text
        self.init_XrfNodePrefixItems()
        self.components.XrfNodePrefix.text = txt



    #######################################################################


##################################################
if __name__ == '__main__':
    app = model.Application(wxXRF)
    app.MainLoop()
