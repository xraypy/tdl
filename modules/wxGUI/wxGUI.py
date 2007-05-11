#!/usr/bin/python
############################################################################
# Tom Trainor (fftpt@uaf.edu)
# This should get run from tdl.py to set up paths correctly
#
# --------------
# Modifications:
# --------------
# 4-16-2006 T2:
# - handle startup files
#
#
############################################################################

from PythonCard import model
import wx, string, sys, os, time
from wx import stc

tdl   = None
libs  = []
intro = None
debug = False
files = []
sys_vars = {}
rsrc_path = '.'

#####################################################

class wxGUI(model.Background):


    def on_initialize(self, event):

        # including sizer setup, do it here
        self.setupSizers()

        print "initialize"
        self.indent = 0
        self.prompt = 'tdl>'
        self.components.Prompt.text = self.prompt
        self.isreading = False
        self.input_text = ''
        self.cmd_history=['']
        self.cmd_from_hist=False
        self.cmd_count = 1
        self.cmd_idx = 0
        #redir stdio
        #sys.stdin  = self.readline
        #self.tdl = tdl
        self.rsrc_path = rsrc_path
        self.shell = tdl.shell(libs=libs,stdin=self,stdout=self,
                               GUI='WXAgg',debug=debug)
        for var in sys_vars.keys():
            self.shell.tdl.setVariable(var,sys_vars[var])
        for f,warn in files:
            if os.path.exists(f) and os.path.isfile(f):
                self.shell.tdl_execute("load('%s')" % f)
            elif warn:
                print "\n  ***Cannot find tdl file: %s" % f
        self.run_tdl()
        # when done we've quit  
        sys.__stdout__.write("\nExiting TDL\n")
        time.sleep(.5)
        self.close()
        sys.exit()



    ################################################
    def run_tdl(self,fname=''):
        ret = self.shell.cmdloop()

    ################################################
    # stdIO
    def write(self, text):
        #self.PostLineToShellText(text)
        self.components.ShellText.appendText(text)

    def flush(self):
        self.input_text = ''
        
    def readline(self):
        return(self.ReadInputLine())
    
    def raw_input(self, prompt=''):
        """Return string based on user input."""
        if prompt:
            self.UpdatePrompt(prompt)
        line = self.ReadInputLine()         
        return(line)

    #########################################
    def UpdatePrompt(self,prompt):
        self.prompt = prompt
        self.components.Prompt.text = prompt

    def ReadInputLine(self):
        self.isreading = True
        while self.isreading:
            time.sleep(.01)
            wx.YieldIfNeeded()
        input_text = str(self.input_text)
        self.input_text = ''
        line = input_text.strip() + '\n'
        #history
        if len(line)>1:
            # element zero always blank line
            # element 1 most recent etc..
            if self.cmd_from_hist:
                del self.cmd_history[self.cmd_idx]
            else:
                self.cmd_count = self.cmd_count + 1
            self.cmd_history.insert(1,line[:-1])
        self.cmd_idx = 0
        self.cmd_from_hist=False
        return (line)
    
    def PostLineToShellText(self,line):
        # is this too slow????
        # get the last line
        txtlines  = string.split(self.components.ShellText.text,'\n')
        if txtlines:
            last_line = txtlines[len(txtlines)-2]
        else:
            last_line = ''
                
        # auto-indent block
        if len(string.strip(last_line)) == 0:
            self.indent = self.indent - 4
        else:
            tmp = string.lstrip(last_line)
            indent = len(last_line) - len(tmp)
            if indent != self.indent:
                self.indent = indent
            if last_line[-1] == ':':
                self.indent = self.indent + 4
        if self.indent < 0: self.indent = 0
        padding = " " * self.indent
        line = padding + line
        self.components.ShellText.appendText(line)
    
    #def PostLineToShellText(self,line):
    #   self.components.ShellText.appendText(line)

    def get_shell(self):
        #print "hello from get shell"
        return self.shell
    
    ###########################################################
    #             Menus                                       #
    ###########################################################

    def on_menuFileExit_select(self,event):        
        self.close()
        sys.exit()

    def on_menuShellStart_select(self, event):
        # load shell is defined in model.Background (model.py)
        # it starts an instance of pycrust
        self.loadShell()
        if self.application.shell is not None:
            self.application.shellFrame.visible = not self.application.shellFrame.visible

    def on_menuWindowXRF_select(self, event):
        from wxXRF import wxXRF
        # the child has to get self.shell from the parent
        filename = os.path.join(self.rsrc_path,'wxGUI','wxXRF.rsrc.py')
        self.XrfWindow = model.childWindow(self, wxXRF,filename=filename)
        self.XrfWindow.position = (200, 5)
        self.XrfWindow.visible = True
    
    ###########################################################
    #             EVENTS                                      #
    ###########################################################

    def on_ShellCmd_keyDown(self, event):
        # print "keyPress", event.keyCode,
        # print event.shiftDown, event.controlDown, event.altDown
        #sys.__stdout__.write(str(event.keyCode))

        keyCode = event.keyCode

        if keyCode == 317: # uparrow
            self.cmd_idx=self.cmd_idx+1
            if self.cmd_idx > self.cmd_count-1:
                self.cmd_idx = self.cmd_count -1
            self.components.ShellCmd.text = self.cmd_history[self.cmd_idx]
            #self.components.ShellCmd.setInsertionPointEnd()
            self.cmd_from_hist=True

        if keyCode == 319: # downarrow
            self.cmd_idx=self.cmd_idx+-1
            if self.cmd_idx<0:
                self.cmd_idx=0
                self.cmd_from_hist=False
            else:
                self.cmd_from_hist=True             
            self.components.ShellCmd.text = self.cmd_history[self.cmd_idx]
            #self.components.ShellCmd.setInsertionPointEnd()
            
        elif keyCode == wx.WXK_RETURN:
            self.input_text = self.components.ShellCmd.text
            self.input_text = self.input_text + '\n'
            #tmp = self.prompt + self.input_text
            tmp = self.input_text
            self.PostLineToShellText(tmp)
            self.components.ShellCmd.text = ''
            self.isreading = False

        else:
            event.skip()


####
    def setupSizers( self ):
        comp = self.components
        # Create base sizers
        base_sizer = wx.BoxSizer( wx.VERTICAL )
        base_sizer_V = wx.BoxSizer( wx.VERTICAL )

        # here are the sizers for the shell and sep bar
        # note 80:1 ratio of vertical dims
        shell_sizer_V = wx.BoxSizer(wx.VERTICAL)
        shell_sizer_V.Add( comp.ShellText, 80, wx.ALL | wx.EXPAND | wx.ALIGN_LEFT, 0 )
        shell_sizer_V.Add( comp.Sep1, 1, wx.ALL | wx.EXPAND | wx.ALIGN_LEFT, 0 )
        
        # here make horz sizer for prompt and cmd line
        cmd_sizer_H = wx.BoxSizer( wx.HORIZONTAL )
        cmd_sizer_H.Add( comp.Prompt, 1, wx.ALL | wx.EXPAND | wx.ALIGN_LEFT, 0 )
        cmd_sizer_H.Add( comp.ShellCmd, 20, wx.ALL | wx.EXPAND | wx.ALIGN_LEFT, 0 )
        
        # Now add both of these to the base vert sizer
        base_sizer_V.Add(shell_sizer_V,50, wx.ALL | wx.EXPAND,0 )
        base_sizer_V.Add(cmd_sizer_H, 2, wx.ALL | wx.EXPAND,0)

        # Now add this to the base sizer        
        base_sizer.Add(base_sizer_V, 50, wx.ALL | wx.EXPAND ,5)

        # Magic
        base_sizer.Fit( self )
        base_sizer.SetSizeHints( self )
        self.panel.SetSizer( base_sizer )
        self.panel.SetAutoLayout( 1 )
        self.panel.Layout()
        self.visible = True

####


################################################################
if __name__ == '__main__':
    app = model.Application(wxGUI)
    app.MainLoop()

