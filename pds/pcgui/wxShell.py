#!/usr/bin/python
"""
GUI shell program

Authors/Modifications:
----------------------
* Tom Trainor (tptrainor@alaska.edu)


Notes:
-----
This should normally be run from the pds script
to set up paths correctly

Note this module fails with wx2.8 (works fine with wx2.6).
In the later versions hitting enter in a field acts like a tab
in the old version - so moved focus from one widget to another.
That resulted in the 'enter' key never getting picked up by the
Process_Event method.  This is due to a change in the default
behaviour of wx in the latter versions.  The fix is to edit the
following in the PythonCard module to add wx.WANTS_CHAR and
remove wx.TAB_TRAVERSAL from the style argument to
wx.Panel.__init__() and wx.Frame.__init__():

* model.py, line 671:
    'style | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WANTS_CHARS,'
* model.py, line 1345:
    'style=wx.WANTS_CHARS,'
* widget.py, line 427:
    'style=wx.WANTS_CHARS | wx.NO_FULL_REPAINT_ON_RESIZE)'

"""
#######################################################################

from   PythonCard import model, dialog
import wx, string, sys, os, time
from   wx import stc
from  .wxUtil import wxUtil
from  .app_menu import menuApps
from tdl import pds
from  ..shellutil import mod_import

#######################################################################

intro = None
debug = False
files = []
args  = None
rsrc_path = '.'

#######################################################################
class wxShell(model.Background,menuApps,wxUtil):

    def on_initialize(self, event):
        
        # including sizer setup, do it here
        self.setupSizers()

        #print "initialize"
        self.indent = 0
        self.prompt = 'pds>'
        #self.components.Prompt.text = self.prompt
        self.save_file = None
        
        self.isreading = False
        self.input_text = ''
        self.cmd_history=['']
        self.cmd_from_hist=False
        self.cmd_count = 1
        self.cmd_idx = 0

        #redir stdio
        #sys.stdin  = self.readline
        self.rsrc_path = rsrc_path
        self.shell = pds.shell.Shell(args=args,
                                     stdin=self,stdout=self,
                                     GUI='WXAgg',debug=debug)
        for f,warn in files:
            if os.path.exists(f) and os.path.isfile(f):
                self.shell.do_load(f)
            elif warn:
                print "\n  ***Cannot find file: %s" % f

        # run it, when done we've quit  
        self.run()
        sys.__stdout__.write("\nExiting \n")
        time.sleep(.5)
        self.close()
        sys.exit()

    ################################################
    def run(self,fname=''):
        ret = self.shell.loop()

    ################################################
    def get_shell(self):
        #print "hello from get shell"
        return self.shell

    ################################################
    # stdIO
    ################################################
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
            #time.sleep(.01)
            time.sleep(.025)
            wx.YieldIfNeeded()
            #wx.Yield()
        #print "here"
        input_text = str(self.input_text)
        self.input_text = ''
        #line = input_text.strip() + '\n'
        line = input_text
        
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
    
    """
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
    """
    
    def PostLineToShellText(self,line):
       self.components.ShellText.appendText(line)
       #self.components.ShellText.setInsertionPointEnd()

    ###########################################################
    #             Menus                                       #
    ###########################################################
    
    # note see app_menu.py for more...  

    def on_menuFileExit_select(self,event):
        print "QUIT"
        #self.exec_line("quit")
        #self.shell.do_quit()
        #self.close()
        sys.exit()

    def on_menuFile_CD_select(self,event):        
        cdir = self.eval_line("pwd()")
        result = dialog.directoryDialog(self, 'Open', cdir)
        #print result
        if result.accepted:
            dir = result.path
            if os.path.abspath(dir) != os.path.abspath(cdir):
                dir = dir.replace("\\","\\\\")
                line = "cd('%s')" % dir
                self.eval_line(line)
        else:
            self.post_message("Change Dir cancelled.")
            return

    def on_menuFile_Restore_select(self,event):        
        cdir = self.eval_line("pwd()")
        result = dialog.fileDialog(self, 'Open', cdir, '',"*")
        if result.accepted:
            path        = result.paths[0]
            dir,fname   = os.path.split(path)
            if os.path.abspath(dir) != os.path.abspath(cdir):
                dir = dir.replace("\\","\\\\")
                line = "cd('%s')" % dir
                self.eval_line(line)
        else:
            self.post_message("File selection cancelled.")
            return
        
        #print dir, fname
        self.save_file = path
        line = "restore %s" % path  
        self.exec_line(line)

    def on_menuFile_SaveAs_select(self,event):
        cdir = self.eval_line("pwd()")
        result = dialog.fileDialog(self, 'Save', cdir, '',"*")
        if result.accepted:
            path        = result.paths[0]
            dir,fname   = os.path.split(path)
            if os.path.abspath(dir) != os.path.abspath(cdir):
                dir = dir.replace("\\","\\\\")
                line = "cd('%s')" % dir
                self.eval_line(line)
        else:
            self.post_message("File selection cancelled.")
            return
        
        #print dir, fname
        self.save_file = path
        line = "save %s" % path  
        self.exec_line(line)

    def on_menuFile_Save_select(self,event):        
        if self.save_file == None:
            self.exec_line("save")
        else:
            line = "save %s" % self.save_file  
            self.exec_line(line)
                
    def on_menuShellStart_select(self, event):
        # load shell is defined in model.Background (model.py)
        # it starts an instance of pycrust
        self.loadShell()
        if self.application.shell is not None:
            self.application.shellFrame.visible = not self.application.shellFrame.visible

    def on_menuOptionsEditSiteStartup_select(self,event):
        self.exec_line("import os")
        self.exec_line("edit os.path.join(__pds__.__path__,'startup.pds')")

    def on_menuOptionsEditHomeStartup_select(self,event):
        self.exec_line("import os")
        self.exec_line("edit os.path.join(__home__,'.pds')")
        
    def on_menuHelpUse_select(self, event):
        import wxShellHelp
        wxShellHelp = mod_import(wxShellHelp)
        dir       = os.path.dirname(wxShellHelp.__file__)
        filename  = os.path.join(dir,'wxShellHelp.rsrc.py')
        wxShellHelp = wxShellHelp.wxShellHelp
        self.wxShellHelp = model.childWindow(self,wxShellHelp,
                                             filename=filename)
        self.wxShellHelp.position = (200, 5)
        self.wxShellHelp.visible = True

    def on_menuHelpDocumentation_select(self,event):
        self.exec_line("web 'http://cars9.uchicago.edu/iffwiki/tdl/Pds'")

    ###########################################################
    #             EVENTS                                      #
    ###########################################################

    def on_ShellCmd_keyDown(self, event):
        #print "keyDown", event.keyCode, "\n"
        #print event.shiftDown, event.controlDown, event.altDown
        self.Process_Event(event)

    def Process_Event(self, event):
        """
        """
        #print dir(event)
        keyCode = event.GetKeyCode()
        #print keyCode
        
        #if keyCode == wx.WXK_CONTROL:
        # 372 (old) 310 (new) is enter on the numeric keypad
        if (keyCode == wx.WXK_RETURN) or (keyCode == wx.WXK_NUMPAD_ENTER):
            self.input_text = self.components.ShellCmd.text
            #self.input_text = self.components.ShellCmd2.text
            self.input_text = self.input_text + '\n'
            #tmp = self.prompt + self.input_text
            tmp = self.input_text
            self.PostLineToShellText(tmp)
            self.components.ShellCmd.text = ''
            self.isreading = False
        # uparrow, 317 in old wx, 315 in newer: 
        elif keyCode == wx.WXK_UP: 
            self.cmd_idx=self.cmd_idx+1
            if self.cmd_idx > self.cmd_count-1:
                self.cmd_idx = self.cmd_count -1
            self.components.ShellCmd.text = self.cmd_history[self.cmd_idx]
            #self.components.ShellCmd.setInsertionPointEnd()
            self.cmd_from_hist=True
        # downarrow, 319 in old wx, 317 in new
        elif keyCode == wx.WXK_DOWN: 
            self.cmd_idx=self.cmd_idx+-1
            if self.cmd_idx<0:
                self.cmd_idx=0
                self.cmd_from_hist=False
            else:
                self.cmd_from_hist=True             
            self.components.ShellCmd.text = self.cmd_history[self.cmd_idx]
            #self.components.ShellCmd.setInsertionPointEnd()
            
        else:
            event.Skip()

    #####################################################
    def setupSizers( self ):
        p1 = 40  # 80
        p2 = 0   # 1
        p3 = 1   # 1
        p4 = 15  # 20
        p5 = 10  # 50
        p6 = 0   #1   # 2 
        p7 = 30  # 50 
        comp = self.components
        # Create base sizers
        base_sizer = wx.BoxSizer( wx.VERTICAL )
        base_sizer_V = wx.BoxSizer( wx.VERTICAL )

        # here are the sizers for the shell and sep bar
        # note 80:1 ratio of vertical dims
        shell_sizer_V = wx.BoxSizer(wx.VERTICAL)
        shell_sizer_V.Add( comp.ShellText, p1, wx.ALL | wx.EXPAND | wx.ALIGN_LEFT, 0 )
        shell_sizer_V.Add( comp.Sep1, p2, wx.ALL | wx.EXPAND | wx.ALIGN_LEFT, 0 )
        
        # here make horz sizer for prompt and cmd line
        cmd_sizer_H = wx.BoxSizer( wx.HORIZONTAL )
        cmd_sizer_H.Add( comp.Prompt, p3, wx.ALL | wx.EXPAND | wx.ALIGN_LEFT, 0 )
        cmd_sizer_H.Add( comp.ShellCmd, p4, wx.ALL | wx.EXPAND | wx.ALIGN_LEFT, 0 )
        
        # Now add both of these to the base vert sizer
        base_sizer_V.Add(shell_sizer_V, p5, wx.ALL | wx.EXPAND,0 )
        base_sizer_V.Add(cmd_sizer_H, p6, wx.ALL | wx.EXPAND,0)

        # Now add this to the base sizer        
        base_sizer.Add(base_sizer_V, p7, wx.ALL | wx.EXPAND ,5)

        # Magic
        base_sizer.Fit( self )
        base_sizer.SetSizeHints( self )
        self.panel.SetSizer( base_sizer )
        self.panel.SetAutoLayout( 1 )
        self.panel.Layout()
        self.visible = True

    #####################################################
    """
    # some tests

    #def on_ShellCmd_textUpdate(self, event):
    #    print "textUpdate", event.keyCode, "\n"
    #def on_ShellCmd_keyUp(self, event):
    #    print "keyUp", event.keyCode, "\n"
    #def on_ShellCmd_keyPress(self, event):
    #    print "keyPress", event.keyCode, "\n"

    # Use below to test other key bindings
    def _init_keys(self):
        #wx.EVT_KEY_DOWN(self.components.ShellCmd, self._hit_enter)
        #wx.EVT_KEY_UP(self.components.ShellCmd, self._hit_enter)
        #wx.EVT_CHAR(self.components.ShellCmd, self._hit_enter)
        #wx.EVT_CHAR(self.components.ShellText, self._hit_enter)
        #wx.EVT_KEY_DOWN(self.components.ShellCmd2, self._hit_enter)

    def _hit_enter(self,event):
        #keyCode = event.GetKeyCode()
        #print keyCode
        #print dir(self)
        #print self.WindowStyle
        #print self.WindowStyleFlag
        # print wx.WXK_RETURN
        #sys.__stdout__.write(str(event.keyCode))
        #sys.__stdout__.write(str(wx.WXK_RETURN))
        self.Process_Event(event)

    # global key down event
    #def on_keyDown(self,event):
    #    keycode = event.GetKeyCode()
    #    #print 'on keydown'
    #    #print keycode
    #    self.Process_Event(event)
    """

################################################################
if __name__ == '__main__':
    app = model.Application(wxGui)
    app.MainLoop()

