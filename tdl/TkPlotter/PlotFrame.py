#!/usr/bin/python
##
## MPlot PlotFrame: a Tk.Frame for 2D line plotting, using matplotlib
##

import os
import time
import Tkinter as Tk
import Pmw
from PlotPanel import PlotPanel
import types

class Command:
    " Generic Command execution (simpler than lambdas)"
    def __init__(self, func, *args, **kw):
        self.func = func
        self.args = args
        self.kw   = kw
    def __call__(self, *args, **kw):
        args = self.args + args
        kw.update(self.kw)
        apply(self.func,args, kw)

class PlotFrame(Pmw.MegaWidget):
    """
    MatPlotlib 2D plot as a Tk.Frame, using PlotPanel
    """

    help_msg =  """MPlot PlotFrame quick help:

 Left-Click:   to display X,Y coordinates
 Left-Drag:    to zoom in on plot region
 Right-Click:  display popup menu with choices:
                Zoom out 1 level       (that is, to previous view)
                Zoom all the way out   (to full data range)
                --------------------
                Configure Plot
                Save Plot Image

Also, these key bindings can be used
(For Mac OSX, replace 'Ctrl' with 'Apple'):

  Ctrl-S:     save plot image to file
  Ctrl-C:     copy plot image to clipboard
  Ctrl-K:     Configure Plot 
  Ctrl-Q:     quit
"""

    about_msg =  """TkPlot  version 0.1 
att Newville <newville@cars.uchicago.edu>"""

    usecommandarea=1
    appname       = 'TkPlot'
    appversion    = '0.1'
    copyright     = 'Copyright 2006, The University of Chicago'
    contactname   = 'Matt Newville'
    contactphone  = '(630) 252-0431'
    contactemail  = 'newville@cars.uchicago.edu'
    contactweb    = 'http://cars.uchicago.edu/~newville'
    frameHeight   = 375
    frameWidth    = 660
    padx            = 5
    pady            = 5
    usecommandarea  = 0
    balloonhelp     = 0
    busyCursor = 'watch'
    
    def __init__(self, root=None, size=(650,375), exit_callback=None, **kwds):
        self.exit_callback = exit_callback
        if root is None: root = Tk.Tk()
        self.root = root
        optiondefs = (
            ('padx',           1,                   Pmw.INITOPT),
            ('pady',           1,                   Pmw.INITOPT),
            ('framewidth',     1,                   Pmw.INITOPT),
            ('frameheight',    1,                   Pmw.INITOPT),
            ('usecommandarea', self.usecommandarea, Pmw.INITOPT))
        self.defineoptions(kwds, optiondefs)
        self.initializeTk(self.root)
        Pmw.initialise(self.root)
        self.root.title(self.appname)
        self.root.geometry('%dx%d' % (self.frameWidth, self.frameHeight))
        Pmw.MegaWidget.__init__(self, parent=self.root)

        # create the interface
        self.createBalloon()
        self.createAboutBox()
        self.createMenuBar()
        self.createStatusBar()
        self.createMainFrame()
        
        self.busyWidgets = ( self.root, )
        self.createInterface()
       
        # create a table to hold the cursors for
        # widgets which get changed when we go busy
        self.preBusyCursors = None
        
        # pack the container and set focus to ourselves
        self._hull.pack(side=Tk.TOP, fill=Tk.BOTH, expand=Tk.YES)
        self.focus_set()
        Tk.Wm.protocol(self.root, 'WM_DELETE_WINDOW', self.onExit)
        
    def onExit(self,event=None):
        ' called when the PlotFrame is closed'
        # print 'PlotFrame onExit '
        if self.exit_callback is not None: self.exit_callback()
        self.root.destroy()


   #######
    def initializeTk(self, root):
        bg    = '#F8F8F4'
        altbg = '#FEFEF0'
        ent   = '#F8F8F8'
        sel   = bg # 'lemonchiffon2' 
        bsel  = bg # 'lemonchiffon2'
        root.tk_setPalette(bg)
        root.option_add('*foreground',                 'black')
        root.option_add('*EntryField.Entry.background', ent)
        root.option_add('*Entry.background',            ent)        
        root.option_add('*Listbox*background',          ent)
        root.option_add('*StatusBar.Entry.background', bg)
        root.option_add('*MenuBar.background',          bg)
        root.option_add('*Listbox*selectBackground',    sel)
        root.option_add('*Checkbutton.selectColor',     'white')
        root.option_add('*Checkbutton.background',      altbg)
        root.option_add('*Radiobutton.selectColor',     bsel)        
        root.option_add('*Checkbutton.background',      bg)
        root.option_add('*Checkbutton.foreground',      'black')
        root.option_add('*Checkbutton.activecolor',     'white')
        root.option_add('*Checkbutton.highlightcolor',   'white')
        root.option_add('*Checkbutton.selectcolor',     'blue')        
        root.option_add('*Radiobutton.selectColor',     altbg)
        root.option_add('*Font', 'Arial 11 bold')

    def busyStart(self, newcursor=None):
        if not newcursor:
            newcursor = self.busyCursor
        newPreBusyCursors = {}
        for component in self.busyWidgets:
            newPreBusyCursors[component] = component['cursor']
            component.configure(cursor=newcursor)
            component.update_idletasks()
        self.preBusyCursors = (newPreBusyCursors, self.preBusyCursors)
        
    def busyEnd(self):
        if not self.preBusyCursors:
            return
        oldPreBusyCursors = self.preBusyCursors[0]
        self.preBusyCursors = self.preBusyCursors[1]
        for component in self.busyWidgets:
            try:
                component.configure(cursor=oldPreBusyCursors[component])
            except KeyError:
                pass
            component.update_idletasks()
              
    def createAboutBox(self):
        Pmw.aboutversion(self.appversion)
        Pmw.aboutcopyright(self.copyright)
        ss = '%s <%s>\n%s' % (self.contactname, self.contactemail,self.contactweb)
        Pmw.aboutcontact(ss)
        self.about = Pmw.AboutDialog(self._hull, 
                                     applicationname=self.appname)
        self.about.withdraw()
        return None
       
    def showAbout(self):
        # Create the dialog to display about and contact information.
        self.about.show()
        self.about.focus_set()
       
    def toggleBalloon(self):
	if self.toggleBalloonVar.get():
            self.__balloon.configure(state = 'both')
	else:
            self.__balloon.configure(state = 'status')

    def createBalloon(self):
        # Create the balloon help manager for the frame.
        # Create the manager for the balloon help
        self.__balloon = self.createcomponent('balloon', (), None,
                                              Pmw.Balloon, (self._hull,))

    def balloon(self):
        return self.__balloon

    def createMainFrame(self):
        # Create data area where data entry widgets are placed.
        self.main =  Tk.Frame(self._hull, relief=Tk.GROOVE,bd=1)
        self.main.pack(side=Tk.TOP, fill=Tk.BOTH, expand=Tk.YES,
                           padx=self['padx'], pady=self['pady'])

    def createStatusBar(self):
        # Create the status bar area for help, status messages.
        self.statusBar = Tk.Frame(self._hull, relief=Tk.SUNKEN,bd=1)
        self.__balloon.configure(statuscommand = self.setStatusText)
        self.statusBar.pack(side=Tk.BOTTOM, expand=Tk.NO, fill=Tk.X)
                   
    def bind(self, child, balloonHelpMsg, statusHelpMsg=None):
        # Bind a help message and/or status message to a widget.
        self.__balloon.bind(child, balloonHelpMsg, statusHelpMsg)

    def buttonBox(self):
        # Retrieve the button box.
        return self.__buttonBox

    def buttonAdd(self, buttonName, helpMessage=None,
                  statusMessage=None, **kw):
        # Add a button to the button box.
        newBtn = self.__buttonBox.add(buttonName)
        newBtn.configure(kw)
        if helpMessage:
             self.bind(newBtn, helpMessage, statusMessage)
        return newBtn

    def plot(self,x,y,**kw):
        """plot after clearing current plot """        
        self.plotpanel.plot(x,y,**kw)
        
    def oplot(self,x,y,**kw):
        """generic plotting method, overplotting any existing plot """
        self.plotpanel.oplot(x,y,**kw)

    def show_map(self,z, **kw):
        """generic plotting method, overplotting any existing plot """
        self.plotpanel.show_map(z,**kw)

    def update_line(self,t,x,y,**kw):
        """overwrite data for trace t """
        self.plotpanel.update_line(t,x,y,**kw)

    def set_xylims(self,xylims,**kw):
        """overwrite data for trace t """
        self.plotpanel.set_xylims(xylims,**kw)

    def get_xylims(self):
        """overwrite data for trace t """
        return self.plotpanel.get_xylims()

    def clear(self):
        """clear plot """
        self.plotpanel.clear()

    def unzoom_all(self,event=None):
        """zoom out full data range """
        self.plotpanel.unzoom_all(event=event)

    def unzoom(self,event=None):
        """zoom out 1 level, or to full data range """
        self.plotpanel.unzoom(event=event)
        
    def set_title(self,s):
        "set plot title"
        self.plotpanel.set_title(s)
        
    def set_xlabel(self,s):
        "set plot xlabel"        
        self.plotpanel.set_xlabel(s)

    def set_ylabel(self,s):
        "set plot xlabel"
        self.plotpanel.set_ylabel(s)        

    def save_figure(self,event=None):
        """ save figure image to file"""
        self.plotpanel.save_figure(event=event)

    def configure(self,event=None):
        self.plotpanel.configure(event=event)

    ####
    ## create GUI 
    ####
    def createInterface(self):
        self.plotpanel = PlotPanel(self.main)
        self.plotpanel.messenger = self.setStatusText
        self.plotpanel.report_xy = self.setXYtext
        self.plotpanel.pack(side=Tk.TOP,expand=Tk.YES,fill=Tk.BOTH)

        self.status_y  = Tk.Label(self.statusBar,text='Y: ',width=25, anchor=Tk.W)
        self.status_x  = Tk.Label(self.statusBar,text='X: ',width=25, anchor=Tk.W)
        self.statusmsg = Tk.Label(self.statusBar,text=''  ,width=60, anchor=Tk.W)

        self.status_y.pack(side=Tk.RIGHT)
        self.status_x.pack(side=Tk.RIGHT)
        self.statusmsg.pack(side=Tk.RIGHT)
 
    def setStatusText(self,val):
        # print 'status message ', val
        self.statusmsg.configure(text=val)

    def setXYtext(self,val):
        x,y = val.split('&')
        self.status_x.configure(text="X: %12.6g" % float(x))
        self.status_y.configure(text="Y: %12.6g" % float(y))
        self.cursor = (float(x),float(y))

    def createMenuBar(self):
        self.menuBar = Pmw.MenuBar(self._hull,
                                   hull_relief=Tk.RAISED,
                                   hull_borderwidth=1,
                                   balloon=self.balloon())

        self.menuBar.pack(fill=Tk.X)
        self.menuBar.addmenu('Help', 'About %s' % self.appname, side='right')
        self.menuBar.addmenu('File', 'File commands and Quit')
        self.menuBar.addmenu('Configure', 'Plot Options')
       
        self.menuBar.addmenuitem('Help', 'command',
                                 'Get information on application',
                                 label='About...', command=self.showAbout)
        self.toggleBalloonVar = Tk.IntVar()
        self.toggleBalloonVar.set(0)

        self.help = self.balloon()
        self.help.configure(state = 'status')
        self.menuBar.addmenuitem('Help', 'checkbutton',
                                 'Toggle balloon help',
                                 label='Balloon help',
                                 variable = self.toggleBalloonVar,
                                 command=self.toggleBalloon)

        self.menuBar.addmenuitem('File', 'command', 'Save PNG Image of Plot',
                                label='Save ',
                                command=self.save_figure)

        self.menuBar.addmenuitem('File', 'command', 'Quit this application',
                                label='Quit ',
                                command=self.onExit)

        self.menuBar.addmenuitem('Configure', 'command', 'zoom out',
                                label='Zoom Out ',
                                command=self.unzoom_all)

        self.menuBar.addmenuitem('Configure', 'command', 'configure',
                                label='Plot Options ',
                                command=self.configure)

