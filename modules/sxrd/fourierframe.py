"""
Functions and classes used in pisurf for plotting electron density maps

Authors/modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""

################################################################################
import wx
import os
import numpy as Num
import matplotlib
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure

################################################################################
class wxFourierFrame(wx.Frame):
    def __init__(self, parent, rho,cell,zmin, title, size):
        wx.Frame.__init__(self, parent, title= title, size=size)
        
        self.CreateStatusBar()
        self.menubar = wx.MenuBar()
        self.file = wx.Menu()
        menusavefig = self.file.Append(wx.ID_ANY, '&Save',\
                                       "Save the current figure") 
        menuexit = self.file.Append(wx.ID_EXIT, \
                                    '&Exit', u"Close \u03c0-surf Fourier Frame")
        self.menubar.Append(self.file, '&File')
        self.SetMenuBar(self.menubar)
        
        self.plot = PlotPanel(self)
        self.plotsizer = wx.BoxSizer()
        self.plotsizer.Add(self.plot, 1, wx.EXPAND)

        self.choice = ChoicePanel(self,rho,cell,zmin)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.plotsizer, 1, wx.EXPAND)
        self.sizer.Add(self.choice, 0,wx.EXPAND)
        self.SetSizer(self.sizer)
        self.SetMinSize((380, 400))
        self.Bind(wx.EVT_MENU, self._OnSave,menusavefig)
        self.Bind(wx.EVT_MENU, self._OnClose,menuexit)
    def _OnSave(self, event):
        dirname = ''
        dlg = wx.FileDialog(self, "Save the current Figure to a file",\
                            dirname, ".png", "*.png", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            self.plot.figure.savefig(filename)
            
        dlg.Destroy()

    def _OnClose(self, event):
        self.Close()
################################################################################          
class PlotPanel(wx.Panel):
    def __init__( self, parent):
        wx.Panel.__init__( self, parent)

        # initialize matplotlib stuff
        self.figure = Figure()
        self.canvas = FigureCanvasWxAgg( self, -1, self.figure )
        self.sizer = wx.BoxSizer()
        self.sizer.Add(self.canvas,1,wx.EXPAND)
        self.SetSizer(self.sizer)

################################################################################
class ChoicePanel(wx.Panel):
    def __init__( self, parent,rho,cell,zmin,pos=(0,500),size=(600,100)):
        wx.Panel.__init__( self, parent)
        self.figure = self.GetParent().plot.figure.add_subplot(111)
        self.view = int(3)
        wx.StaticBox(self, -1, 'Viewing Options', pos=(10,30), size=(350,90))
        wx.StaticText(self, label = 'viewing direction', pos=(30,50), \
                      size =(80,20))
        wx.StaticText(self, label = ' choose layer ', pos=(190,50), \
                      size = (100,20))
        self.viewingoptions = ['xy','xz','yz','z']
        self. viewchoice = wx.ComboBox(self, -1, \
                                       value=self.viewingoptions[self.view],\
                                       pos=(30,70), size=(80,20), \
                                       choices = self.viewingoptions, \
                                       style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.SelectView, self.viewchoice)
        self.DownButton = wx.Button(self, -1, ' < ', pos=(140,70), size=(40,20))
        self.Bind(wx.EVT_BUTTON, self.OnDown, self.DownButton)
        self.UpButton = wx.Button(self, -1, ' > ', pos=(300,70), size=(40,20))
        self.Bind(wx.EVT_BUTTON, self.OnUp, self.UpButton)

        self.layeroptions = ['sum']
        self.layer = int(0)
        self.layerchoice = wx.ComboBox(self, -1,\
                                       value=self.layeroptions[self.layer],\
                                       pos=(190,70), size=(100,20), \
                                       choices = self.layeroptions,\
                                       style=wx.CB_READONLY)
        
        self.Bind(wx.EVT_COMBOBOX, self.SelectLayer, self.layerchoice)

        self.rho = rho
        self.cell = cell
        self.zmin = zmin
        self.plot_rho()
       
    def SelectView(self, event):
        self.view = event.GetSelection()
        self.layer = 0
        self.layeroptions = ['sum']
        if self.view == 3:
            pass
        elif self.view == 0:
            for i in range(self.rho.shape[2]):
                self.layeroptions.append(str(i))
        elif self.view == 1:
            for i in range(self.rho.shape[1]):
                self.layeroptions.append(str(i))
        elif self.view == 2:
            for i in range(self.rho.shape[0]):
                self.layeroptions.append(str(i))
        self.layerchoice.Clear()
        for item in self.layeroptions:
            self.layerchoice.Append(item)
        self.layerchoice.SetSelection(0)
            
        self.plot_rho()
        
    def OnDown(self, event):
        layer = self.layerchoice.GetSelection()
        if layer > 0: layer = layer -1
        else: layer = len(self.layeroptions)-1
        self.layer = layer
        self.layerchoice.SetSelection(layer)
        self.plot_rho()
    
    def SelectLayer(self,event):
        self.layer = event.GetSelection()
        self.plot_rho()
        
    def OnUp(self, event):
        layer = self.layerchoice.GetSelection()
        if layer < len(self.layeroptions)-1: layer = layer +1
        else: layer = 0
        self.layer = layer
        self.layerchoice.SetSelection(layer)
        self.plot_rho()

    def plot_rho(self):
        rho = self.rho
        Max = Num.max(rho)
        Min = Num.min(rho)
        if self.layer == 0: plane = None
        else: plane = self.layer-1
        view = self.view
        cell = self.cell
        zmin = self.zmin
        fig = self.figure
        fig.clear()
        
        if view == 3:
            z = (Num.arange(rho.shape[2], dtype = float)/rho.shape[2]*cell[2])\
                +zmin
            image = Num.sum(rho, axis = 1)
            graph = Num.sum(image, axis = 0)
            fig.set_title('z-projection')
            fig.set_xlabel('z (Angstroem)')
            fig.set_ylabel('electron density')
            fig.plot(z,graph)
        else:
            if view == 0:
                x = cell[0]
                ymin = 0
                y = cell[1]
                if plane == None:
                    image = Num.sum(rho, axis = 2)
                    tit = 'xy-projection: sum'
                    Max = Num.max(image)
                    Min = Num.min(image)
                else:
                    rho_rolled = Num.rollaxis(rho,2,0)
                    image = rho_rolled[plane]
                    tit = 'xy-projection: slice '+str(plane)
                xlab = 'x (Angstroem)'
                ylab = 'y (Angstroem)'
            elif view == 1:
                x = cell[0]
                ymin = zmin
                y = cell[2] + zmin
                if plane == None:
                    image = Num.sum(rho, axis = 1)
                    tit = 'xz-projection: sum'
                    Max = Num.max(image)
                    Min = Num.min(image)
                else:
                    rho_rolled = Num.rollaxis(rho,1,0)
                    image = rho_rolled[plane]
                    tit = 'xz-projection: slice '+str(plane)
                xlab = 'x (Angstroem)'
                ylab = 'z (Angstroem)'
            elif view == 2:
                x = cell[1]
                ymin = zmin
                y = cell[2] + zmin
                if plane == None:
                    image = Num.sum(rho, axis = 0)
                    tit = 'yz-projection: sum'
                    Max = Num.max(image)
                    Min = Num.min(image)
                else:  
                    image = rho[plane]
                    tit = 'yz-projection: slice '+str(plane)
                xlab = 'y (Angstroem)'
                ylab = 'z (Angstroem)'
            image = Num.transpose(image)
            fig.set_xlabel(xlab)
            fig.set_ylabel(ylab)
            fig.set_title(tit)
            fig.imshow(image, interpolation='bilinear', cmap=matplotlib.cm.hot,\
                       origin='lower', vmin = Min, vmax = Max,\
                       extent = [0,x,ymin,y] )
        fig.figure.canvas.draw()
################################################################################        
def createfourierframe(parent, rho,cell,zmin):
    frame = wxFourierFrame(parent,rho,cell,zmin,\
                           u"\u03c0-surf Fourier Frame", (600,600))
    return frame
