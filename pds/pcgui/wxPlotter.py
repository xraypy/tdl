'''
Integrated Scan Plotter
Author: Craig Biwer (cbiwer@uchicago.edu)
Last modified: 11.28.2011
'''

import wx
import copy
import matplotlib
from matplotlib import pyplot
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas

class Plotter(wx.Frame):
    def __init__(self, *args, **kwargs):
        wx.Frame.__init__(self, args[0], -1, title="Specfile Plotter", size=(1024, 664))
        #print self.GetParent().filename
        
        #Create the menu bar and status bar
        self.createMenus()
        self.createStatusBar()
        
        #Whole window
        self.splitWindow = wx.SplitterWindow(self)
        #The left tree panels
        self.treePanes = wx.SplitterWindow(self.splitWindow)
        self.basicPane = wx.Panel(self.treePanes)
        self.transformPane = wx.Panel(self.treePanes)
        #The plot and option panels
        self.plotPanes = wx.SplitterWindow(self.splitWindow)
        self.allPlotPane = wx.Panel(self.plotPanes, style=wx.SUNKEN_BORDER)
        self.optionsPane = wx.Panel(self.plotPanes, style=wx.SUNKEN_BORDER)
        
        #Split the windows
        self.splitWindow.SetMinimumPaneSize(16)
        self.splitWindow.SplitVertically(self.treePanes, self.plotPanes, 180)
        self.treePanes.SetMinimumPaneSize(16)
        self.treePanes.SplitHorizontally(self.basicPane, self.transformPane, 290)
        self.plotPanes.SetMinimumPaneSize(16)
        self.plotPanes.SplitHorizontally(self.allPlotPane, self.optionsPane, 440)
        
        #Tree of scans
        self.basicTree = myTreeCtrl(self.basicPane)
        self.basicRoot = self.basicTree.AddRoot('Project 1')
        #Tree of transforms
        self.transformTree = myTreeCtrl(self.transformPane)
        self.transformRoot = self.transformTree.AddRoot('test too')
        #Plotting canvas
        self.plotFig = Figure(figsize=(1, 1))
        self.plotCanvas = FigCanvas(self.allPlotPane, -1, self.plotFig)
        #Options panel
        self.holding = wx.StaticText(self.optionsPane, label='Holding...')
        
        #Size the tree of scans
        self.basicSizer = wx.BoxSizer(wx.VERTICAL)
        self.basicSizer.Add(self.basicTree, proportion=1, flag=wx.EXPAND | wx.LEFT, border = 4)
        #Size the tree of transforms
        self.transformSizer = wx.BoxSizer(wx.VERTICAL)
        self.transformSizer.Add(self.transformTree, proportion=1, flag=wx.EXPAND | wx.LEFT, border = 4)
        #Size the plotting canvas
        self.plotSizer = wx.BoxSizer(wx.VERTICAL)
        self.plotSizer.Add(self.plotCanvas, proportion=1, flag=wx.EXPAND, border = 4)
        #Size the options panel
        self.optionsSizer = wx.BoxSizer(wx.VERTICAL)
        self.optionsSizer.Add(self.holding, proportion=1, flag=wx.EXPAND)
        
        #Set the sizes
        self.basicPane.SetSizerAndFit(self.basicSizer)
        self.transformPane.SetSizerAndFit(self.transformSizer)
        self.allPlotPane.SetSizerAndFit(self.plotSizer)
        self.optionsPane.SetSizerAndFit(self.optionsSizer)
        
        
        
        #Center the window and show it
        self.CenterOnParent()
        self.Show()
        
        print self.basicPane.GetSize()
        print self.transformPane.GetSize()
        
    #Create the menu bar
    def createMenus(self):
        self.menuBar = wx.MenuBar()
        self.fileMenu = wx.Menu()
        self.saveWork = self.fileMenu.Append(-1, 'Save Session...')
        self.loadWork = self.fileMenu.Append(-1, 'Load Session...')
        self.menuBar.Append(self.fileMenu, 'File')
        
        self.editMenu = wx.Menu()
        self.copyParams = self.editMenu.Append(-1, 'Scan parameters...')
        self.menuBar.Append(self.editMenu, 'Edit')
        
        self.SetMenuBar(self.menuBar)
        
    #Create the status bar
    def createStatusBar(self):
        self.statusBar = self.CreateStatusBar()
        #self.statusBar.SetFieldsCount(5)
        #self.statusBar.SetStatusWidths([238, -1, 238, 56, 91])
    
    #
    #def
    
    
#Tree class identical to a regular wx.TreeCtrl except for an overridden sort function
class myTreeCtrl(wx.TreeCtrl):
    def __init__(self, *args, **kwargs):
        wx.TreeCtrl.__init__(self, *args, **kwargs)
        
    def OnCompareItems(self, item1, item2):
        return cmp(self.GetItemPyData(item1), self.GetItemPyData(item2))
