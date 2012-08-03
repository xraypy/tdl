'''
Specfile Integrator
Author: Craig Biwer (cbiwer@uchicago.edu)
Last modified: 7.12.2012
'''

import os
import math
import sys
import time
import copy

import wx
import matplotlib
from matplotlib import pyplot
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
from matplotlib.widgets import Cursor
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from wx.lib.splitter import MultiSplitterWindow
from wx.tools.Editra.src.eclib import pstatbar

from tdl.modules.ana_upgrade import file_locker
from tdl.modules.ana_upgrade import image_data
from tdl.modules.ana_upgrade import ctr_data
from tdl.modules.ana_upgrade import hdf_data
from tdl.modules.geom import gonio_psic
from pds.pcgui import wxHDFToTree
from pds.pcgui.wxUtil import wxUtil
from pds.shellutil import mod_import

import wx.lib.inspection

# This class is the interface for integrating data
class Integrator(wx.Frame, wx.Notebook, wxUtil):

        def __init__(self, *args, **kwargs):
            self.newFile = True
            self.firstOpen = True
            self.hdfRoot = None
            self.hdfObject = None
            self.hdfTreeObject = None
            self.customTreeObject = None
            
            # Make the window
            wx.Frame.__init__(self, args[0], -1, title="HDF Integrator",
                              size=(1024, 664))
            
            self.menuBar = wx.MenuBar()
            self.fileMenu = wx.Menu()
            self.loadHDF = self.fileMenu.Append(-1, 'Load Project File...')
            self.fileMenu.AppendSeparator()
            self.saveAttributes = self.fileMenu.Append(-1, 'Save attributes...')
            self.fileMenu.AppendSeparator()
            self.saveHKL = self.fileMenu.Append(-1, 'Save HKLFFerr...',
                                                'Write H, K, L, F, and Ferr\n' +
                                                'values to a file')
            self.saveHKLE = self.fileMenu.Append(-1, 'Save HKLEFFerr...',
                                                 'Write H, K, L, E, F, and\n' + 
                                                 'Ferr values to a file')
            self.menuBar.Append(self.fileMenu, 'File')
            
            self.editMenu = wx.Menu()
            self.copyParams = self.editMenu.Append(-1, 'Copy parameters...')
            self.menuBar.Append(self.editMenu, 'Edit')
            
            self.viewMenu = wx.Menu()
            self.viewAreaCorrection = \
                    self.viewMenu.Append(-1, 'View Area Correction')
            self.menuBar.Append(self.viewMenu, 'View')
            
            self.SetMenuBar(self.menuBar)
            
            # Make the status bar
            self.statusBar = pstatbar.ProgressStatusBar(self)
            self.statusBar.SetFieldsCount(5)
            self.statusBar.SetStatusWidths([238, -1, 238, 56, 91])
            
            self.statusSizer = wx.BoxSizer(wx.VERTICAL)
            
            self.integrateCancel = wx.Button(self.statusBar, label='Cancel')
            self.integrateCancel.Bind(wx.EVT_BUTTON, self.integrateStop)
            self.integrateCancel.Hide()
            
            # Window layout
            self.splitWindow = wx.lib.splitter.MultiSplitterWindow(self)
            
            self.treeBook = wx.Notebook(self.splitWindow,
                                        style=wx.SUNKEN_BORDER)
            self.treePanel = wx.Panel(self.treeBook)
            self.treeInfo = wx.Panel(self.treeBook)
            self.treeBook.AddPage(self.treePanel, "All")
            self.treeBook.AddPage(self.treeInfo, "Info")
            
            self.dataWindow = wx.SplitterWindow(self.splitWindow)
            self.graphPanel = wx.Panel(self.dataWindow, style=wx.SUNKEN_BORDER)
            self.rodPanel = wx.Panel(self.dataWindow, style=wx.SUNKEN_BORDER)
            self.infoPanel = wx.Panel(self.splitWindow, style=wx.SUNKEN_BORDER)
            self.splitWindow.SetOrientation(wx.HORIZONTAL)
            self.splitWindow.AppendWindow(self.treeBook, 180)
            self.splitWindow.AppendWindow(self.dataWindow, 500)
            self.splitWindow.AppendWindow(self.infoPanel)
            self.splitWindow.SetMinimumPaneSize(32)
            self.dataWindow.SetMinimumPaneSize(32)
            self.dataWindow.SplitHorizontally(self.graphPanel,
                                              self.rodPanel, 391)
            
            self.paramBook = wx.Notebook(self.infoPanel, style=wx.SUNKEN_BORDER)
            self.basicPage = wx.Panel(self.paramBook)
            self.morePage = wx.Panel(self.paramBook)
            self.paramBook.AddPage(self.basicPage, "Basic")
            self.paramBook.AddPage(self.morePage, "More")
            
            # Initialize the data tree structure
            self.hdfTree = wxHDFToTree.myTreeCtrl(self.treePanel)
            
            # Make the custom selector window
            self.customSelection = customSelector(self)
            
            # Populate the sizers
            self.fig4 = Figure(figsize=(1, 1))
            self.canvas4 = FigCanvas(self.graphPanel, -1, self.fig4)
            
            self.rodFig = Figure(figsize=(1, 1))
            self.rodCanvas = FigCanvas(self.rodPanel, -1, self.rodFig)
            self.rodCanvas.mpl_connect('motion_notify_event',
                                       self.updateCursorStatus)
            self.rodCanvas.mpl_connect('button_release_event',
                                       self.clickToSelect)
            
            self.treeSizer = wx.BoxSizer(wx.VERTICAL)
            self.treeSizer.Add(self.hdfTree, proportion=1,
                               flag=wx.EXPAND | wx.ALL, border=4)
            
            self.graphControlSizer = wx.BoxSizer(wx.HORIZONTAL)
            
            # Bad Point checkbox
            self.badPointLbl = wx.StaticText(self.graphPanel,
                                             label='Bad Point: ')
            self.badPointToggle = wx.CheckBox(self.graphPanel)
            
            # Image max field
            self.imageMaxLbl = wx.StaticText(self.graphPanel,
                                             label='Image Max: ')
            self.imageMaxField = wx.TextCtrl(self.graphPanel, size=(56, -1),
                                             style=wx.TE_PROCESS_ENTER)
            self.keepMaxLbl = wx.StaticText(self.graphPanel, label='Freeze: ')
            self.keepMaxToggle = wx.CheckBox(self.graphPanel)
            self.imageMaxValue = wx.StaticText(self.graphPanel,
                                               label='ROI Max: ')
            
            self.graphControlSizer.Add(self.badPointLbl,
                                       flag=wx.ALIGN_CENTER | wx.LEFT,
                                       border=8)
            self.graphControlSizer.Add(self.badPointToggle,
                                       flag=wx.ALIGN_CENTER)
            self.graphControlSizer.Add(wx.StaticLine(self.graphPanel,
                                                     style=wx.LI_VERTICAL),
                                       flag=wx.EXPAND | wx.LEFT, border=16)
            self.graphControlSizer.Add(self.imageMaxLbl,
                                       flag=wx.ALIGN_CENTER | wx.LEFT,
                                       border=16)
            self.graphControlSizer.Add(self.imageMaxField, flag=wx.ALIGN_CENTER)
            self.graphControlSizer.Add(self.keepMaxLbl,
                                       flag=wx.ALIGN_CENTER | wx.LEFT, border=8)
            self.graphControlSizer.Add(self.keepMaxToggle, flag=wx.ALIGN_CENTER)
            self.graphControlSizer.Add(self.imageMaxValue,
                                       flag=wx.ALIGN_CENTER | wx.LEFT, border=8)
            #self.graphControlSizer.AddStretchSpacer(1)
            
            self.graphSizer = wx.BoxSizer(wx.VERTICAL)
            self.graphSizer.Add(self.canvas4, proportion=1,
                                flag=wx.EXPAND | wx.ALL, border=4)
            self.graphSizer.Add(self.graphControlSizer,
                                flag=wx.ALIGN_CENTER | wx.BOTTOM, border=4)
            
            self.rodSizer = wx.BoxSizer(wx.VERTICAL)
            self.rodSizer.Add(self.rodCanvas, proportion=1,
                              flag=wx.EXPAND | wx.ALL, border=4)
            
            # The top is the point parameters, the middle the scan parameters,
            # and the bottom the point information and integrate button
            self.infoSizer1 = wx.BoxSizer(wx.VERTICAL)
            # This spaces the parameter label and the tab box appropriately
            self.paramSizer1 = wx.BoxSizer(wx.VERTICAL)
            # This spaces the basic tab contents
            self.basicSizer1 = wx.BoxSizer(wx.VERTICAL)
            self.basicSizer2 = wx.BoxSizer(wx.HORIZONTAL)
            self.basicSizer3 = wx.BoxSizer(wx.VERTICAL)
            self.basicSizer4 = wx.BoxSizer(wx.VERTICAL)
            
            # The top label
            self.pointParamLbl = wx.StaticText(self.infoPanel,
                                               label='Point Parameters:')
            
            ####################################################################
            # START POINT BASIC PARAM SIZERS
            ####################################################################
            
            # basicSizer1 (vertical) holds basicSizer2, followed by
            # the 'Apply above to:' sizer
            # basicSizer2 (horizontal) holds basicSizer3 (vertical,
            # the first label and the field from each row) followed by
            # basicSizer4 (vertical, the 'Apply to:' label and
            # buttons from each row)
            
            # The basic panel contents
            # xxSizer1 contains the first label and the input
            # field (from each row; xx is the attribute name)
            # xxSizer2 contains the second label and the apply
            # buttons (from each row; xx is the attribute name)
            self.colNbgrSizer1 = wx.BoxSizer(wx.HORIZONTAL)
            self.colNbgrSizer2 = wx.BoxSizer(wx.HORIZONTAL)
            self.colPowerSizer1 = wx.BoxSizer(wx.HORIZONTAL)
            self.colPowerSizer2 = wx.BoxSizer(wx.HORIZONTAL)
            self.colWidthSizer1 = wx.BoxSizer(wx.HORIZONTAL)
            self.colWidthSizer2 = wx.BoxSizer(wx.HORIZONTAL)
            self.rowNbgrSizer1 = wx.BoxSizer(wx.HORIZONTAL)
            self.rowNbgrSizer2 = wx.BoxSizer(wx.HORIZONTAL)
            self.rowPowerSizer1 = wx.BoxSizer(wx.HORIZONTAL)
            self.rowPowerSizer2 = wx.BoxSizer(wx.HORIZONTAL)
            self.rowWidthSizer1 = wx.BoxSizer(wx.HORIZONTAL)
            self.rowWidthSizer2 = wx.BoxSizer(wx.HORIZONTAL)
            self.flagSizer1 = wx.BoxSizer(wx.HORIZONTAL)
            self.flagSizer2 = wx.BoxSizer(wx.HORIZONTAL)
            self.roiSizer1 = wx.BoxSizer(wx.HORIZONTAL)
            self.roiSizer2 = wx.BoxSizer(wx.HORIZONTAL)
            self.rotateSizer1 = wx.BoxSizer(wx.HORIZONTAL)
            self.rotateSizer2 = wx.BoxSizer(wx.HORIZONTAL)
            self.applySizer = wx.BoxSizer(wx.HORIZONTAL)
            
            # Create all the fields here; this allows for easier tab traversal
            # since the order of focus depends on the order of creation
            self.colNbgrField = wx.TextCtrl(self.basicPage,
                                            style=wx.PROCESS_ENTER)
            self.colPowerField = wx.TextCtrl(self.basicPage,
                                             style=wx.PROCESS_ENTER)
            self.colWidthField = wx.TextCtrl(self.basicPage,
                                             style=wx.PROCESS_ENTER)
            self.rowNbgrField = wx.TextCtrl(self.basicPage,
                                            style=wx.PROCESS_ENTER)
            self.rowPowerField = wx.TextCtrl(self.basicPage,
                                             style=wx.PROCESS_ENTER)
            self.rowWidthField = wx.TextCtrl(self.basicPage,
                                             style=wx.PROCESS_ENTER)
            self.flagField = wx.TextCtrl(self.basicPage, style=wx.PROCESS_ENTER)
            self.roiField = wx.TextCtrl(self.basicPage, style=wx.PROCESS_ENTER)
            self.rotateField = wx.TextCtrl(self.basicPage,
                                           style=wx.PROCESS_ENTER)
            
            # The 'Apply above to:' row (placed here for tab traversal order)
            self.applyLbl1 = wx.StaticText(self.basicPage,
                                           label='Apply above to: ')
            self.applyScan = wx.Button(self.basicPage, label='Scan')
            self.applyCustom = wx.Button(self.basicPage, label='Custom...')
            
            self.applySizer.Add(self.applyLbl1,
                                flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                                border=4)
            self.applySizer.Add(self.applyScan, proportion=1,
                                flag=wx.EXPAND | wx.ALL)
            self.applySizer.Add(self.applyCustom, proportion=1,
                                flag=wx.EXPAND | wx.ALL)
            
            # Number of columns to use for background
            self.colNbgrLbl1 = wx.StaticText(self.basicPage,
                                             label='# bgr col: ')
            self.colNbgrScan = wx.Button(self.basicPage, label='Scan')
            self.colNbgrCustom = wx.Button(self.basicPage, label='Custom')
            self.colNbgrFreeze = wx.CheckBox(self.basicPage)
            
            self.colNbgrSizer1.Add(self.colNbgrLbl1,
                                   flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                                   border=4)
            self.colNbgrSizer1.Add(self.colNbgrField, proportion=1,
                                   flag=wx.EXPAND | wx.ALL)
            self.colNbgrSizer2.Add(self.colNbgrScan, proportion=3,
                                   flag=wx.EXPAND | wx.ALL)
            self.colNbgrSizer2.Add(self.colNbgrCustom, proportion=3,
                                   flag=wx.EXPAND | wx.ALL)
            self.colNbgrSizer2.AddStretchSpacer()
            self.colNbgrSizer2.Add(self.colNbgrFreeze, flag=wx.ALIGN_CENTER)
            self.colNbgrSizer2.AddStretchSpacer()
            
            # Power of column background subtraction
            self.colPowerLbl1 = wx.StaticText(self.basicPage,
                                              label='bgr col power: ')
            self.colPowerScan = wx.Button(self.basicPage, label='Scan')
            self.colPowerCustom = wx.Button(self.basicPage, label='Custom')
            self.colPowerFreeze = wx.CheckBox(self.basicPage)
            
            self.colPowerSizer1.Add(self.colPowerLbl1,
                                    flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                                    border=4)
            self.colPowerSizer1.Add(self.colPowerField, proportion=1,
                                    flag=wx.EXPAND | wx.ALL)
            self.colPowerSizer2.Add(self.colPowerScan, proportion=3,
                                    flag=wx.EXPAND | wx.ALL)
            self.colPowerSizer2.Add(self.colPowerCustom, proportion=3,
                                    flag=wx.EXPAND | wx.ALL)
            self.colPowerSizer2.AddStretchSpacer()
            self.colPowerSizer2.Add(self.colPowerFreeze, flag=wx.ALIGN_CENTER)
            self.colPowerSizer2.AddStretchSpacer()
            
            # Width of column background subtraction peaks
            self.colWidthLbl1 = wx.StaticText(self.basicPage,
                                              label='bgr col width: ')
            self.colWidthScan = wx.Button(self.basicPage, label='Scan')
            self.colWidthCustom = wx.Button(self.basicPage, label='Custom')
            self.colWidthFreeze = wx.CheckBox(self.basicPage)
            
            self.colWidthSizer1.Add(self.colWidthLbl1,
                                    flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                                    border=4)
            self.colWidthSizer1.Add(self.colWidthField, proportion=1,
                                    flag=wx.EXPAND | wx.ALL)
            self.colWidthSizer2.Add(self.colWidthScan, proportion=3,
                                    flag=wx.EXPAND | wx.ALL)
            self.colWidthSizer2.Add(self.colWidthCustom, proportion=3,
                                    flag=wx.EXPAND | wx.ALL)
            self.colWidthSizer2.AddStretchSpacer()
            self.colWidthSizer2.Add(self.colWidthFreeze, flag=wx.ALIGN_CENTER)
            self.colWidthSizer2.AddStretchSpacer()
            
            # Number of rows to use for background
            self.rowNbgrLbl1 = wx.StaticText(self.basicPage,
                                             label='# bgr row: ')
            self.rowNbgrScan = wx.Button(self.basicPage, label='Scan')
            self.rowNbgrCustom = wx.Button(self.basicPage, label='Custom')
            self.rowNbgrFreeze = wx.CheckBox(self.basicPage)
            
            self.rowNbgrSizer1.Add(self.rowNbgrLbl1,
                                   flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                                   border=4)
            self.rowNbgrSizer1.Add(self.rowNbgrField, proportion=1,
                                   flag=wx.EXPAND | wx.ALL)
            self.rowNbgrSizer2.Add(self.rowNbgrScan, proportion=3,
                                   flag=wx.EXPAND | wx.ALL)
            self.rowNbgrSizer2.Add(self.rowNbgrCustom, proportion=3,
                                   flag=wx.EXPAND | wx.ALL)
            self.rowNbgrSizer2.AddStretchSpacer()
            self.rowNbgrSizer2.Add(self.rowNbgrFreeze, flag=wx.ALIGN_CENTER)
            self.rowNbgrSizer2.AddStretchSpacer()
            
            # Power of row background subtraction
            self.rowPowerLbl1 = wx.StaticText(self.basicPage,
                                              label='bgr row power: ')
            self.rowPowerScan = wx.Button(self.basicPage, label='Scan')
            self.rowPowerCustom = wx.Button(self.basicPage, label='Custom')
            self.rowPowerFreeze = wx.CheckBox(self.basicPage)
            
            self.rowPowerSizer1.Add(self.rowPowerLbl1,
                                    flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                                    border=4)
            self.rowPowerSizer1.Add(self.rowPowerField, proportion=1,
                                    flag=wx.EXPAND | wx.ALL)
            self.rowPowerSizer2.Add(self.rowPowerScan, proportion=3,
                                    flag=wx.EXPAND | wx.ALL)
            self.rowPowerSizer2.Add(self.rowPowerCustom, proportion=3,
                                    flag=wx.EXPAND | wx.ALL)
            self.rowPowerSizer2.AddStretchSpacer()
            self.rowPowerSizer2.Add(self.rowPowerFreeze, flag=wx.ALIGN_CENTER)
            self.rowPowerSizer2.AddStretchSpacer()
            
            # Width of row background subtraction peaks
            self.rowWidthLbl1 = wx.StaticText(self.basicPage,
                                              label='bgr row width: ')
            self.rowWidthScan = wx.Button(self.basicPage, label='Scan')
            self.rowWidthCustom = wx.Button(self.basicPage, label='Custom')
            self.rowWidthFreeze = wx.CheckBox(self.basicPage)
            
            self.rowWidthSizer1.Add(self.rowWidthLbl1,
                                    flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                                    border=4)
            self.rowWidthSizer1.Add(self.rowWidthField, proportion=1,
                                    flag=wx.EXPAND | wx.ALL)
            self.rowWidthSizer2.Add(self.rowWidthScan, proportion=3,
                                    flag=wx.EXPAND | wx.ALL)
            self.rowWidthSizer2.Add(self.rowWidthCustom, proportion=3,
                                    flag=wx.EXPAND | wx.ALL)
            self.rowWidthSizer2.AddStretchSpacer()
            self.rowWidthSizer2.Add(self.rowWidthFreeze, flag=wx.ALIGN_CENTER)
            self.rowWidthSizer2.AddStretchSpacer()
            
            # What method to use for background subtraction
            self.flagLbl1 = wx.StaticText(self.basicPage, label='Flag: ')
            self.flagScan = wx.Button(self.basicPage, label='Scan')
            self.flagCustom = wx.Button(self.basicPage, label='Custom')
            self.flagFreeze = wx.CheckBox(self.basicPage)
            
            self.flagSizer1.Add(self.flagLbl1,
                                flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                                border=4)
            self.flagSizer1.Add(self.flagField, proportion=1,
                                flag=wx.EXPAND | wx.ALL)
            self.flagSizer2.Add(self.flagScan, proportion=3,
                                flag=wx.EXPAND | wx.ALL)
            self.flagSizer2.Add(self.flagCustom, proportion=3,
                                flag=wx.EXPAND | wx.ALL)
            self.flagSizer2.AddStretchSpacer()
            self.flagSizer2.Add(self.flagFreeze, flag=wx.ALIGN_CENTER)
            self.flagSizer2.AddStretchSpacer()
            
            # The ROI to use
            self.roiLbl1 = wx.StaticText(self.basicPage, label='ROI: ')
            self.roiScan = wx.Button(self.basicPage, label='Scan')
            self.roiCustom = wx.Button(self.basicPage, label='Custom')
            self.roiFreeze = wx.CheckBox(self.basicPage)
            
            self.roiSizer1.Add(self.roiLbl1,
                               flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                               border=4)
            self.roiSizer1.Add(self.roiField, proportion=2,
                               flag=wx.EXPAND | wx.ALL)
            self.roiSizer2.Add(self.roiScan, proportion=3,
                               flag=wx.EXPAND | wx.ALL)
            self.roiSizer2.Add(self.roiCustom, proportion=3,
                               flag=wx.EXPAND | wx.ALL)
            self.roiSizer2.AddStretchSpacer()
            self.roiSizer2.Add(self.roiFreeze, flag=wx.ALIGN_CENTER)
            self.roiSizer2.AddStretchSpacer()
            
            # How much (if any) to rotate the image
            self.rotateLbl1 = wx.StaticText(self.basicPage,
                                            label='Rotate image: ')
            self.rotateScan = wx.Button(self.basicPage, label='Scan')
            self.rotateCustom = wx.Button(self.basicPage, label='Custom')
            self.rotateFreeze = wx.CheckBox(self.basicPage)
            
            self.rotateSizer1.Add(self.rotateLbl1,
                                  flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                                  border=4)
            self.rotateSizer1.Add(self.rotateField, proportion=1,
                                  flag=wx.EXPAND | wx.ALL)
            self.rotateSizer2.Add(self.rotateScan, proportion=3,
                                  flag=wx.EXPAND | wx.ALL)
            self.rotateSizer2.Add(self.rotateCustom, proportion=3,
                                  flag=wx.EXPAND | wx.ALL)
            self.rotateSizer2.AddStretchSpacer()
            self.rotateSizer2.Add(self.rotateFreeze, flag=wx.ALIGN_CENTER)
            self.rotateSizer2.AddStretchSpacer()
            
            # Add all the label / input field pairs
            self.basicSizer3.Add((1,1), proportion=1, flag=wx.EXPAND | wx.ALL)
            self.basicSizer3.Add(self.colNbgrSizer1, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer3.Add(self.colPowerSizer1, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer3.Add(self.colWidthSizer1, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer3.Add(self.rowNbgrSizer1, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer3.Add(self.rowPowerSizer1, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer3.Add(self.rowWidthSizer1, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer3.Add(self.flagSizer1, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer3.Add(self.roiSizer1, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer3.Add(self.rotateSizer1, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            
            # The sizer for the headings above the scan / custom 
            # buttons and freeze check box
            self.headingSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.headingSizer.AddStretchSpacer(prop=5)
            self.headingSizer.Add(wx.StaticText(self.basicPage,
                                                label='Apply to:'),
                                  flag=wx.ALIGN_CENTER | wx.LEFT, border=4)
            self.headingSizer.AddStretchSpacer(prop=6)
            self.headingSizer.Add(wx.StaticText(self.basicPage,
                                                label='Freeze:'),
                                  flag=wx.ALIGN_CENTER | wx.LEFT, border=4)
            self.headingSizer.AddStretchSpacer(prop=1)
            
            # Add all the label / button pairs
            self.basicSizer4.Add(self.headingSizer, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer4.Add(self.colNbgrSizer2, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer4.Add(self.colPowerSizer2, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer4.Add(self.colWidthSizer2, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer4.Add(self.rowNbgrSizer2, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer4.Add(self.rowPowerSizer2, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer4.Add(self.rowWidthSizer2, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer4.Add(self.flagSizer2, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer4.Add(self.roiSizer2, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer4.Add(self.rotateSizer2, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            
            # Add the label / input field and label / button
            # sizers next to each other
            self.basicSizer2.Add(self.basicSizer3, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer2.Add(self.basicSizer4, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            
            # Add the total parameter sizer and the 'Apply above to:' sizer
            self.basicSizer1.Add(self.basicSizer2, proportion=9,
                                 flag=wx.EXPAND | wx.ALL)
            self.basicSizer1.Add(self.applySizer, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            
            self.basicPage.SetSizerAndFit(self.basicSizer1)
            
            # Add a text field to the 'More' page for notes
            self.histSizer = wx.BoxSizer(wx.VERTICAL)
            self.histLabel = wx.StaticText(self.morePage, label='Notes: ')
            self.histBox = wx.TextCtrl(self.morePage,
                                       style=wx.PROCESS_ENTER | wx.TE_MULTILINE)
            self.histButtonSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.histApplyLbl = wx.StaticText(self.morePage, label='Apply to:')
            self.histScan = wx.Button(self.morePage, label='Scan')
            self.histCustom = wx.Button(self.morePage, label='Custom')
            self.histButtonSizer.Add(self.histApplyLbl,
                                     flag=wx.RIGHT | wx.ALIGN_CENTER, border=4)
            self.histButtonSizer.Add(self.histScan, proportion=1,
                                     flag=wx.EXPAND | wx.RIGHT, border=4)
            self.histButtonSizer.Add(self.histCustom, proportion=1,
                                     flag=wx.EXPAND)
            self.histSizer.Add(self.histLabel, flag=wx.TOP | wx.LEFT, border=4)
            self.histSizer.Add(self.histBox, proportion=1,
                               flag=wx.EXPAND | wx.ALL, border=4)
            self.histSizer.Add(self.histButtonSizer,
                               flag=wx.EXPAND | wx.BOTTOM | wx.LEFT | wx.RIGHT,
                               border=4)
            self.morePage.SetSizerAndFit(self.histSizer)
            
            ####################################################################
            # END POINT BASIC PARAM SIZERS
            ####################################################################
            
            self.paramSizer1.Add(self.pointParamLbl)
            self.paramSizer1.Add(self.paramBook, proportion=1,
                                 flag=wx.EXPAND | wx.LEFT, border=4)
            
            ####################################################################
            # START SCAN PARAM SIZERS
            ####################################################################

            # This spaces the parameter label and the field rows appropriately
            self.corrSizer1 = wx.BoxSizer(wx.VERTICAL)
            self.corrSizer2 = wx.BoxSizer(wx.VERTICAL)
            # Each field row contains a label, an entry field,
            # and an 'Apply to Scan' button
            self.scaleSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.beamSlitSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.detSlitSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.sampleAngleSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.sampleDiameterSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.samplePolygonSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.badMapSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.applyCorrSizer = wx.BoxSizer(wx.HORIZONTAL)
            
            # Create a panel to hold all of the scan-level correction values
            # This is used to allow for tab-traversal
            self.corrPanel = wx.Panel(self.infoPanel,
                                      style=wx.SIMPLE_BORDER | wx.TAB_TRAVERSAL)
            self.corrPanel.SetBackgroundColour(wx.WHITE)
            
            self.scaleField = wx.TextCtrl(self.corrPanel,
                                          style=wx.PROCESS_ENTER)
            self.beamSlitField = wx.TextCtrl(self.corrPanel,
                                             style=wx.PROCESS_ENTER)
            self.detSlitField = wx.TextCtrl(self.corrPanel,
                                            style=wx.PROCESS_ENTER)
            self.sampleAngleField = wx.TextCtrl(self.corrPanel,
                                                style=wx.PROCESS_ENTER)
            self.sampleDiameterField = wx.TextCtrl(self.corrPanel,
                                                   style=wx.PROCESS_ENTER)
            self.samplePolygonField = wx.TextCtrl(self.corrPanel,
                                                  style=wx.PROCESS_ENTER)
            self.badMapField = wx.TextCtrl(self.corrPanel,
                                           style=wx.PROCESS_ENTER)
            
            # The 'Apply to Custom' button
            self.applyCorrLbl = wx.StaticText(self.corrPanel,
                                              label='Apply above to: ')
            self.applyCorrCustom = wx.Button(self.corrPanel, label='Custom...')
            
            self.applyCorrSizer.Add(self.applyCorrLbl,
                                    flag=wx.ALIGN_CENTER | wx.RIGHT, border=4)
            self.applyCorrSizer.Add(self.applyCorrCustom, proportion=1,
                                    flag=wx.EXPAND | wx.ALL)
            
            # How much to scale I
            self.scaleLbl1 = wx.StaticText(self.corrPanel, label='Scale: ')
            
            self.scaleSizer.Add(self.scaleLbl1, flag=wx.ALIGN_CENTER | wx.RIGHT,
                                border=4)
            self.scaleSizer.Add(self.scaleField, proportion=1,
                                flag=wx.EXPAND | wx.ALL)
            
            # The beam slits
            self.beamSlitLbl1 = wx.StaticText(self.corrPanel,
                                              label='Beam slits: ')
            
            self.beamSlitSizer.Add(self.beamSlitLbl1,
                                   flag=wx.ALIGN_CENTER | wx.RIGHT, border=4)
            self.beamSlitSizer.Add(self.beamSlitField, proportion=1,
                                   flag=wx.EXPAND | wx.ALL)
            
            # The detector slits
            self.detSlitLbl1 = wx.StaticText(self.corrPanel,
                                             label='Detector slits: ')
            
            self.detSlitSizer.Add(self.detSlitLbl1,
                                  flag=wx.ALIGN_CENTER | wx.RIGHT, border=4)
            self.detSlitSizer.Add(self.detSlitField, proportion=1,
                                  flag=wx.EXPAND | wx.ALL)
            
            # The sample angles
            self.sampleAngleLbl1 = wx.StaticText(self.corrPanel,
                                                 label='Sample angles: ')
            
            self.sampleAngleSizer.Add(self.sampleAngleLbl1,
                                      flag=wx.ALIGN_CENTER | wx.RIGHT,
                                      border=4)
            self.sampleAngleSizer.Add(self.sampleAngleField, proportion=1,
                                      flag=wx.EXPAND | wx.ALL)
            
            # The sample diameter
            self.sampleDiameterLbl1 = wx.StaticText(self.corrPanel,
                                                    label='Sample diameter: ')
            
            self.sampleDiameterSizer.Add(self.sampleDiameterLbl1,
                                         flag=wx.ALIGN_CENTER | wx.RIGHT,
                                         border=4)
            self.sampleDiameterSizer.Add(self.sampleDiameterField, proportion=1,
                                         flag=wx.EXPAND | wx.ALL)
            
            # The sample polygon
            self.samplePolygonLbl1 = wx.StaticText(self.corrPanel,
                                                   label='Sample polygon: ')
            
            self.samplePolygonSizer.Add(self.samplePolygonLbl1,
                                        flag=wx.ALIGN_CENTER | wx.RIGHT,
                                        border=4)
            self.samplePolygonSizer.Add(self.samplePolygonField, proportion=1,
                                        flag=wx.EXPAND | wx.ALL)
            
            # The bad pixel map
            self.badMapLbl1 = wx.StaticText(self.corrPanel,
                                            label='Bad pixel map: ')
            
            self.badMapSizer.Add(self.badMapLbl1,
                                 flag=wx.ALIGN_CENTER | wx.RIGHT, border=4)
            self.badMapSizer.Add(self.badMapField, proportion=1,
                                 flag=wx.EXPAND | wx.ALL)
            
            # The top label
            self.scanParamLbl = wx.StaticText(self.infoPanel,
                                              label='Scan Parameters:')
            
            # Add each row to the larger sizer
            self.corrSizer2.Add(self.scaleSizer, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            self.corrSizer2.Add(self.beamSlitSizer, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            self.corrSizer2.Add(self.detSlitSizer, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            self.corrSizer2.Add(self.sampleAngleSizer, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            self.corrSizer2.Add(self.sampleDiameterSizer, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            self.corrSizer2.Add(self.samplePolygonSizer, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            self.corrSizer2.Add(self.badMapSizer, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            self.corrSizer2.Add(self.applyCorrSizer, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            
            # Calculate the appropriate sizes
            self.corrPanel.SetSizerAndFit(self.corrSizer2)
            
            # Add the label and the panel to be displayed
            self.corrSizer1.Add(self.scanParamLbl, flag=wx.LEFT, border=4)
            self.corrSizer1.Add(self.corrPanel, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            
            ####################################################################
            # END SCAN PARAM SIZERS
            ####################################################################
            
            ####################################################################
            # START DATA SIZERS
            ####################################################################
            
            # The panel to hold the point information
            self.dataPanel = wx.Panel(self.infoPanel, style=wx.SIMPLE_BORDER)
            self.dataPanel.SetBackgroundColour(wx.WHITE)
            
            # The label for the section
            self.dataLbl = wx.StaticText(self.infoPanel,
                                         label='Point Information:')
            
            # The 'Integrate Scan' and 'Integrate Custom' sizer and buttons
            self.integrateSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.integrateScan = wx.Button(self.infoPanel,
                                           label='Integrate Scan')
            self.integrateCustom = wx.Button(self.infoPanel,
                                             label='Integrate Custom...')
            self.integrateSizer.Add(self.integrateScan, proportion=1,
                                    flag=wx.EXPAND)
            self.integrateSizer.Add(self.integrateCustom, proportion=1,
                                    flag=wx.EXPAND)
            
            # The sizers to hold the entire section and
            # the point information, respectively
            self.dataSizer1 = wx.BoxSizer(wx.VERTICAL)
            self.dataSizer2 = wx.BoxSizer(wx.VERTICAL)
            
            # The sizers for each row of point information
            self.hklSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.iSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.fSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.abSecSizer = wx.BoxSizer(wx.HORIZONTAL)
            
            # The labels for H, K, and L
            self.hLbl = wx.StaticText(self.dataPanel, label='H: ')
            self.kLbl = wx.StaticText(self.dataPanel, label='K: ')
            self.lLbl = wx.StaticText(self.dataPanel, label='L: ')
            
            # Adding the HKL labels to their sizer
            self.hklSizer.Add(self.hLbl, proportion=1,
                              flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                              border=4)
            self.hklSizer.Add(self.kLbl, proportion=1,
                              flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                              border=4)
            self.hklSizer.Add(self.lLbl, proportion=1,
                              flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                              border=4)
            
            # The labels for I, Ierr, and Ibgr
            self.iLbl = wx.StaticText(self.dataPanel, label='I: ')
            self.iErrLbl = wx.StaticText(self.dataPanel, label='Ierr: ')
            self.iBgrLbl = wx.StaticText(self.dataPanel, label='Ibgr: ')
            
            # Adding the I labels to their sizer
            self.iSizer.Add(self.iLbl, proportion=1,
                            flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=4)
            self.iSizer.Add(self.iErrLbl, proportion=1,
                            flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=4)
            self.iSizer.Add(self.iBgrLbl, proportion=1,
                            flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=4)
            
            # The labels for F, Ferr, and Ctot
            self.fLbl = wx.StaticText(self.dataPanel, label='F: ')
            self.fErrLbl = wx.StaticText(self.dataPanel, label='Ferr: ')
            self.ctotLbl = wx.StaticText(self.dataPanel, label='Ctot: ')
            
            # Adding the F labels to their sizer
            self.fSizer.Add(self.fLbl, proportion=1,
                            flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=4)
            self.fSizer.Add(self.fErrLbl, proportion=1,
                            flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=4)
            self.fSizer.Add(self.ctotLbl, proportion=1,
                            flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=4)
            
            # The labels for alpha, beta, and seconds
            self.aLbl = wx.StaticText(self.dataPanel, label='Alpha: ')
            self.bLbl = wx.StaticText(self.dataPanel, label='Beta: ')
            self.secLbl = wx.StaticText(self.dataPanel, label='Seconds: ')
            
            # Adding the alpha, beta, and seconds labels to their sizer
            self.abSecSizer.Add(self.aLbl, proportion=1,
                                flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                                border=4)
            self.abSecSizer.Add(self.bLbl, proportion=1,
                                flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                                border=4)
            self.abSecSizer.Add(self.secLbl, proportion=1,
                                flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT,
                                border=4)
            
            # Adding the row sizers to the sizer for the panel
            self.dataSizer2.Add(self.hklSizer, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            self.dataSizer2.Add(self.iSizer, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            self.dataSizer2.Add(self.fSizer, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            self.dataSizer2.Add(self.abSecSizer, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            
            # Calculate the appropriate sizes
            self.dataPanel.SetSizerAndFit(self.dataSizer2)
            
            # Add the label, point information, and 
            # integrate button to the overall sizer
            self.dataSizer1.Add(self.dataLbl, flag=wx.LEFT, border=4)
            self.dataSizer1.Add(self.dataPanel, proportion=1,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            self.dataSizer1.Add(self.integrateSizer,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            
            ####################################################################
            # END DATA SIZERS
            ####################################################################
            
            self.infoSizer1.Add(self.paramSizer1, proportion=12,
                                flag=wx.EXPAND | wx.LEFT, border=4)
            self.infoSizer1.Add(self.corrSizer1, proportion=9, flag=wx.EXPAND)
            self.infoSizer1.Add(self.dataSizer1, proportion=6, flag=wx.EXPAND)
            
            self.treePanel.SetSizerAndFit(self.treeSizer)
            self.graphPanel.SetSizerAndFit(self.graphSizer)
            self.rodPanel.SetSizerAndFit(self.rodSizer)
            self.infoPanel.SetSizerAndFit(self.infoSizer1)
            
            ####################################################################
            # START TOOLTIPS
            ####################################################################
            
            self.colNbgrLbl1.SetToolTipString(labelTips['colNbgr'])
            self.colPowerLbl1.SetToolTipString(labelTips['colPower'])
            self.colWidthLbl1.SetToolTipString(labelTips['colWidth'])
            self.rowNbgrLbl1.SetToolTipString(labelTips['rowNbgr'])
            self.rowPowerLbl1.SetToolTipString(labelTips['rowPower'])
            self.rowWidthLbl1.SetToolTipString(labelTips['rowWidth'])
            self.flagLbl1.SetToolTipString(labelTips['flag'])
            self.roiLbl1.SetToolTipString(labelTips['roi'])
            self.rotateLbl1.SetToolTipString(labelTips['rotate'])
            
            self.scaleLbl1.SetToolTipString(labelTips['scale'])
            self.beamSlitLbl1.SetToolTipString(labelTips['beamSlit'])
            self.detSlitLbl1.SetToolTipString(labelTips['detSlit'])
            self.sampleAngleLbl1.SetToolTipString(labelTips['sampleAngle'])
            self.sampleDiameterLbl1.SetToolTipString(\
                    labelTips['sampleDiameter'])
            self.samplePolygonLbl1.SetToolTipString(labelTips['samplePolygon'])
            self.badMapLbl1.SetToolTipString(labelTips['badMap'])
            
            ####################################################################
            # END TOOLTIPS
            ####################################################################
            
            ####################################################################
            # START BINDINGS
            ####################################################################
            
            # Window bindings
            self.statusBar.Bind(wx.EVT_SIZE, self.onSize)
            self.Bind(wx.EVT_CLOSE, self.onClose)
            
            # Menu bindings
            self.Bind(wx.EVT_MENU, self.saveHKLEFFerr, self.saveHKLE)
            self.Bind(wx.EVT_MENU, self.saveHKLFFerr, self.saveHKL)
            self.Bind(wx.EVT_MENU, self.saveAttrFile, self.saveAttributes)
            self.Bind(wx.EVT_MENU, self.loadFileDialog, self.loadHDF)
            
            self.Bind(wx.EVT_MENU, self.copyFromTo, self.copyParams)
            
            self.Bind(wx.EVT_MENU, self.showAreaCorrection,
                      self.viewAreaCorrection)
            
            # Escape and delete key bindings
            self.Bind(wx.EVT_CHAR_HOOK, self.escapeKey)
            self.Bind(wx.EVT_CHAR_HOOK, self.deleteItem)
            
            # Freeze ROI key binding
            self.Bind(wx.EVT_CHAR_HOOK, self.freezeROI)
            
            # Bind switching items in the tree to
            # updating the graphs and data fields
            self.Bind(wx.EVT_TREE_SEL_CHANGED, self.newSelected, self.hdfTree)
            
            # Bind the bad point toggle to updating the current selection
            self.badPointToggle.Bind(wx.EVT_CHECKBOX, self.updateBadPoint)
            self.Bind(wx.EVT_CHAR_HOOK, self.updateBadPoint)
            
            # Bind the image max field to updating the current selection
            self.imageMaxField.Bind(wx.EVT_TEXT_ENTER, self.updateImageMax)
            self.imageMaxField.Bind(wx.EVT_KILL_FOCUS, self.updateImageMax)
            
            # Bind losing focus on a field (tab or click away)
            # to updating that field
            self.colNbgrField.Bind(wx.EVT_KILL_FOCUS, self.updateItem)
            self.colPowerField.Bind(wx.EVT_KILL_FOCUS, self.updateItem)
            self.colWidthField.Bind(wx.EVT_KILL_FOCUS, self.updateItem)
            self.rowNbgrField.Bind(wx.EVT_KILL_FOCUS, self.updateItem)
            self.rowPowerField.Bind(wx.EVT_KILL_FOCUS, self.updateItem)
            self.rowWidthField.Bind(wx.EVT_KILL_FOCUS, self.updateItem)
            self.flagField.Bind(wx.EVT_KILL_FOCUS, self.updateItem)
            self.roiField.Bind(wx.EVT_KILL_FOCUS, self.updateItem)
            self.rotateField.Bind(wx.EVT_KILL_FOCUS, self.updateItem)
            
            self.histBox.Bind(wx.EVT_KILL_FOCUS, self.updateItem)
            
            self.scaleField.Bind(wx.EVT_KILL_FOCUS, self.applyToScan)
            self.beamSlitField.Bind(wx.EVT_KILL_FOCUS, self.applyToScan)
            self.detSlitField.Bind(wx.EVT_KILL_FOCUS, self.applyToScan)
            self.sampleAngleField.Bind(wx.EVT_KILL_FOCUS, self.applyToScan)
            self.sampleDiameterField.Bind(wx.EVT_KILL_FOCUS, self.applyToScan)
            self.samplePolygonField.Bind(wx.EVT_KILL_FOCUS, self.applyToScan)
            self.badMapField.Bind(wx.EVT_KILL_FOCUS, self.applyToScan)
            
            # Bind pressing enter in a field to
            # updating that field (focus remains)
            self.colNbgrField.Bind(wx.EVT_TEXT_ENTER, self.updateItem)
            self.colPowerField.Bind(wx.EVT_TEXT_ENTER, self.updateItem)
            self.colWidthField.Bind(wx.EVT_TEXT_ENTER, self.updateItem)
            self.rowNbgrField.Bind(wx.EVT_TEXT_ENTER, self.updateItem)
            self.rowPowerField.Bind(wx.EVT_TEXT_ENTER, self.updateItem)
            self.rowWidthField.Bind(wx.EVT_TEXT_ENTER, self.updateItem)
            self.flagField.Bind(wx.EVT_TEXT_ENTER, self.updateItem)
            self.roiField.Bind(wx.EVT_TEXT_ENTER, self.updateItem)
            self.rotateField.Bind(wx.EVT_TEXT_ENTER, self.updateItem)
            
            self.scaleField.Bind(wx.EVT_TEXT_ENTER, self.applyToScan)
            self.beamSlitField.Bind(wx.EVT_TEXT_ENTER, self.applyToScan)
            self.detSlitField.Bind(wx.EVT_TEXT_ENTER, self.applyToScan)
            self.sampleAngleField.Bind(wx.EVT_TEXT_ENTER, self.applyToScan)
            self.sampleDiameterField.Bind(wx.EVT_TEXT_ENTER, self.applyToScan)
            self.samplePolygonField.Bind(wx.EVT_TEXT_ENTER, self.applyToScan)
            self.badMapField.Bind(wx.EVT_TEXT_ENTER, self.applyToScan)
            
            # Bind the scan buttons to copying the chosen parameter
            # to all the points in a scan
            self.colNbgrScan.Bind(wx.EVT_BUTTON, self.applyToScan)
            self.colPowerScan.Bind(wx.EVT_BUTTON, self.applyToScan)
            self.colWidthScan.Bind(wx.EVT_BUTTON, self.applyToScan)
            self.rowNbgrScan.Bind(wx.EVT_BUTTON, self.applyToScan)
            self.rowPowerScan.Bind(wx.EVT_BUTTON, self.applyToScan)
            self.rowWidthScan.Bind(wx.EVT_BUTTON, self.applyToScan)
            self.flagScan.Bind(wx.EVT_BUTTON, self.applyToScan)
            self.roiScan.Bind(wx.EVT_BUTTON, self.applyToScan)
            self.rotateScan.Bind(wx.EVT_BUTTON, self.applyToScan)
            self.applyScan.Bind(wx.EVT_BUTTON, self.applyToScan)
            
            self.histScan.Bind(wx.EVT_BUTTON, self.applyToScan)
            
            # Bind the custom buttons to copying the chosen
            # parameter(s) to a chosen subset
            self.colNbgrCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            self.colPowerCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            self.colWidthCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            self.rowNbgrCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            self.rowPowerCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            self.rowWidthCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            self.flagCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            self.roiCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            self.rotateCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            self.applyCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            self.applyCorrCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            
            self.histCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            
            # Bind the 'Integrate Scan' button to integrating the
            # currently selected scan
            self.integrateScan.Bind(wx.EVT_BUTTON, self.integrateSelectedScan)
            # Bind the integrate button to integrating a chosen subset
            self.integrateCustom.Bind(wx.EVT_BUTTON, self.applyToCustom)
            
            ####################################################################
            # END BINDINGS
            ####################################################################
            
            # Show the window
            self.statusSizer.Add(self.splitWindow, proportion=664,
                                 flag=wx.EXPAND)
            self.statusSizer.Add(self.statusBar, proportion=31, flag=wx.EXPAND)
            self.SetSizer(self.statusSizer)
            self.CenterOnScreen()
            time3 = time.time()
            self.Show()
            #wx.lib.inspection.InspectionTool().Show()
            #self.Bind(wx.EVT_CLOSE, self.on_close)
            #print 'time 1: ', time1, '(', time1-time1, ')'
            #print 'time 2: ', time2, '(', time2-time1, ')'
            #print 'time 25: ', time25, '(', time25-time1, ')'
            #print 'time 3: ', time3, '(', time3-time1, ')'

        # Open a dialog to choose an HDF file to load
        def loadFileDialog(self, event):
            loadDialog = wx.FileDialog(self, message='Load file...',
                                       defaultDir=os.getcwd(), defaultFile='',
                                       wildcard='HDF files (*.ph5)|*.ph5|'+\
                                                'All file(*.*)|*',
                                       style=wx.OPEN)
            if loadDialog.ShowModal() == wx.ID_OK:
                if not os.path.isfile(loadDialog.GetPath()):
                    print 'Error: File does not exist'
                    return
                try:
                    self.hdfObject.close()
                    self.lockFile.release()
                    print 'Lock released'
                except:
                    pass
                
                # Since hdfObject is now closed, clear all variables
                # to prevent attempted read / writes
                self.newFile = True
                self.firstOpen = True
                self.hdfTree.DeleteAllItems()
                try:
                    self.customSelection.customTree.DeleteAllItems()
                except:
                    pass
                del self.hdfRoot
                del self.hdfObject
                del self.hdfTreeObject
                del self.customTreeObject
                self.hdfRoot = None
                self.hdfObject = None
                self.hdfTreeObject = None
                self.customTreeObject = None
                
                print 'Loading ' + loadDialog.GetPath()
                self.loadFile(loadDialog.GetPath())
            loadDialog.Destroy()
        
        # Load an HDF file into an hdf_data object and parse
        # it into a tree.
        def loadFile(self, fileName):
            try:
                self.hdfObject = hdf_data.HdfDataFile(fileName)
            except file_locker.FileLockException as e:
                print 'Error: ' + str(e)
                return
            except:
                return
            #self.set_data('hdfObject', self.hdfObject)
            self.hdfTreeObject = wxHDFToTree.hdfToTree()
            self.hdfTreeObject.populateTree(self.hdfTree, self.hdfObject)
            self.hdfRoot = self.hdfTree.GetRootItem()
            self.hdfTree.Expand(self.hdfRoot)
            self.hdfTree.SelectItem(self.hdfRoot)
            
            self.customTreeObject = wxHDFToTree.hdfToTree()
            self.customTreeObject.populateTree(self.customSelection.customTree,
                                            self.hdfObject)
            self.customSelection.customRoot = \
                            self.customSelection.customTree.GetRootItem()
        
        # Return focus to the tree when the escape key is pressed
        def escapeKey(self, event):
            event.Skip()
            if event.GetKeyCode() == wx.WXK_ESCAPE:
                self.hdfTree.SetFocus()
        
        # Deletes whatever is selected when the delete
        # key is pressed while the tree has focus
        def deleteItem(self, event):
            event.Skip()
            if event.GetKeyCode() == wx.WXK_DELETE:
                if not self.FindFocus() == self.hdfTree:
                    return
                thisItem = self.hdfTree.GetSelection()
                thisParent = self.hdfTree.GetItemParent(thisItem)
                itemText = \
                        self.hdfTree.GetItemText(thisItem)#.split()[1][:-1]
                confirmation = wx.MessageDialog(self,
                                                message='Are you sure you ' + \
                                                        'want to delete ' + \
                                                        #'scan ' + \
                                                        itemText + '?',
                                                caption='Warning!',
                                                style=wx.YES_NO | wx.NO_DEFAULT)
                if confirmation.ShowModal() == wx.ID_NO:
                    confirmation.Destroy()
                    return
                confirmation.Destroy()
                
                self.hdfTreeObject.deleteItem(self.hdfTree, self.hdfObject,
                                              thisItem)
                self.customSelection.customTree.DeleteAllItems()
                self.customTreeObject.populateTree(\
                                self.customSelection.customTree, self.hdfObject)
                self.customSelection.customRoot = \
                                self.customSelection.customTree.GetRootItem()
        
        # Toggle the 'Freeze ROI' box when the 'F2' key is pressed
        def freezeROI(self, event):
            event.Skip()
            if event.GetKeyCode() == wx.WXK_F2:
                if not self.roiField.GetValue() == '':
                    self.roiFreeze.SetValue(not self.roiFreeze.GetValue())
        
        # Updates the status bar text with the L and F coordinates
        # of the cursor
        def updateCursorStatus(self, event):
            if event.inaxes:
                self.statusBar.SetStatusText('L: %3.2f, F: %3.2f' % 
                                             (event.xdata, event.ydata), 1)
            else:
                self.statusBar.SetStatusText('', 1)
        
        # Selects the point nearest the location of the click
        def clickToSelect(self, event):
            if event.inaxes:
                closestL = lambda a, l:min(l, key=lambda x:abs(x-a))
                ofMe = self.hdfTree.GetSelection()
                myParent = self.hdfTree.GetItemParent(ofMe)
                childrenList = []
                dataLookup = {}
                item, cookie = self.hdfTree.GetFirstChild(myParent)
                while item:
                    iterData = self.hdfTree.GetItemPyData(item)
                    childrenList.append(iterData)
                    dataLookup[iterData] = item
                    item, cookie = self.hdfTree.GetNextChild(myParent, cookie)
                allLs = {}
                if self.hdfObject.get_all('type',
                        [childrenList[0]])[childrenList[0]].startswith('Escan'):
                    allLs = self.hdfObject.get_all('Energy', childrenList)
                elif self.hdfObject.get_all('type',
                        [childrenList[0]])[childrenList[0]].startswith('ascan'):
                    get_this = self.hdfObject.get_all('info',
                            [childrenList[0]])[childrenList[0]].split()[1]
                    allLs = self.hdfObject.get_all(get_this, childrenList)
                else:
                    allLs = self.hdfObject.get_all('L', childrenList)
                allLs = dict([L,child] for child,L in allLs.iteritems())
                thisPoint = closestL(event.xdata, allLs.keys())
                #item, cookie = self.hdfTree.GetFirstChild(myParent)
                #while item and \
                #        self.hdfTree.GetItemPyData(item) != allLs[thisPoint]:
                #    item, cookie = self.hdfTree.GetNextChild(myParent, cookie)
                try:
                    self.hdfTree.SelectItem(dataLookup[allLs[thisPoint]])#item)
                except:
                    print 'Error: can\'t locate point.'
                self.hdfTree.SetFocus()
        
        # Copy parameters from one set of points to another by closest L value
        def copyFromTo(self, event):
            closestL = lambda a, l:min(l, key=lambda x:abs(x-a))
            self.customSelection.customTree.Expand(\
                            self.customSelection.customRoot)
            self.customSelection.CenterOnParent()
            self.customSelection.SetTitle('Copy from...')
            userReply = self.customSelection.ShowModal()
            if userReply == wx.ID_CANCEL:
                self.hdfTree.SetFocus()
                return
            treeSelections = self.customSelection.customTree.GetSelections()
            fromThese = []
            for selected in treeSelections:
                fromThese.extend(\
                        self.customTreeObject.getRelevantChildren(\
                                self.customSelection.customTree,
                                self.hdfObject,
                                selected))
            fromThese = list(set(fromThese))
            fromDict = {}
            #for parentNumber, myNumber in fromThese:
            #    fromL = scanData[parentNumber][myNumber].get('L', 0.)
            #    fromDict[fromL] = (parentNumber, myNumber)
            fromDict = self.hdfObject.get_all('L', fromThese)
            fromDict = dict([L,child] for child,L in fromDict.iteritems())
            possibleLs = fromDict.keys()
            
            self.customSelection.SetTitle('Apply to...')
            userReply = self.customSelection.ShowModal()
            self.hdfTree.SetFocus()
            self.customSelection.SetTitle('Choose points...')
            if userReply == wx.ID_CANCEL:
                return
            treeSelections = self.customSelection.customTree.GetSelections()
            toThese = []
            for selected in treeSelections:
                toThese.extend(\
                        self.customTreeObject.getRelevantChildren(\
                                self.customSelection.customTree,
                                self.hdfObject,
                                selected))
            toThese = list(set(toThese))
            
            for itemData in toThese:
                #print itemData
                toL = self.hdfObject.get_all('L', [itemData])
                toL = toL[itemData]
                copyFrom = closestL(toL, possibleLs)
                copyData = fromDict[copyFrom]
                self.hdfObject.set_all(('det_0', 'image_changed'),
                                       'True', [itemData])
                self.hdfObject.set_all(('det_0', 'F_changed'),
                                       True, [itemData])
                for key in ('cnbgr', 'cpow', 'cwidth', 'rnbgr', 'rpow',
                            'rwidth', 'bgrflag', 'roi', 'rotangle',
                            'bad_point', 'scale', 'beam_slits', 'det_slits',
                            'sample_angles', 'sample_polygon',
                            'sample_diameter', 'bad_pixel_map'):
                    holding = self.hdfObject.get_all(('det_0', key), [copyData])
                    holding = holding[copyData]
                    self.hdfObject.set_all(('det_0', key), holding, [itemData])
            ofMe = self.hdfTree.GetSelection()
            itemData = self.hdfTree.GetItemPyData(ofMe)
            if itemData is None:
                return
            myParent = self.hdfTree.GetItemParent(ofMe)
            self.updateFourPlot(myParent, itemData)
            self.updateRodPlot(myParent, itemData)
            self.updateFields(itemData)
        
        # Build the appropriate data structure to pass to
        # gonio_psic instead of the unparsed G array
        def buildPsicG(self, itemData):
            psicG = []
            # cell
            psicG.append([self.hdfObject[itemData]['real_a'],
                          self.hdfObject[itemData]['real_b'],
                          self.hdfObject[itemData]['real_c'],
                          self.hdfObject[itemData]['real_alpha'],
                          self.hdfObject[itemData]['real_beta'],
                          self.hdfObject[itemData]['real_gamma'],
                          self.hdfObject[itemData]['lambda']])
            # or0
            psicG.append({'h': [self.hdfObject[itemData]['or0_h'],
                                self.hdfObject[itemData]['or0_k'],
                                self.hdfObject[itemData]['or0_L']],
                          'delta': self.hdfObject[itemData]['or0_del'],
                          'eta': self.hdfObject[itemData]['or0_eta'],
                          'chi': self.hdfObject[itemData]['or0_chi'],
                          'phi': self.hdfObject[itemData]['or0_phi'],
                          'nu': self.hdfObject[itemData]['or0_nu'],
                          'mu': self.hdfObject[itemData]['or0_mu'],
                          'lam': self.hdfObject[itemData]['or0_lambda']})
            # or1
            psicG.append({'h': [self.hdfObject[itemData]['or1_h'],
                                self.hdfObject[itemData]['or1_k'],
                                self.hdfObject[itemData]['or1_L']],
                          'delta': self.hdfObject[itemData]['or1_del'],
                          'eta': self.hdfObject[itemData]['or1_eta'],
                          'chi': self.hdfObject[itemData]['or1_chi'],
                          'phi': self.hdfObject[itemData]['or1_phi'],
                          'nu': self.hdfObject[itemData]['or1_nu'],
                          'mu': self.hdfObject[itemData]['or1_mu'],
                          'lam': self.hdfObject[itemData]['or1_lambda']})
            # n
            psicG.append(list(self.hdfObject[itemData]['haz']))
            return psicG
        
        # Plot the area correction
        def showAreaCorrection(self, event):
            ofMe = self.hdfTree.GetSelection()
            itemData = self.hdfTree.GetItemPyData(ofMe)
            if itemData is None:
                return
            # Because the G array is already parsed, we need to build
            # the appropriate data structures to pass to gonio_psic
            psicG = self.buildPsicG(itemData)
            
            psic = gonio_psic.psic_from_spec(psicG, preparsed=True)
            angles = self.hdfObject[itemData]
            #test = angles.keys()
            #test.sort()
            #print test
            #print angles['Beta']
            #print angles['beta']
            #print angles['Alpha']
            #print angles['alpha']
            psic.set_angles(phi = angles.get('phi', None),
                            chi = angles.get('chi', None),
                            eta = angles.get('eta', None),
                            mu = angles.get('mu', None),
                            nu = angles.get('nu', None),
                            delta = angles.get('del', None))
            #print psic
            beam_slits = eval(self.hdfObject[itemData]['det_0']['beam_slits'])
            det_slits = eval(self.hdfObject[itemData]['det_0']['det_slits'])
            if det_slits == {}:
                det_slits = None
            sample = {'dia': eval(self.hdfObject[itemData]['det_0']\
                                    ['sample_diameter']),
                      'angles': eval(self.hdfObject[itemData]['det_0']\
                                    ['sample_angles']),
                      'polygon': eval(self.hdfObject[itemData]['det_0']\
                                    ['sample_polygon'])}
            cor = ctr_data.CtrCorrectionPsic(gonio=psic, beam_slits=beam_slits,
                                             det_slits=det_slits, sample=sample)
            cor.ctot_stationary(plot=True)
        
        # Update a point when the bad point checkbox is toggled
        def updateBadPoint(self, event):
            event.Skip()
            try:
                if event.GetKeyCode() != wx.WXK_F4:
                    return
                else:
                    self.badPointToggle.SetValue(not \
                                                 self.badPointToggle.GetValue())
            except:
                pass
            try:
                ofMe = self.hdfTree.GetSelection()
                itemData = self.hdfTree.GetItemPyData(ofMe)
            except:
                print 'Error getting tree selection'
                self.badPointToggle.SetValue(False)
                return
            if itemData is None:
                self.badPointToggle.SetValue(False)
                return
            #myText = self.hdfTree.GetItemText(ofMe)
            #myNumber = myText.split()[1]
            myParent = self.hdfTree.GetItemParent(ofMe)
            #parentText = self.hdfTree.GetItemText(myParent)
            #parentNumber = parentText.split()[1][:-1]
            #scanData[parentNumber][myNumber]['badPoint'] = \
            #        self.badPointToggle.GetValue()
            self.hdfObject[itemData]['det_0']['bad_point'] = \
                    str(self.badPointToggle.GetValue())
            self.updateF(itemData)
            self.updateRodPlot(myParent, itemData)
            self.updateLabels(itemData)
        
        # Update a point when the image max is changed
        def updateImageMax(self, event):
            ofMe = self.hdfTree.GetSelection()
            myParent = self.hdfTree.GetItemParent(ofMe)
            itemData = self.hdfTree.GetItemPyData(ofMe)
            if itemData is None:
                self.imageMaxField.Clear()
                return
            newValue = self.imageMaxField.GetValue()
            currentValue = self.hdfObject[itemData]['det_0']['image_max']
            if event.GetEventType() == wx.wxEVT_KILL_FOCUS:
                self.imageMaxField.SetValue(str(currentValue))
                return
            try:
                if str(int(newValue)) == currentValue:
                    return
                elif int(newValue) <= 0 and int(newValue) != -1:
                    self.imageMaxField.SetValue(currentValue)
                else:
                    self.hdfObject[itemData]['det_0']['image_max'] = \
                            str(int(newValue))
                    self.updateFourPlot(myParent, itemData)
                    self.imageMaxValue.SetLabel('ROI Max: ' + \
                            str(self.hdfObject[itemData]['det_0']\
                                    ['real_image_max']))
            except:
                self.imageMaxField.SetValue(str(currentValue))
        
        # When a text field loses focus, this function is
        # called to update the scanData dictionary
        # accordingly. After checking that a valid point is
        # selected, this checks if any change was made. If
        # not, nothing happens; if a change is detected,
        # this attempts to save the new parameter value and
        # update the plots accordingly. Failing this, it
        # simply rewrites the current value of the parameter
        # to the field and no change is made.
        def updateItem(self, event):
            #print 'Updating item'
            event.Skip()
            try:
                ofMe = self.hdfTree.GetSelection()
                itemData = self.hdfTree.GetItemPyData(ofMe)
            except:
                print 'Error getting tree selection'
                self.clearFields()
                return
            if itemData is None:
                self.clearFields()
                return
            possibilities = {self.colNbgrField: ('cnbgr', int),
                                self.colPowerField: ('cpow', float),
                                self.colWidthField: ('cwidth', int),
                                self.rowNbgrField: ('rnbgr', int),
                                self.rowPowerField: ('rpow', float),
                                self.rowWidthField: ('rwidth', int),
                                self.flagField: ('bgrflag', int),
                                self.roiField: ('roi', list),
                                self.rotateField: ('rotangle', float)}
            myText = self.hdfTree.GetItemText(ofMe)
            myNumber = myText.split()[1]
            myParent = self.hdfTree.GetItemParent(ofMe)
            parentText = self.hdfTree.GetItemText(myParent)
            parentNumber = parentText.split()[1][:-1]
            #print possibilities[event.GetEventObject()]
            whatField = event.GetEventObject()
            if whatField == self.histBox:
                self.hdfObject[itemData]['hist'] = str(self.histBox.GetValue())
                return
            updateThis, ofType = possibilities[whatField]
            #print event.GetEventObject()#.GetValue()
            try:
                if eval(self.hdfObject[itemData]['det_0'][updateThis]) == \
                        ofType(eval(whatField.GetValue())):
                    pass
                else:
                    self.hdfObject[itemData]['det_0'][updateThis] = \
                            str(ofType(eval(whatField.GetValue())))
                    whatField.SetValue(self.hdfObject[itemData]\
                                                     ['det_0'][updateThis])
                    self.hdfObject[itemData]['det_0']['image_changed'] = 'True'
                    self.hdfObject[itemData]['det_0']['F_changed'] = True
                    self.updateFourPlot(myParent, itemData)
                    self.updateRodPlot(myParent, itemData)
                    self.updateLabels(itemData)
            except:
                whatField.SetValue(self.hdfObject[itemData]\
                                                 ['det_0'][updateThis])
        
        # Apply a given parameter or all parameters
        # to every point in the current scan
        def applyToScan(self, event):
            #print 'Copying to scan'
            ofMe = self.hdfTree.GetSelection()
            itemData = self.hdfTree.GetItemPyData(ofMe)
            if itemData is None:
                self.clearFields()
                return
            possibilities = \
                    {self.colNbgrScan: ('cnbgr', int, self.colNbgrField),
                     self.colPowerScan: ('cpow', float, self.colPowerField),
                     self.colWidthScan: ('cwidth', int, self.colWidthField),
                     self.rowNbgrScan: ('rnbgr', int, self.rowNbgrField),
                     self.rowPowerScan: ('rpow', float, self.rowPowerField),
                     self.rowWidthScan: ('rwidth', int, self.rowWidthField),
                     self.flagScan: ('bgrflag', int, self.flagField),
                     self.roiScan: ('roi', list, self.roiField),
                     self.rotateScan: ('rotangle', float, self.rotateField),
                     self.scaleField: ('scale', float, self.scaleField),
                     self.beamSlitField: ('beam_slits', dict,
                                          self.beamSlitField),
                     self.detSlitField: ('det_slits', dict, self.detSlitField),
                     self.sampleAngleField: ('sample_angles', dict,
                                             self.sampleAngleField),
                     self.sampleDiameterField: ('sample_diameter', float,
                                                self.sampleDiameterField),
                     self.samplePolygonField: ('sample_polygon', list,
                                               self.samplePolygonField)}
            myParent = self.hdfTree.GetItemParent(ofMe)
            whatButton = event.GetEventObject()
            toChange = []
            if whatButton == self.applyScan:
                toChange = possibilities.values()
            elif whatButton == self.badMapField:
                if whatButton.GetValue() == '':
                    whatButton.SetValue(str(self.hdfObject[itemData]['det_0']\
                                                          ['bad_pixel_map']))
                    return
                elif whatButton.GetValue() == \
                        str(self.hdfObject[itemData]['det_0']['bad_pixel_map']):
                    return
                else:
                    updateThese = []
                    item, cookie = self.hdfTree.GetFirstChild(myParent)
                    while item:
                        updateThese.append(self.hdfTree.GetItemPyData(item))
                        item, cookie = \
                                self.hdfTree.GetNextChild(myParent, cookie)
                    del item
                    allBPM = self.hdfObject.get_all(('det_0', 'bad_pixel_map'),
                                                    updateThese)
                    updateThese = [item for item in updateThese if \
                                str(allBPM[item]) != str(whatButton.GetValue())]
                    self.hdfObject.set_all(('det_0', 'bad_pixel_map'),
                                           str(whatButton.GetValue()),
                                           updateThese)
                    self.hdfObject.set_all(('det_0', 'pixel_map_changed'),
                                           'True', updateThese)
                    self.hdfObject.set_all(('det_0', 'image_changed'),
                                           'True', updateThese)
                    self.hdfObject.set_all(('det_0', 'F_changed'),
                                           True, updateThese)
                    #self.hdfObject[itemData]['det_0']['pixel_map_changed'] = \
                    #                                                    'False'
                    #self.hdfObject.write_point(self.hdfObject[itemData])
                    #self.hdfObject.read_point(itemData)
                    self.updateFourPlot(myParent, itemData)
                    self.updateF(itemData)
                    self.updateRodPlot(myParent, itemData)
                    return
            elif whatButton == self.histScan:
                updateThese = []
                item, cookie = self.hdfTree.GetFirstChild(myParent)
                while item:
                    updateThese.append(self.hdfTree.GetItemPyData(item))
                    item, cookie = self.hdfTree.GetNextChild(myParent, cookie)
                self.hdfObject.set_all('hist', str(self.histBox.GetValue()),
                                       updateThese)
                return
            else:
                toChange = [possibilities[whatButton]]
            if whatButton in [self.scaleField, self.beamSlitField,
                              self.detSlitField, self.sampleAngleField,
                              self.sampleDiameterField,
                              self.samplePolygonField]:
                fOnly = True
            else:
                fOnly = False
            updateThese = []
            item, cookie = self.hdfTree.GetFirstChild(myParent)
            while item:
                updateThese.append(self.hdfTree.GetItemPyData(item))
                item, cookie = self.hdfTree.GetNextChild(myParent, cookie)
            del item
            for updateThis, ofType, whatField in toChange:
                try:
                    allValues = self.hdfObject.get_all(('det_0', updateThis),
                                                       updateThese)
                    justThese = [item for item in updateThese if \
                            str(allValues[item]) != \
                            str(ofType(eval(whatField.GetValue())))]
                    self.hdfObject.set_all(('det_0', updateThis),
                                        str(ofType(eval(whatField.GetValue()))),
                                        justThese)
                    if not fOnly:
                        self.hdfObject.set_all(('det_0', 'image_changed'),
                                               'True', justThese)
                    self.hdfObject.set_all(('det_0', 'F_changed'),
                                           True, justThese)
                except:
                    whatField.SetValue(str(\
                            self.hdfObject[itemData]['det_0'][updateThis]))
            if fOnly:
                self.updateF(itemData)
            self.updateRodPlot(myParent, itemData)
            self.updateLabels(itemData)
        
        # Interrupt the current integration
        def integrateStop(self, event):
            self.integrateContinue = False
        
        # Apply a given parameter or all parameters
        # to a custom selection of scans and points
        def applyToCustom(self, event):
            #print 'Copying to custom'
            ofMe = self.hdfTree.GetSelection()
            itemData = self.hdfTree.GetItemPyData(ofMe)
            if itemData is None and \
                                event.GetEventObject() != self.integrateCustom:
                return
            while itemData is None:
                ofMe = self.hdfTree.GetFirstChild(ofMe)[0]
                itemData = self.hdfTree.GetItemPyData(ofMe)
            myParent = self.hdfTree.GetItemParent(ofMe)
            if self.firstOpen:
                self.firstOpen = False
                self.customSelection.customTree.Expand(\
                                            self.customSelection.customRoot)
                pickMe = self.customTreeObject.reverseLookup[itemData]
                self.customSelection.customTree.SelectItem(pickMe)
                pickMe = self.customSelection.customTree.GetItemParent(pickMe)
                while pickMe != self.customSelection.customRoot:
                    self.customSelection.customTree.Expand(pickMe)
                    pickMe = \
                        self.customSelection.customTree.GetItemParent(pickMe)
            self.customSelection.CenterOnParent()
            userReply = self.customSelection.ShowModal()
            self.hdfTree.SetFocus()
            if userReply == wx.ID_CANCEL:
                return

            treeSelections = self.customSelection.customTree.GetSelections()
            updateThese = []
            for selected in treeSelections:
                updateThese.extend(self.customTreeObject.getRelevantChildren(\
                                        self.customSelection.customTree,
                                        self.hdfObject,
                                        selected))
            updateThese = list(set(updateThese))
            updateThese.sort()
            
            # If integrate custom was pressed, integrate the selected
            # scans and return. Otherwise, skip this section and apply
            # the appropriate attribute(s) to the selected scans.
            if event.GetEventObject() == self.integrateCustom:
                self.integrateContinue = True
                integrateProgress = 0
                self.statusBar.SetRange(len(updateThese))
                self.statusBar.SetProgress(0)
                self.onSize(None)
                self.integrateCancel.Show()
                self.statusBar.Start()
                for iterData in updateThese:
                    iterItem = self.hdfTreeObject.reverseLookup[iterData]
                    iterString = \
                            self.hdfTreeObject.statusString(self.hdfTree,
                                                            self.hdfObject,
                                                            iterItem)
                    self.statusBar.SetStatusText("Integrating " + iterString, 2)
                    if not self.integrateContinue:
                        print 'Integration aborted on ' + iterString.lower()
                        break
                    try:
                        if eval(self.hdfObject[iterData]['det_0']\
                                                        ['image_changed']):
                            self.integratePoint(iterData)
                        if self.hdfObject[iterData]['det_0']['F_changed']:
                            self.updateF(iterData)
                    except:
                        print 'Error reading ' + iterString.lower()
                        raise
                    integrateProgress += 1
                    self.statusBar.SetProgress(integrateProgress)
                    while wx.GetApp().Pending():
                        wx.GetApp().Dispatch()
                        wx.GetApp().Yield(True)
                self.integrateContinue = False
                self.hdfTree.SetFocus()
                if self.hdfTree.GetSelection() == ofMe:
                    #self.updateFourPlot(myParent, itemData)
                    self.updateRodPlot(myParent, itemData)
                self.statusBar.SetStatusText('', 2)
                self.integrateCancel.Hide()
                self.statusBar.Stop()
                return
            elif event.GetEventObject() == self.histCustom:
                self.hdfObject.set_all('hist', str(self.histBox.GetValue()),
                                       updateThese)
                return
            
            possibilities1 = \
                    {self.colNbgrCustom: ('cnbgr', int, self.colNbgrField),
                     self.colPowerCustom: ('cpow', float, self.colPowerField),
                     self.colWidthCustom: ('cwidth', int, self.colWidthField),
                     self.rowNbgrCustom: ('rnbgr', int, self.rowNbgrField),
                     self.rowPowerCustom: ('rpow', float, self.rowPowerField),
                     self.rowWidthCustom: ('rwidth', int, self.rowWidthField),
                     self.flagCustom: ('bgrflag', int, self.flagField),
                     self.roiCustom: ('roi', list, self.roiField),
                     self.rotateCustom: ('rotangle', float, self.rotateField)}
                                
            possibilities2 = \
                    {self.scaleField: ('scale', float, self.scaleField),
                     self.beamSlitField: ('beam_slits', dict, 
                                          self.beamSlitField),
                     self.detSlitField: ('det_slits', dict, self.detSlitField),
                     self.sampleAngleField: ('sample_angles', dict,
                                                    self.sampleAngleField),
                     self.sampleDiameterField: ('sample_diameter', float,
                                                    self.sampleDiameterField),
                     self.samplePolygonField: ('sample_polygon', list,
                                                    self.samplePolygonField),
                     self.badMapField: ('bad_pixel_map', str, self.badMapField)}
            whatButton = event.GetEventObject()
            toChange = []
            fOnly = False
            if whatButton == self.applyCustom:
                toChange = possibilities1.values()
            elif whatButton == self.applyCorrCustom:
                toChange = possibilities2.values()
                fOnly = True
            else:
                toChange = [possibilities1[whatButton]]
            for updateThis, ofType, whatField in toChange:
                try:
                    if updateThis.startswith('bad_pixel_map'):
                        allBPM = self.hdfObject.get_all(\
                                                    ('det_0', 'bad_pixel_map'),
                                                    updateThese)
                        justThese = [item for item in updateThese if \
                                    str(allBPM[item]) != \
                                               str(self.badMapField.GetValue())]
                        self.hdfObject.set_all(('det_0', 'bad_pixel_map'),
                                               str(self.badMapField.GetValue()),
                                               justThese)
                        self.hdfObject.set_all(('det_0', 'pixel_map_changed'),
                                               'True', justThese)
                        self.hdfObject.set_all(('det_0', 'image_changed'),
                                               'True', justThese)
                        self.hdfObject.set_all(('det_0', 'F_changed'),
                                               True, justThese)
                        if justThese:
                            fOnly = False
                    else:
                        allValues = self.hdfObject.get_all(('det_0',
                                                            updateThis),
                                                           updateThese)
                        justThese = [item for item in updateThese if \
                                     str(allValues[item]) != \
                                     str(ofType(eval(whatField.GetValue())))]
                        self.hdfObject.set_all(('det_0', updateThis),
                                       str(ofType(eval(whatField.GetValue()))),
                                       justThese)
                        if not fOnly:
                            self.hdfObject.set_all(('det_0', 'image_changed'),
                                                   'True', justThese)
                        self.hdfObject.set_all(('det_0', 'F_changed'),
                                               True, justThese)
                except:
                    print 'Error updating selection'
                    raise
            if fOnly:
                self.updateF(itemData)
            self.updateRodPlot(myParent, itemData)
            self.updateLabels(itemData)
        
        # Integrate the currently-selected scan
        def integrateSelectedScan(self, event):
            ofMe = self.hdfTree.GetSelection()
            itemData = self.hdfTree.GetItemPyData(ofMe)
            updateThese = []
            if itemData is None:
                firstChild, cookie = self.hdfTree.GetFirstChild(ofMe)
                childData = self.hdfTree.GetItemPyData(firstChild)
                if childData is None:
                    return
                while firstChild:
                    childData = self.hdfTree.GetItemPyData(firstChild)
                    updateThese.append(childData)
                    firstChild, cookie = self.hdfTree.GetNextChild(ofMe, cookie)
            else:
                myParent = self.hdfTree.GetItemParent(ofMe)
                item, cookie = self.hdfTree.GetFirstChild(myParent)
                while item:
                    iterData = self.hdfTree.GetItemPyData(item)
                    updateThese.append(iterData)
                    item, cookie = self.hdfTree.GetNextChild(myParent, cookie)
                
            self.integrateContinue = True
            integrateProgress = 0
            self.statusBar.SetRange(len(updateThese))
            self.statusBar.SetProgress(0)
            self.onSize(None)
            self.integrateCancel.Show()
            self.statusBar.Start()
            for iterData in updateThese:
                iterItem = self.hdfTreeObject.reverseLookup[iterData]
                iterString = \
                        self.hdfTreeObject.statusString(self.hdfTree,
                                                        self.hdfObject,
                                                        iterItem)
                self.statusBar.SetStatusText("Integrating " + iterString, 2)
                if not self.integrateContinue:
                    print 'Integration aborted on ' + iterString.lower()
                    break
                try:
                    if eval(self.hdfObject[iterData]['det_0']\
                                                    ['image_changed']):
                        self.integratePoint(iterData)
                    if self.hdfObject[iterData]['det_0']['F_changed']:
                        self.updateF(iterData)
                except:
                    print 'Error reading ' + iterString.lower()
                    raise
                integrateProgress += 1
                self.statusBar.SetProgress(integrateProgress)
                while wx.GetApp().Pending():
                    wx.GetApp().Dispatch()
                    wx.GetApp().Yield(True)
            self.integrateContinue = False
            self.hdfTree.SetFocus()
            if self.hdfTree.GetSelection() == ofMe and itemData is not None:
                #self.updateFourPlot(myParent, itemData)
                self.updateRodPlot(myParent, itemData)
            self.statusBar.SetStatusText('', 2)
            self.integrateCancel.Hide()
            self.statusBar.Stop()
            return
        
        # Gets the values of the parameters from the
        # hdfObject and writes them to the appropriate
        # entry fields.
        def updateFields(self, itemData):
            #print 'Updating fields'
            self.badPointToggle.SetValue(\
                    eval(self.hdfObject[itemData]['det_0']['bad_point']))
            self.imageMaxField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['image_max']))
            self.imageMaxValue.SetLabel('ROI Max: ' + \
                    str(self.hdfObject[itemData]['det_0']['real_image_max']))
            
            self.colNbgrField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['cnbgr']))
            self.colPowerField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['cpow']))
            self.colWidthField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['cwidth']))
            self.rowNbgrField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['rnbgr']))
            self.rowPowerField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['rpow']))
            self.rowWidthField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['rwidth']))
            self.flagField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['bgrflag']))
            self.roiField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['roi']))
            self.rotateField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['rotangle']))
            
            self.histBox.SetValue(str(self.hdfObject[itemData]['hist']))
            
            self.scaleField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['scale']))
            self.beamSlitField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['beam_slits']))
            self.detSlitField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['det_slits']))
            self.sampleAngleField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['sample_angles']))
            self.sampleDiameterField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['sample_diameter']))
            self.samplePolygonField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['sample_polygon']))
            self.badMapField.SetValue(\
                    str(self.hdfObject[itemData]['det_0']['bad_pixel_map']))
            
            self.updateLabels(itemData)
            
        # Updates the information labels from the hdfObject
        def updateLabels(self, itemData):
            #print 'Updating labels'
            self.hLbl.SetLabel('H: ' + str(self.hdfObject[itemData]['H']))
            self.kLbl.SetLabel('K: ' + str(self.hdfObject[itemData]['K']))
            self.lLbl.SetLabel('L: ' + str(self.hdfObject[itemData]['L']))
            self.iLbl.SetLabel('I: ' + \
                    str(round(self.hdfObject[itemData]['det_0']['I'], 2)))
            self.iErrLbl.SetLabel('Ierr: ' + \
                    str(round(self.hdfObject[itemData]['det_0']['Ierr'], 2)))
            self.iBgrLbl.SetLabel('Ibgr: ' + \
                    str(round(self.hdfObject[itemData]['det_0']['Ibgr'], 2)))
            self.fLbl.SetLabel('F: ' + \
                    str(round(self.hdfObject[itemData]['det_0']['F'], 2)))
            self.fErrLbl.SetLabel('Ferr: ' + \
                    str(round(self.hdfObject[itemData]['det_0']['Ferr'], 2)))
            self.ctotLbl.SetLabel('Ctot: ' + \
                    str(round(self.hdfObject[itemData]['det_0']['ctot'], 2)))
            self.aLbl.SetLabel('Alpha: ' + \
                    str(round(self.hdfObject[itemData]['det_0']['alpha'], 2)))
            self.bLbl.SetLabel('Beta: ' + \
                    str(round(self.hdfObject[itemData]['det_0']['beta'], 2)))
            self.secLbl.SetLabel('Seconds: ' + \
                    str(round(self.hdfObject[itemData]['Seconds'], 2)))
            
        # How to (re)calculate the F value for a point
        # Since the rod plot may need to update regardless
        # of whether or not this function is called
        # (i.e. a new point is selected means the red ball
        # must move, but no new F is necessarily calculated),
        # this function does not call updateRodPlot, leaving
        # that to whichever function initially called updateF
        def updateF(self, itemData):
            #print 'Updating F'
            if eval(self.hdfObject[itemData]['det_0']['bad_point']):
                self.hdfObject[itemData]['det_0']['F'] = 0
                self.hdfObject[itemData]['det_0']['Ferr'] = 0
                self.hdfObject[itemData]['det_0']['F_changed'] = False
                return
            sample_params = \
                    {'dia': \
                            eval(self.hdfObject[itemData]\
                                               ['det_0']['sample_diameter']),
                     'angles': \
                            eval(self.hdfObject[itemData]\
                                               ['det_0']['sample_angles']),
                     'polygon': \
                            eval(self.hdfObject[itemData]
                                               ['det_0']['sample_polygon'])}
            corr_params = \
                    {'scale': eval(self.hdfObject[itemData]['det_0']['scale']),
                     'geom': self.hdfObject[itemData]['geom'],
                     'beam_slits': \
                          eval(self.hdfObject[itemData]['det_0']['beam_slits']),
                     'det_slits': \
                          eval(self.hdfObject[itemData]['det_0']['det_slits']),
                     'sample': sample_params}
            if corr_params['det_slits'] == {}:
                corr_params['det_slits'] = None
            # TPT changed 'numPoints' to 'dims'
            psicG = self.buildPsicG(itemData)
            scan_dict = {'I':[self.hdfObject[itemData]['det_0']['I']],
                         'io':[self.hdfObject[itemData]['io']],
                         'Ierr':[self.hdfObject[itemData]['det_0']['Ierr']],
                         'Ibgr':[self.hdfObject[itemData]['det_0']['Ibgr']],
                         'dims':(1,0),
                         'phi':float(self.hdfObject[itemData].get('phi')),
                         'chi':float(self.hdfObject[itemData].get('chi')),
                         'eta':float(self.hdfObject[itemData].get('eta')),
                         'mu':float(self.hdfObject[itemData].get('mu')),
                         'nu':float(self.hdfObject[itemData].get('nu')),
                         'del':float(self.hdfObject[itemData].get('del')),
                         'G':psicG}
            fDict = ctr_data.image_point_F(scan=scan_dict,
                                            point=0,
                                            corr_params = corr_params,
                                            preparsed=True)
            self.hdfObject[itemData]['det_0']['F'] = fDict['F']
            self.hdfObject[itemData]['det_0']['Ferr'] = fDict['Ferr']
            self.hdfObject[itemData]['det_0']['ctot'] = fDict['ctot']
            self.hdfObject[itemData]['det_0']['alpha'] = fDict['alpha']
            self.hdfObject[itemData]['det_0']['beta'] = fDict['beta']
            self.hdfObject[itemData]['det_0']['F_changed'] = False
            
        # Integrate a point without updating the GUI
        def integratePoint(self, itemData):
            #print 'Integrating scan ' + parentNumber + ' point ' + myNumber
            if eval(self.hdfObject[itemData]['det_0']['pixel_map_changed']):
                self.hdfObject[itemData]['det_0']['pixel_map_changed'] = 'False'
                self.hdfObject.write_point(self.hdfObject[itemData])
                self.hdfObject.read_point(itemData)
            self.hdfObject[itemData]['det_0']['image_changed'] = 'False'
            self.hdfObject[itemData]['det_0']['F_changed'] = True
            imageAna = image_data.ImageAna(
                           self.hdfObject[itemData]['det_0']['corrected_image'],
                           eval(self.hdfObject[itemData]['det_0']['roi']),
                           eval(self.hdfObject[itemData]['det_0']['rotangle']),
                           eval(self.hdfObject[itemData]['det_0']['bgrflag']),
                           eval(self.hdfObject[itemData]['det_0']['cnbgr']),
                           eval(self.hdfObject[itemData]['det_0']['cwidth']),
                           eval(self.hdfObject[itemData]['det_0']['cpow']),
                           eval(self.hdfObject[itemData]['det_0']['ctan']),
                           eval(self.hdfObject[itemData]['det_0']['rnbgr']),
                           eval(self.hdfObject[itemData]['det_0']['rwidth']),
                           eval(self.hdfObject[itemData]['det_0']['rpow']),
                           eval(self.hdfObject[itemData]['det_0']['rtan']),
                           eval(self.hdfObject[itemData]['det_0']['nline']),
                           eval(self.hdfObject[itemData]['det_0']['filter']),
                           eval(self.hdfObject[itemData]['det_0']['compress']),
                           False, # 'plot'
                           None, # 'fig'
                           '', # 'figtitle'
                           im_max=eval(self.hdfObject[itemData]['det_0']\
                                                               ['image_max']))
            # TPT changed getVars to get_vars
            (holding1, # 'clpimg',
              holding2, # 'bgrimg',
              integrated,
              self.hdfObject[itemData]['det_0']['I'],
              self.hdfObject[itemData]['det_0']['Ibgr'],
              self.hdfObject[itemData]['det_0']['Ierr'],
              self.hdfObject[itemData]['det_0']['I_c'],
              self.hdfObject[itemData]['det_0']['I_r'],
              self.hdfObject[itemData]['det_0']['Ibgr_c'],
              self.hdfObject[itemData]['det_0']['Ibgr_r'],
              self.hdfObject[itemData]['det_0']['Ierr_c'],
              self.hdfObject[itemData]['det_0']['Ierr_r']) = imageAna.get_vars()
            self.hdfObject[itemData]['det_0']['integrated'] = str(integrated)
            self.updateF(itemData)
        
        # How to update the L vs I/F plot
        # This update usually occurs last, and so does not call other updates
        def updateRodPlot(self, myParent, itemData):
            #print 'Updating rod plot'
            #numChildren = self.hdfTree.GetChildrenCount(myParent)
            iterList = []
            pendingLList = []
            pendingFList = []
            pendingFerrList = []
            doneLList = []
            doneFList = []
            doneFerrList = []
            item, cookie = self.hdfTree.GetFirstChild(myParent)
            while item:
                iterData = self.hdfTree.GetItemPyData(item)
                iterList.append(iterData)
                item, cookie = self.hdfTree.GetNextChild(myParent, cookie)
            iterImageChanged = \
                    self.hdfObject.get_all(('det_0', 'image_changed'), iterList)
            iterFChanged = self.hdfObject.get_all(('det_0', 'F_changed'),
                                                  iterList)
            if self.hdfObject[itemData]['type'].startswith('Escan'):
                iterLList = self.hdfObject.get_all('Energy', iterList)
            elif self.hdfObject[itemData]['type'].startswith('ascan'):
                get_this = self.hdfObject[itemData]['info'].split()[1]
                iterLList = self.hdfObject.get_all(get_this, iterList)
            else:
                iterLList = self.hdfObject.get_all('L', iterList)
            iterFList = self.hdfObject.get_all(('det_0', 'F'), iterList)
            iterFerrList = self.hdfObject.get_all(('det_0', 'Ferr'), iterList)
            for key in iterImageChanged.keys():
                if not (eval(iterImageChanged[key]) or iterFChanged[key]):
                    doneLList.append(iterLList[key])
                    doneFList.append(iterFList[key])
                    doneFerrList.append(iterFerrList[key])
                else:
                    pendingLList.append(iterLList[key])
                    pendingFList.append(iterFList[key])
                    pendingFerrList.append(iterFerrList[key])
            self.rodFig.clear()
            rodPlot = self.rodFig.add_subplot(111)
            rodPlot.plot(pendingLList, pendingFList, '0.8',
                         marker='.', linestyle='')
            try:
                rodPlot.errorbar(pendingLList, pendingFList, pendingFerrList,
                                 fmt='0.8', linestyle='')
            except:
                pass
            rodPlot.plot(doneLList, doneFList, 'b.')
            try:
                rodPlot.errorbar(doneLList,doneFList,doneFerrList,
                                 fmt ='b', linestyle='')
            except:
                pass
            if self.hdfObject[itemData]['type'].startswith('Escan'):
                rodPlot.plot(self.hdfObject[itemData]['Energy'],
                             self.hdfObject[itemData]['det_0']['F'], 'ro')
            elif self.hdfObject[itemData]['type'].startswith('ascan'):
                rodPlot.plot(self.hdfObject[itemData][get_this],
                             self.hdfObject[itemData]['det_0']['F'], 'ro')
            else:
                rodPlot.plot(self.hdfObject[itemData]['L'],
                             self.hdfObject[itemData]['det_0']['F'], 'ro')
            try:
                if not self.hdfObject[itemData]['type'].startswith('Escan') and\
                   not self.hdfObject[itemData]['type'].startswith('ascan'):
                    rodPlot.semilogy()
            except:
                pass
            doneLList.extend(pendingLList)
            doneFList.extend(pendingFList)
            minL = math.floor(min(doneLList))
            maxL = math.ceil(max(doneLList))
            try:
                minF = min([f for f in doneFList if f > 0])/10.**0.1
            except:
                minF = 0
            maxF = max(doneFList)*(10**0.1)
            if minF == 0 and maxF == 0:
                minF = 0.1
                maxF = 1
            rodPlot.axis([minL, maxL, minF, maxF])
            self.rodCanvas.draw()
        
        # How to update the plot of four graphs
        # (image, roi, row and column sums)
        # Calls updateF if needed, since F depends on the integration
        def updateFourPlot(self, myParent, itemData):
            #print 'Updating four plot'
            
            def updateROIFromClick(eclick, erelease):
                #'eclick and erelease are the press and release events'
                #print " The button you used were: ", \
                #       eclick.button, erelease.button
                x1, y1 = eclick.xdata, eclick.ydata
                x2, y2 = erelease.xdata, erelease.ydata
                thisROI = str(map(int, map(round, [x1, y1, x2, y2])))
                self.hdfObject[itemData]['det_0']['roi'] = thisROI
                # Just to make sure it was written properly:
                thisROI = self.hdfObject[itemData]['det_0']['roi']
                self.roiField.SetValue(thisROI)
                self.hdfObject[itemData]['det_0']['image_changed'] = 'True'
                self.hdfObject[itemData]['det_0']['F_changed'] = True
                self.updateFourPlot(myParent, itemData)
                self.imageMaxValue.SetLabel('ROI Max: ' + \
                       str(self.hdfObject[itemData]['det_0']['real_image_max']))
                self.updateRodPlot(myParent, itemData)
                self.updateLabels(itemData)
                self.hdfTree.SetFocus()

            if eval(self.hdfObject[itemData]['det_0']['pixel_map_changed']):
                self.hdfObject[itemData]['det_0']['pixel_map_changed'] = \
                                                                    'False'
                self.hdfObject.write_point(self.hdfObject[itemData])
                self.hdfObject.read_point(itemData)
            self.fig4.clear()
            if not eval(self.hdfObject[itemData]['det_0']['image_changed']):
                imageAna = \
                        image_data.ImageAna(
                          self.hdfObject[itemData]['det_0']['corrected_image'],
                          eval(self.hdfObject[itemData]['det_0']['roi']),
                          eval(self.hdfObject[itemData]['det_0']['rotangle']),
                          eval(self.hdfObject[itemData]['det_0']['bgrflag']),
                          eval(self.hdfObject[itemData]['det_0']['cnbgr']),
                          eval(self.hdfObject[itemData]['det_0']['cwidth']),
                          eval(self.hdfObject[itemData]['det_0']['cpow']),
                          eval(self.hdfObject[itemData]['det_0']['ctan']),
                          eval(self.hdfObject[itemData]['det_0']['rnbgr']),
                          eval(self.hdfObject[itemData]['det_0']['rwidth']),
                          eval(self.hdfObject[itemData]['det_0']['rpow']),
                          eval(self.hdfObject[itemData]['det_0']['rtan']),
                          eval(self.hdfObject[itemData]['det_0']['nline']),
                          eval(self.hdfObject[itemData]['det_0']['filter']),
                          eval(self.hdfObject[itemData]['det_0']['compress']),
                          False, # 'plot'
                          None, # 'fig'
                          '', # 'figtitle'
                          None, # 'clpimg'
                          None, # 'bgrimg'
                          False, # 'integrated'
                          self.hdfObject[itemData]['det_0']['I'],
                          self.hdfObject[itemData]['det_0']['Ibgr'],
                          self.hdfObject[itemData]['det_0']['Ierr'],
                          self.hdfObject[itemData]['det_0']['I_c'],
                          self.hdfObject[itemData]['det_0']['I_r'],
                          self.hdfObject[itemData]['det_0']['Ibgr_c'],
                          self.hdfObject[itemData]['det_0']['Ibgr_r'],
                          self.hdfObject[itemData]['det_0']['Ierr_c'],
                          self.hdfObject[itemData]['det_0']['Ierr_r'],
                          eval(self.hdfObject[itemData]['det_0']['image_max']))
            else:
                self.hdfObject[itemData]['det_0']['image_changed'] = 'False'
                self.hdfObject[itemData]['det_0']['F_changed'] = True
                imageAna = \
                        image_data.ImageAna(
                           self.hdfObject[itemData]['det_0']['corrected_image'],
                           eval(self.hdfObject[itemData]['det_0']['roi']),
                           eval(self.hdfObject[itemData]['det_0']['rotangle']),
                           eval(self.hdfObject[itemData]['det_0']['bgrflag']),
                           eval(self.hdfObject[itemData]['det_0']['cnbgr']),
                           eval(self.hdfObject[itemData]['det_0']['cwidth']),
                           eval(self.hdfObject[itemData]['det_0']['cpow']),
                           eval(self.hdfObject[itemData]['det_0']['ctan']),
                           eval(self.hdfObject[itemData]['det_0']['rnbgr']),
                           eval(self.hdfObject[itemData]['det_0']['rwidth']),
                           eval(self.hdfObject[itemData]['det_0']['rpow']),
                           eval(self.hdfObject[itemData]['det_0']['rtan']),
                           eval(self.hdfObject[itemData]['det_0']['nline']),
                           eval(self.hdfObject[itemData]['det_0']['filter']),
                           eval(self.hdfObject[itemData]['det_0']['compress']),
                           False, # 'plot'
                           None, # 'fig'
                           '', # 'figtitle'
                           im_max=eval(self.hdfObject[itemData]['det_0']\
                                                               ['image_max']))
                # TPT changed getVars to get_vars
                (holding1, # 'clpimg',
                    holding2, # 'bgrimg',
                    self.hdfObject[itemData]['det_0']['integrated'],
                    self.hdfObject[itemData]['det_0']['I'],
                    self.hdfObject[itemData]['det_0']['Ibgr'],
                    self.hdfObject[itemData]['det_0']['Ierr'],
                    self.hdfObject[itemData]['det_0']['I_c'],
                    self.hdfObject[itemData]['det_0']['I_r'],
                    self.hdfObject[itemData]['det_0']['Ibgr_c'],
                    self.hdfObject[itemData]['det_0']['Ibgr_r'],
                    self.hdfObject[itemData]['det_0']['Ierr_c'],
                    self.hdfObject[itemData]['det_0']['Ierr_r']) = \
                                                            imageAna.get_vars()
            # TPT rename embedPlot to embed_plot
            (im_max,
                colormap, # 'colormap',
                subplot2) = imageAna.embed_plot(self.fig4)
            self.hdfObject[itemData]['det_0']['real_image_max'] = str(im_max)
            if self.hdfObject[itemData]['det_0']['F_changed']:
                self.updateF(itemData)
            toggle_selector.RS = RectangleSelector(subplot2,
                                                   updateROIFromClick,
                                                   drawtype='box',
                                                   useblit=True,
                                                   button=[1],
                                                   minspanx=5,
                                                   minspany=5,
                                                   spancoords='pixels',
                                                   rectprops=\
                                                        {'edgecolor': 'red',
                                                         'alpha': 1,
                                                         'fill': False})
            self.canvas4.draw()
        
        # Clear all input fields of values
        # (e.g. for when there is no scan point selected)
        def clearFields(self):
            #print 'Clearing fields'
            self.badPointToggle.SetValue(False)
            if not self.keepMaxToggle.GetValue():
                self.imageMaxField.Clear()
            
            if not self.colNbgrFreeze.GetValue():
                self.colNbgrField.Clear()
            if not self.colPowerFreeze.GetValue():
                self.colPowerField.Clear()
            if not self.colWidthFreeze.GetValue():
                self.colWidthField.Clear()
            if not self.rowNbgrFreeze.GetValue():
                self.rowNbgrField.Clear()
            if not self.rowPowerFreeze.GetValue():
                self.rowPowerField.Clear()
            if not self.rowWidthFreeze.GetValue():
                self.rowWidthField.Clear()
            if not self.flagFreeze.GetValue():
                self.flagField.Clear()
            if not self.roiFreeze.GetValue():
                self.roiField.Clear()
            if not self.rotateFreeze.GetValue():
                self.rotateField.Clear()
            
            self.histBox.Clear()
            
            self.scaleField.Clear()
            self.beamSlitField.Clear()
            self.detSlitField.Clear()
            self.sampleAngleField.Clear()
            self.sampleDiameterField.Clear()
            self.samplePolygonField.Clear()
            self.badMapField.Clear()
            
            self.hLbl.SetLabel('H: ')
            self.kLbl.SetLabel('K: ')
            self.lLbl.SetLabel('L: ')
            self.iLbl.SetLabel('I: ')
            self.iErrLbl.SetLabel('Ierr: ')
            self.iBgrLbl.SetLabel('Ibgr: ')
            self.fLbl.SetLabel('F: ')
            self.fErrLbl.SetLabel('Ferr: ')
            self.ctotLbl.SetLabel('Ctot: ')
            self.aLbl.SetLabel('Alpha: ')
            self.bLbl.SetLabel('Beta: ')
            self.secLbl.SetLabel('Seconds: ')
        
        # What to do when a new selection is made
        def newSelected(self, event):
            #print 'New selected'
            ofMe = event.GetItem()
            myParent = self.hdfTree.GetItemParent(ofMe)
            itemData = self.hdfTree.GetItemPyData(ofMe)
            #print self.hdfTree.GetItemPyData(ofMe)
            if itemData is None:
                self.fig4.clear()
                self.canvas4.draw()
                self.rodFig.clear()
                self.rodCanvas.draw()
                self.clearFields()
                self.statusBar.SetStatusText(self.hdfTree.GetItemText(ofMe))
            else:
                if eval(self.hdfObject[itemData]['det_0']['pixel_map_changed']):
                    self.hdfObject[itemData]['det_0']['pixel_map_changed'] = \
                                                                        'False'
                    #self.hdfObject.write_point(self.hdfObject[itemData])
                    #self.hdfObject.read_point(itemData)
                if self.keepMaxToggle.GetValue():
                    self.hdfObject[itemData]['det_0']['image_max'] = \
                            str(self.imageMaxField.GetValue())
                possibilities = {self.colNbgrFreeze: (self.colNbgrField,
                                                      'cnbgr', int),
                                    self.colPowerFreeze: (self.colPowerField,
                                                          'cpow', float),
                                    self.colWidthFreeze: (self.colWidthField,
                                                          'cwidth', int),
                                    self.rowNbgrFreeze: (self.rowNbgrField,
                                                         'rnbgr', int),
                                    self.rowPowerFreeze: (self.rowPowerField,
                                                          'rpow', float),
                                    self.rowWidthFreeze: (self.rowWidthField,
                                                          'rwidth', int),
                                    self.flagFreeze: (self.flagField,
                                                      'bgrflag', int),
                                    self.roiFreeze: (self.roiField,
                                                     'roi', list),
                                    self.rotateFreeze: (self.rotateField,
                                                        'rotangle', float)}
                for key in possibilities:
                    if key.GetValue():
                        whatField, updateThis, ofType = possibilities[key]
                        if self.hdfObject[itemData]['det_0'][updateThis] == \
                                       str(ofType(eval(whatField.GetValue()))):
                            pass
                        else:
                            self.hdfObject[itemData]['det_0'][updateThis] = \
                                    str(ofType(eval(whatField.GetValue())))
                            whatField.SetValue(str(self.hdfObject[itemData]\
                                                        ['det_0'][updateThis]))
                            self.hdfObject[itemData]['det_0']\
                                                    ['image_changed'] = 'True'
                            self.hdfObject[itemData]['det_0']\
                                                    ['F_changed'] = True
                self.updateFourPlot(myParent, itemData)
                self.updateRodPlot(myParent, itemData)
                self.updateFields(itemData)
                statusText = self.hdfTreeObject.statusString(self.hdfTree,
                                                             self.hdfObject,
                                                             ofMe)
                self.statusBar.SetStatusText(statusText)

        # Save the current point's attributes into a tab-delimited file
        def saveAttrFile(self, event):
            try:
                ofMe = self.hdfTree.GetSelection()
                itemData = self.hdfTree.GetItemPyData(ofMe)
            except:
                print 'Error getting tree selection'
                return
            if itemData is None:
                return
            
            saveDialog = wx.FileDialog(self, message='Save file...',
                                       defaultDir=os.getcwd(), defaultFile='',
                                       wildcard='txt files (*.txt)|*.txt|'+\
                                                'All files (*.*)|*',
                                       style=wx.SAVE | wx.OVERWRITE_PROMPT)
            if saveDialog.ShowModal() == wx.ID_OK:
                print 'Saving attribute file ' + saveDialog.GetPath()
                try:
                    attributeFile = open(saveDialog.GetPath(), 'w')
                except:
                    print 'Error opening attribute file'
                    saveDialog.Destroy()
                    raise
                try:
                    for key in ['bad_pixel_map', 'beam_slits', 'bgrflag',
                                'cnbgr', 'cpow', 'ctan', 'cwidth',
                                'det_slits', 'rnbgr', 'roi', 'rotangle',
                                'rpow', 'rtan', 'rwidth', 'sample_angles',
                                'sample_diameter', 'sample_polygon', 'scale']:
                        value = self.hdfObject[itemData]['det_0'][key]
                        attributeFile.write(key + '\t' + value + '\n')
                    attributeFile.write('geom\t' + \
                                        self.hdfObject[itemData]['geom'])
                except:
                    print 'Error writing to file'
                    attributeFile.close()
                    saveDialog.Destroy()
                    raise
                attributeFile.close()
            saveDialog.Destroy()
        
        # Write a file to the current directory containing the index, H, K,
        # L, F, and Ferr values for each point
        def saveHKLFFerr(self, event):
            self.customSelection.customTree.Expand(\
                    self.customSelection.customRoot)
            self.customSelection.CenterOnParent()
            userReply = self.customSelection.ShowModal()
            if userReply == wx.ID_CANCEL:
                self.hdfTree.SetFocus()
                return
            treeSelections = self.customSelection.customTree.GetSelections()
            saveThese = []
            for selected in treeSelections:
                saveThese.extend(self.customTreeObject.getRelevantChildren(\
                                        self.customSelection.customTree,
                                        self.hdfObject,
                                        selected))
            saveThese = list(set(saveThese))
            saveThese.sort()
            
            saveDialog = wx.FileDialog(self, message='Save file as...',
                                       defaultDir=os.getcwd(), defaultFile='',
                                       wildcard='lst files (*.lst)|*.lst|' + \
                                                'All files (*.*)|*',
                                       style=wx.SAVE | wx.OVERWRITE_PROMPT)
            if saveDialog.ShowModal() == wx.ID_OK:
                fname = saveDialog.GetPath()
                try:
                    allBadPs = self.hdfObject.get_all(('det_0', 'bad_point'),
                                                      saveThese)
                    allHs = self.hdfObject.get_all('H', saveThese)
                    allKs = self.hdfObject.get_all('K', saveThese)
                    allLs = self.hdfObject.get_all('L', saveThese)
                    allFs = self.hdfObject.get_all(('det_0', 'F'), saveThese)
                    allFerrs = self.hdfObject.get_all(('det_0', 'Ferr'),
                                                      saveThese)
                    f = open(fname, 'w')
                    header = "  #idx %5s %5s %5s %7s %7s\n" % \
                                  ('H','K','L','F','Ferr')
                    f.write(header)
                    for iterData in saveThese:
                        if not eval(allBadPs[iterData]):
                            line = "%6s %3.2f %3.2f %6.3f %6.6g %6.6g\n" % \
                                   (iterData,
                                    allHs[iterData],
                                    allKs[iterData],
                                    allLs[iterData],
                                    allFs[iterData],
                                    allFerrs[iterData])
                            f.write(line)
                    f.close()
                except Exception:
                    oops = wx.MessageDialog(self, 'Error saving file\n' + str(Exception))
                    oops.ShowModal()
                    oops.Destroy()
                    raise
            saveDialog.Destroy()
            
        # Write a file to the current directory containing the index, H, K,
        # L, E, F, and Ferr values for each point
        def saveHKLEFFerr(self, event):
            self.customSelection.customTree.Expand(\
                    self.customSelection.customRoot)
            self.customSelection.CenterOnParent()
            userReply = self.customSelection.ShowModal()
            if userReply == wx.ID_CANCEL:
                self.hdfTree.SetFocus()
                return
            treeSelections = self.customSelection.customTree.GetSelections()
            saveThese = []
            for selected in treeSelections:
                saveThese.extend(self.customTreeObject.getRelevantChildren(\
                                        self.customSelection.customTree,
                                        self.hdfObject,
                                        selected))
            saveThese = list(set(saveThese))
            saveThese.sort()
            
            saveDialog = wx.FileDialog(self, message='Save file as...',
                                       defaultDir=os.getcwd(), defaultFile='',
                                       wildcard='lst files (*.lst)|*.lst|' + \
                                                'All files (*.*)|*',
                                       style=wx.SAVE | wx.OVERWRITE_PROMPT)
            if saveDialog.ShowModal() == wx.ID_OK:
                fname = saveDialog.GetPath()
                try:
                    allBadPs = self.hdfObject.get_all(('det_0', 'bad_point'),
                                                      saveThese)
                    allHs = self.hdfObject.get_all('H', saveThese)
                    allKs = self.hdfObject.get_all('K', saveThese)
                    allLs = self.hdfObject.get_all('L', saveThese)
                    allEs = self.hdfObject.get_all('Energy', saveThese)
                    allFs = self.hdfObject.get_all(('det_0', 'F'), saveThese)
                    allFerrs = self.hdfObject.get_all(('det_0', 'Ferr'),
                                                      saveThese)
                    f = open(fname, 'w')
                    header = "  #idx %5s %5s %5s %5s %7s %7s\n" % \
                                  ('H','K','L','E','F','Ferr')
                    f.write(header)
                    for iterData in saveThese:
                        if not eval(allBadPs[iterData]):
                            line = "%6s %3.2f %3.2f %6.3f %6.5f %6.6g %6.6g\n" % \
                                   (iterData,
                                    allHs[iterData],
                                    allKs[iterData],
                                    allLs[iterData],
                                    allEs[iterData],
                                    allFs[iterData],
                                    allFerrs[iterData])
                            f.write(line)
                    f.close()
                except Exception:
                    oops = wx.MessageDialog(self, 'Error saving file\n' + str(Exception))
                    oops.ShowModal()
                    oops.Destroy()
                    raise
            saveDialog.Destroy()
        
        # Resize the cancel button in the status bar
        def onSize(self, event):
            self.splitWindow.SetSashPosition(1, self.GetSize()[0] - 180 - 344)
            if event != None:
                event.Skip()
            rect = self.statusBar.GetFieldRect(\
                    self.statusBar.GetFieldsCount() - 2)
            self.integrateCancel.SetPosition((rect.x + 2, rect.y + 2))
            self.integrateCancel.SetSize((rect.width - 4, rect.height - 4))
        
        # Close the HDF file, destroy the window and the custom tree window
        def onClose(self, event):
            if self.hdfObject is not None:
                self.hdfObject.close()
            self.customSelection.Destroy()
            del self.customSelection
            self.Destroy()
            del self
            
            
# Tree class identical to a regular wx.TreeCtrl
# except for an overridden sort function
'''class myTreeCtrl(wx.TreeCtrl):
    def __init__(self, *args, **kwargs):
        wx.TreeCtrl.__init__(self, *args, **kwargs)
        
    def OnCompareItems(self, item1, item2):
        return cmp(self.GetItemPyData(item1), self.GetItemPyData(item2))'''

# Opens a new tree for custom selection
class customSelector(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, args[0], -1, title="Choose points...",
                           size=(200, 800),
                           style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        
        self.customTree = wxHDFToTree.myTreeCtrl(self,
                                     style=wx.TR_MULTIPLE | wx.TR_HAS_BUTTONS)
        #self.customRoot = self.customTree.AddRoot(\
        #        os.path.split(self.GetParent().filename)[1])
            
        self.okButton = wx.Button(self, wx.ID_OK, label='Ok')
        self.cancelButton = wx.Button(self, wx.ID_CANCEL, label='Cancel')
        
        self.choosingSizer1 = wx.BoxSizer(wx.VERTICAL)
        self.choosingSizer2 = wx.BoxSizer(wx.HORIZONTAL)
        
        self.choosingSizer2.Add(self.okButton, proportion=1,
                                flag=wx.EXPAND | wx.ALL, border=4)
        self.choosingSizer2.Add(self.cancelButton, proportion=1,
                                flag=wx.EXPAND | wx.ALL, border=4)
        
        self.choosingSizer1.Add(self.customTree, proportion=1,
                                flag=wx.EXPAND | wx.ALL, border=4)
        self.choosingSizer1.Add(self.choosingSizer2, flag=wx.EXPAND)
        
        self.SetSizer(self.choosingSizer1)

# Holds the RectangleSelector used to pick an ROI
def toggle_selector(event):
    print ' Key pressed.'
    '''
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print ' RectangleSelector deactivated.'
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print ' RectangleSelector activated.'
        toggle_selector.RS.set_active(True)
    '''

labelTips = {
    #'I':"""Enter label for intensity. Options:\n%s""",

    #'Inorm':"""Enter label for intensity normalization. Options:\n%s""",

    #'Ierr':"""Enter label for intensity errors. Options:\n%s""",

    #'Ibgr':"""Enter label for intensity background. Options:\n%s""",

    'roi':"""Enter the image roi. The format is [x1,y1,x2,y2] where values
        are in pixels corresponding to two corners of the box. Use the set button to
        select the roi from the image plot zoom value (Fig 1) -> raw scan data plot.""",

    'rotate':"""Enter a rotation angle for the image (counter clockwise degrees)""",

    'flag':"""Enter flag for image background method:
        0 determine row and column backgrounds after summation
        1 determine 2D background using 'c'olumn direction 
        2 determine 2D background using 'r'ow direction
        3 determine 2D background from the average 'r'ow and 'c'olumn directions""",

    'colNbgr':"""Number of background points for linear part of column direction
        (y-direction) bgr fit. If nbgr = 0, no linear fit is included""",

    'colWidth':"""Peak width for the column (y) direction bgr fit.
        The background function should fit features that are in general broader
        than the width value. Estimate cwidth using the peak width in the
        column (y) direction. Note: width = 0 corresponds to no polynomial bgr""",

    'colPower':"""Power of polynomial used in row (x) direction bgr fit.
        pow = 0 results in linear background only (see nbgr).  Larger values of pow
        result in steeper polynomials.""",

    #'bgr col tan':"""Flag for use of tangents in column (y) direction bgr determination.
    #    If 'True' then the local slope is removed when fitting a polynomial to a point.
    #    Helpful for data sitting on a broad sloping background""",

    'rowNbgr':"""Number of background points for linear part of row direction
        (x-direction) bgr fit. If nbgr = 0, no linear fit is included""",

    'rowWidth':"""Peak width for the row (x) direction bgr fit.
        The background function should fit features that are in general broader
        than the width value. Estimate rwidth using the peak width in the
        row (x) direction. Note: width = 0 corresponds to no polynomial bgr""",

    'rowPower':"""Power of polynomial used in row (x) direction bgr fit.
        pow = 0 results in linear background only (see nbgr).  Larger values of pow
        result in steeper polynomials.""",

    #'bgr row tan':"""Flag for use of tangents in row direction (x) bgr determination.
    #    If 'True' then the local slope is removed when fitting a polynomial to a point.
    #    Helpful for data sitting on a broad sloping background""",

    #'bgr nline':"""For 2D backgrounds, number of lines to average for polynomial fits.
    #    This option helps to smear out sharp features in the background""",

    #'bgr filter':"""Flag (True/False) to use a spline filter for 2D background determination""",

    #'bgr compress':"""Compression factor to apply during background fitting.  This will
    #    speed up the background determinations and add some smoothing""",

    'beamSlit':"""Enter the incident beam slit settings: beam_slits = {'horz':.6,'vert':.8}
        horz = beam horz width in mm (total width in lab-z / horz scattering plane)
        vert = beam vert hieght in mm (total width in lab-x / vert scattering plane)
        Note these dimensions should be the values at the sample position (ie measured
        with sample position scans)
        If beam slits are 'None' or {} no area correction will be done""",

    'detSlit':"""Enter the detector slit settings: det_slits = {'horz':1.,'vert':1.}
        horz = det horz width in mm (total width in lab-z / horiz scattering plane)
        vert = det vert hieght in mm (total width in lab-x / vert scattering plane)
        If detector slits are 'None' or {} only a spill-off correction will be computed""",

    'geom':"""Enter goniometer geometry.  Options: psic""",

    'sampleDiameter':"""Enter the diameter (in mm) of a round sample mounted on center
        of the goniometer.  If this is <= 0 then use the sample polygon for computing
        the sample correction (or no sample description if a polygon is also not specified). """,

    'samplePolygon':"""A list of vectors that describe a general polygon sample shape, e.g.:    
            polygon =  [[1.,1.], [.5,1.5], [-1.,1.],[-1.,-1.],[0.,.5],[1.,-1.]]
        Each entry is an [x,y] or [x,y,z] vector pointing to an apex of the sample
        polygon. These vectors should be given in general lab frame coordinates
        (x is lab frame vertical, y is positive along the direction of the beam,
        the origin is the rotation center). The vectors need to be specified at a 
        given set of angles (see 'sample angles').

        If the sample vectors are given at the flat phi and chi values and with
        the correct sample hieght (sample Z set so the sample surface is on the
        rotation center), then the z values of the sample vectors will be zero.
        If 2D vectors are passed we therefore assume these are [x,y,0].  If this
        is the case then make sure:
            angles = {'phi':flatphi,'chi':flatchi,'eta':0.,'mu':0.}

        The easiest way to determine the sample coordinate vectors is to take a picture
        of the sample with a camera mounted so that is looks directly down the omega
        axis and the gonio angles set at the sample flat phi and chi values with
        eta = mu = 0. Then find the sample rotation center and measure the position
        of each corner (in mm) with up being the +x direction, and downstream
        being the +y direction.  """,

    'sampleAngle':"""Sample angles used for the description of the sample polygon.
        Use the format (note if any angles are left out they are assumed zero):
            angles = {'phi':123.5,'chi':0.3,'eta':0.,'mu':0.}
        See 'sample polygon' for more info.""",

    'scale':"""Enter scale factor.  The scale factor multiplies by all the intensity
        values. e.g.if Io ~ 1million cps then using 1e6 as the scale makes the normalized
        intensity close to cps.  ie y = scale*y/norm""",
        
    'badMap':"""The file name of the bad pixel map.
        Enter only the file name, not the entire directory path. The file must be in the 
        current directory."""
}

'''# Load the queue
            for scan in self.inScans:
                self.scanQueue.put(scan)
            
            # Initialize the progress bar window
            self.progressBox = wx.ProgressDialog("Reading scans...",
                                                 "Progress:", len(self.inScans), 
                                                 style=wx.PD_CAN_ABORT)
            self.progressContinue = True
            # Create child threads
            for t in xrange(1):
                myThread = ReaderThread(self, self.scanQueue,
                                        self.progressContinue,
                                        self.imageDirectory, self.filename,
                                        self.fullFilename)
                myThread.setDaemon(True)
                myThread.start()
            
            # Debugging, wait for the scans to be processed
            #print self.scanQueue.qsize()
            time1 = time.time()
            self.scanQueue.join()
            time2 = time.time()
            time25 = time.time()
            #print self.scanQueue.qsize()
            
            # Loop through the processed scans and populate the tree
            self.scanCount = 0
            for scan in self.inScans:
                if not self.progressContinue:
                    print 'Broken on scan', scan[0]
                    break
                number = scan[0]
                type = scan[5]
                subScanRoot = \
                        self.hdfTree.AppendItem(self.hdfRoot,
                                                 'Scan '+number+': %s' % type)
                customScanRoot = \
                        self.customSelection.customTree.AppendItem(\
                                self.customSelection.customRoot,
                                'Scan '+number+': %s' % type)
                self.hdfTree.SetPyData(subScanRoot, int(number))
                self.hdfTree.SetItemPyData(subScanRoot, '4')
                self.customSelection.customTree.SetPyData(customScanRoot,
                                                          int(number))
                self.customSelection.customTree.SetItemPyData(customScanRoot,
                                                              '4')
                points = len(scanData[number].keys())
                for j in xrange(points):
                    point = self.hdfTree.AppendItem(subScanRoot,
                                                     'Point ' + str(j + 1))
                    self.hdfTree.SetItemPyData(point, '5')
                    customPoint = \
                            self.customSelection.customTree.AppendItem(\
                                    customScanRoot, 'Point ' + str(j + 1))
                    self.customSelection.customTree.SetItemPyData(customPoint,
                                                                  '5')
                self.hdfTree.SortChildren(self.hdfRoot)
                self.customSelection.customTree.SortChildren(\
                        self.customSelection.customRoot)
                self.scanCount += 1
                self.progressContinue, holding = \
                        self.progressBox.Update(self.scanCount)
                while wx.GetApp().Pending():
                    wx.GetApp().Dispatch()
                    wx.GetApp().Yield(True)
            
            self.hdfTree.Expand(self.hdfTree.GetRootItem())
            self.hdfTree.SetFocus()
            
            # Make sure the threads stop
            self.progressContinue = False
            self.progressBox.Destroy()'''
'''#Thread class for reading scans from the queue and parsing them
    class ReaderThread(threading.Thread):
    def __init__(self, window, queue, continueBool, imgDir, fName, fullFName):
        #Initializations
        threading.Thread.__init__(self)
        self.scanQueue = queue
        self.progWind = window
        self.progressContinue = continueBool
        self.imageDirectory = imgDir
        self.filename = fName
        self.fullFilename = fullFName
        
    #What to do per scan
    def run(self):
        global scanData
        while self.progressContinue:
            #Get the scan
            scan = self.scanQueue.get()
            
            #Grab the scan number, starting and stopping line numbers from the specfile, and the type
            number = scan[0]
            start = scan[9][6]
            stop = scan[9][7]
            type = scan[5]
            
            #Grab the lock to prevent simultaneous writes to scanData
            queueLock.acquire()
            
            #Determine the image path and initialize the tree data per point
            scanData[number] = {}
            linecache.checkcache(self.fullFilename)
            for i in xrange(0, stop-start):
                imagePath = self.imageDirectory + \
                                    '\\S%0*d\\' % (3, int(number)) + \
                                    self.filename + '_S%0*d_%0*d.tif' % \
                                                    (3, int(number), 3, i)
                scanData[number][str(i+1)] = {'rawData': linecache.getline(self.fullFilename, i+start), 
                                                        'imageFile': imagePath,
                                                        'imageData': imagePath, 
                                                        'H': scan[1],
                                                        'K': scan[2],
                                                        'G': map(float,scan[11].split()),
                                                        'labels': scan[12].split(),
                                                        'type': type,
                                                        'roi':[], 
                                                        'rotangle': 0, 
                                                        'bgrflag': 1, 
                                                        'cnbgr': 5, 
                                                        'cwidth': 15, 
                                                        'cpow': 2., 
                                                        'ctan': False, 
                                                        'rnbgr': 5, 
                                                        'rwidth': 15, 
                                                        'rpow': 0., 
                                                        'rtan': False, 
                                                        'nline': 1, 
                                                        'filter': False, 
                                                        'compress': 1, 
                                                        'plot': False, 
                                                        'fig': None, 
                                                        'figtitle': '', 
                                                        'clpimg': None, 
                                                        'bgrimg': None, 
                                                        'integrated': False, 
                                                        'I': 0, 
                                                        'Ibgr': 0, 
                                                        'Ierr': 0, 
                                                        'I_c': 0, 
                                                        'I_r': 0, 
                                                        'Ibgr_c': 0, 
                                                        'Ibgr_r': 0, 
                                                        'Ierr_c': 0, 
                                                        'Ierr_r': 0, 
                                                        'F': 0, 
                                                        'Ferr': 0, 
                                                        'alpha': 0, 
                                                        'beta': 0, 
                                                        'Seconds': 0, 
                                                        'badPoint': False,
                                                        'imageMax': -1, 
                                                        'scale': 1.e6, 
                                                        'geom': 'psic',
                                                        'beamSlit': {'horz':.1,'vert':1.5},
                                                        'detSlit': {},
                                                        'sampleAngle': {},
                                                        'sampleDiameter': 10.,
                                                        'samplePolygon': [],
                                                        'badPixelMap': None,
                                                        'pixelMapChanged': True,
                                                        'imageChanged':True,
                                                        'fChanged':True
                                                        
                                                }
                #for label in scanData[number][str(i+1)]['labels']:
                for j in xrange(len(scanData[number][str(i+1)]['labels'])):
                    labelName = scanData[number][str(i+1)]['labels'][j]
                    scanData[number][str(i+1)][labelName] = \
                        float(scanData[number][str(i+1)]['rawData'].split()[j])
            
            #Release the lock and signal the processing is complete
            queueLock.release()
            self.scanQueue.task_done()
        return'''
'''#Save non-image data as a .xml file
        def saveSession(self, event):
            saveDialog = wx.FileDialog(self, message='Save file as...', defaultDir=self.directory,
                                        defaultFile='', wildcard='xml files (*.xml)|*.xml|All files (*.*)|*',
                                        style=wx.SAVE | wx.OVERWRITE_PROMPT)
            if saveDialog.ShowModal() == wx.ID_OK:
                fname = saveDialog.GetPath()
                try:
                    xmlRoot = xmlTree.Element('r'+self.hdfTree.GetItemText(self.hdfRoot))
                    for parentKey in scanData:
                        xmlParent = xmlTree.SubElement(xmlRoot, 's'+parentKey)
                        for childKey in scanData[parentKey]:
                            xmlChild = xmlTree.SubElement(xmlParent, 'p'+childKey)
                            for dataKey in scanData[parentKey][childKey]:
                                if dataKey not in ['imageData', 'clpimg', 'bgrimg']:
                                    xmlChild.attrib[dataKey] = str(scanData[parentKey][childKey][dataKey])
                    xmlTree.ElementTree(xmlRoot).write(fname)
                except:
                    print 'oops'
            saveDialog.Destroy()'''
            
'''#Load non-image data from a .xml file
        def loadSession(self, event):
            global scanData
            loadDialog = wx.FileDialog(self, message='Load file...', defaultDir=self.directory,
                                        defaultFile='', wildcard='xml files (*.xml)|*.xml|All files (*.*)|*',
                                        style=wx.OPEN)
            if loadDialog.ShowModal() == wx.ID_OK:
                fname = loadDialog.GetPath()
                try:
                    loadTree = xmlTree.parse(fname)
                    
                    #Load scanData
                    newScanData = {}
                    currentScan = '0'
                    self.hdfTree.DeleteAllItems()
                    self.hdfRoot = self.hdfTree.AddRoot(loadTree.getroot().tag[1:])
                    newSubRoot = self.hdfRoot
                    newPoint = self.hdfRoot
                    #newScanTree = myTreeCtrl(self.treePanel)
                    #newScanRoot = newScanTree.AddRoot(loadTree.getroot().tag[1:])
                    #newSubRoot = newScanRoot
                    #newPoint = newScanRoot
                    newCustomTree = customSelector(self)
                    newCustomTree.customTree.Delete(newCustomTree.customRoot)
                    newCustomTree.customRoot = newCustomTree.customTree.AddRoot(loadTree.getroot().tag[1:])
                    newCustomSub = newCustomTree.customRoot
                    newCustomPoint = newCustomTree.customRoot
                    for scan in loadTree.iter():
                        #print scan.tag
                        if scan.tag[0] == 's':
                            currentScan = scan.tag[1:]
                            newScanData[str(currentScan)] = {}
                            newSubRoot = self.hdfTree.AppendItem(self.hdfRoot, 'Scan ' + scan.tag[1:] + ': ' + list(scan)[0].attrib.get('type', ''))
                            self.hdfTree.SetPyData(newSubRoot, int(scan.tag[1:]))
                            newCustomSub = newCustomTree.customTree.AppendItem(newCustomTree.customRoot, 'Scan ' + scan.tag[1:] + ': ' + list(scan)[0].attrib.get('type', ''))
                            newCustomTree.customTree.SetPyData(newCustomSub, int(scan.tag[1:]))
                        elif scan.tag[0] == 'p':
                            newScanData[str(currentScan)][str(scan.tag[1:])] = scan.attrib
                            newScanData[str(currentScan)][str(scan.tag[1:])]['imageData'] = \
                                newScanData[str(currentScan)][str(scan.tag[1:])]['imageFile']
                            newScanData[str(currentScan)][str(scan.tag[1:])]['clpimg'] = None
                            newScanData[str(currentScan)][str(scan.tag[1:])]['bgrimg'] = None
                            newScanData[str(currentScan)][str(scan.tag[1:])]['imageChanged'] = True
                            newPoint = self.hdfTree.AppendItem(newSubRoot, 'Point ' + scan.tag[1:])
                            self.hdfTree.SetPyData(newPoint, int(scan.tag[1:]))
                            self.hdfTree.SortChildren(newSubRoot)
                            newCustomPoint = newCustomTree.customTree.AppendItem(newCustomSub, 'Point ' + scan.tag[1:])
                            newCustomTree.customTree.SetPyData(newCustomPoint, int(scan.tag[1:]))
                            newCustomTree.customTree.SortChildren(newCustomSub)
                    self.hdfTree.SortChildren(self.hdfRoot)
                    newCustomTree.customTree.SortChildren(newCustomTree.customRoot)
                    
                    for parentKey in newScanData.keys():
                        for childKey in newScanData[parentKey].keys():
                            for dataKey in newScanData[parentKey][childKey].keys():
                                try:
                                    newScanData[parentKey][childKey][dataKey] = \
                                        eval(newScanData[parentKey][childKey][dataKey])
                                except:
                                    pass

                    self.customSelection = newCustomTree
                    scanData = newScanData
                    self.hdfTree.Expand(self.hdfTree.GetRootItem())
                except:
                    print 'Nope'
            loadDialog.Destroy()
            self.customSelection.customTree.Expand(self.customSelection.customRoot)
            self.hdfTree.SelectItem(self.hdfRoot)'''
            
'''# Called by the parent filter window, adds scans to the current session
        def appendScans(self, outScans):
            appendThese = []
            for scan in outScans:
                if scan[0] not in scanData.keys():
                    appendThese.append(scan)
            if appendThese == []:
                return
            #Load the queue
            for scan in appendThese:
                self.scanQueue.put(scan)
            
            #Initialize the progress bar window
            self.progressBox = wx.ProgressDialog("Reading scans...", "Progress:", len(appendThese), \
                                            style=wx.PD_CAN_ABORT)
            self.progressContinue = True
            #Create child threads
            for t in xrange(1):
                myThread = ReaderThread(self, self.scanQueue, self.progressContinue, \
                                                self.imageDirectory, self.filename, self.fullFilename)
                myThread.setDaemon(True)
                myThread.start()
                
            self.scanQueue.join()
            
            self.scanCount = 0
            for scan in appendThese:
                if not self.progressContinue:
                    print 'Broken on scan', scan[0]
                    break
                number = scan[0]
                type = scan[5]
                subScanRoot = self.hdfTree.AppendItem(self.hdfRoot, 'Scan ' + number + ': %s' % type)
                customScanRoot = self.customSelection.customTree.AppendItem(self.customSelection.customRoot, 'Scan ' + number + ': %s' % type)
                self.hdfTree.SetPyData(subScanRoot, int(number))
                self.customSelection.customTree.SetPyData(customScanRoot, int(number))
                points = len(scanData[number].keys())
                for j in xrange(points):
                    self.hdfTree.AppendItem(subScanRoot, 'Point ' + str(j + 1))
                    self.customSelection.customTree.AppendItem(customScanRoot, 'Point ' + str(j + 1))
                self.hdfTree.SortChildren(self.hdfRoot)
                self.customSelection.customTree.SortChildren(self.customSelection.customRoot)
                self.scanCount += 1
                self.progressContinue, holding = self.progressBox.Update(self.scanCount)
                while wx.GetApp().Pending():
                    wx.GetApp().Dispatch()
                    wx.GetApp().Yield(True)'''
        
'''#Send the integrated scans to a plotter window
        def exportToPlotter(self, event):
            import wxPlotter
            plotter = mod_import(wxPlotter)
            plotter = wxPlotter.Plotter(self)'''