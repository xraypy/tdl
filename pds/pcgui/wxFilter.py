'''
Filter GUI
Author: Craig Biwer (cbiwer@uchicago.edu)
7/7/2012
'''

import h5py
import math
import os
import time
import wx
import wx.lib.agw.customtreectrl as treemix
import wx.lib.mixins.listctrl as listmix
import wx.lib.agw.ultimatelistctrl as ULC

import tdl.modules.specfile.filtertools as ft
import tdl.modules.specfile.mastertoproject as mtp
from pds.pcgui.wxUtil import wxUtil

POSSIBLE_ATTRIBUTES = ['bad_pixel_map', 'beam_slits', 'bgrflag',
                              'cnbgr', 'colormap', 'cpow', 'ctan', 'cwidth',
                              'det_slits', 'geom', 'rnbgr', 'roi', 'rotangle',
                              'rpow', 'rtan', 'rwidth', 'sample_angles',
                              'sample_diameter', 'sample_polygon', 'scale']

class filterGUI(wx.Frame, wxUtil):
    '''The GUI window for filtering.'''
    def __init__(self, *args, **kwargs):
        wx.Frame.__init__(self, args[0], -1, title='HDF Project File Builder',
                          size=(1440, 890))
        
        # Set up shell
        self.shell = None
        self.init_shell()
        
        # The file being filtered
        self.filterFile = None
        #self.set_data('filter_file', self.filterFile)
        
        # All the scans parsed from the file
        self.scanItems = []
        
        # The attribute dictionary
        self.attrDict = {}
        
        # The possible specfiles
        self.allSpecs = {}
        # The possible HK pairs
        self.allHK = {}
        # The possible L values
        self.LMin, self.LMax = (float('inf'), float('-inf'))
        # The possible scan types
        self.allTypes = {}
        # The various lists to be passed to the filters
        self.possibleYears = []
        # The earliest and latest scan dates (epoch format)
        self.dateMin, self.dateMax = (float('inf'), float('-inf'))
        
        # The total filter results
        self.activeFilters = []
        # The results of the spec filter
        self.specResult = None
        self.specCases = None
        # The results of the HK filter
        self.hkResult = None
        self.hkCases = None
        # The results of the L filter
        self.LResult = None
        self.LCases = None
        # The results of the type filter
        self.typeResult = None
        self.typeCases = None
        # The results of the date filter
        self.dateResult = None
        self.dateCases = None
        
        # The project to be built
        self.projectDict = {}
        
        # Make the window
        self.fullWindow = wx.Panel(self)
        
        ###############################################################
        # In the window:
        # fullSizer holds three vertical sizers:
        #   leftSizer holds the table and associated buttons on the left 
        #   middleSizer holds info in the middle
        #   rightSizer holds the project layout
        ###############################################################
        
        self.fullSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.leftSizer = wx.BoxSizer(wx.VERTICAL)
        self.leftPanel = wx.Panel(self.fullWindow)
        self.middleSizer = wx.BoxSizer(wx.VERTICAL)
        self.middlePanel = wx.Panel(self.fullWindow)
        self.rightSizer = wx.BoxSizer(wx.VERTICAL)
        self.rightPanel = wx.Panel(self.fullWindow)
        
        ###############################################################
        # On the left:
        # filterButtonSizer holds the filtering buttons
        # tableSizer holds the table in which the scans are displayed
        # managementSizer holds the 'Keep Selected' and reset buttons
        ###############################################################
        
        # The filter buttons
        self.filterButtonSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # Filter by specfile / scan number
        self.specButton = wx.Button(self.leftPanel, label='Filter Specfiles ' + 
                                                           'and Scan #s',
                                    size=(-1, 32))
        # Filter by HK pair
        self.hkButton = wx.Button(self.leftPanel, label='Filter HKs',
                                  size=(-1, 32))
        # Filter by L range
        self.LButton = wx.Button(self.leftPanel, label='Filter L Range',
                                 size=(-1, 32))
        # Filter by type
        self.typeButton = wx.Button(self.leftPanel, label='Filter Types',
                                    size=(-1, 32))
        # Filter by date range
        self.dateButton = wx.Button(self.leftPanel, label='Filter Date Range',
                                    size=(-1, 32))
        
        # Populate the sizer
        self.filterButtonSizer.Add(self.specButton, proportion=210,#168,
                                   flag=wx.EXPAND | wx.TOP | wx.BOTTOM,
                                   border=8)
        self.filterButtonSizer.Add(self.hkButton, proportion=79,
                                   flag=wx.EXPAND | wx.TOP | wx.BOTTOM,
                                   border=8)
        self.filterButtonSizer.Add(self.LButton, proportion=121,
                                   flag=wx.EXPAND | wx.TOP | wx.BOTTOM,
                                   border=8)
        self.filterButtonSizer.Add(self.typeButton, proportion=67,
                                   flag=wx.EXPAND | wx.TOP | wx.BOTTOM,
                                   border=8)
        self.filterButtonSizer.AddStretchSpacer(57)
        self.filterButtonSizer.Add(self.dateButton, proportion=176,
                                   flag=wx.EXPAND | wx.TOP | wx.BOTTOM,
                                   border=8)
        
        # The data table
        self.tableSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # The multicolumn list serving as the data table
        self.dataTable = TableDataCtrl(self.leftPanel, style=wx.LC_REPORT |
                                                             wx.LC_HRULES |
                                                             wx.LC_VRULES)
        self.dataTable.InsertColumn(0, heading='Specfile', width=165)
        self.dataTable.InsertColumn(1, heading='#', width=40)
        self.dataTable.InsertColumn(2, heading='H Val', width=40)
        self.dataTable.InsertColumn(3, heading='K Val', width=38)
        self.dataTable.InsertColumn(4, heading='L Start', width=60)
        self.dataTable.InsertColumn(5, heading='L Stop', width=60)
        self.dataTable.InsertColumn(6, heading='Scan Type', width=66)
        self.dataTable.InsertColumn(7, heading='Aborted', width=56)
        self.dataTable.InsertColumn(8, heading='Date')
        
        # Add the table to the sizer so it scales properly
        self.tableSizer.Add(self.dataTable, proportion=1, flag=wx.EXPAND |
                                                               wx.BOTTOM,
                                                               border=2)
        
        # The 'Keep Selected' and reset buttons
        self.managementSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # 'Keep Selected' button
        self.keepButton = wx.Button(self.leftPanel, label='Keep Selected')
        
        # Reset Button
        self.resetButton = wx.Button(self.leftPanel,
                                     label='Reset Filters and Reread Data')
        
        # Populate the sizer
        self.managementSizer.Add(self.keepButton, proportion=2,
                                 flag=wx.EXPAND | wx.BOTTOM, border=2)
        self.managementSizer.AddStretchSpacer(5)
        self.managementSizer.Add(self.resetButton, proportion=3,
                                 flag=wx.EXPAND | wx.BOTTOM, border=2)
        # Add a border line
        self.managementSizer.Add(wx.StaticLine(self.leftPanel, size=(2, 24)),
                                 flag=wx.LEFT, border=16)
        
        # Arrange the left panel
        self.leftSizer.Add(self.filterButtonSizer, proportion=0, flag=wx.EXPAND)
        self.leftSizer.Add(self.tableSizer, proportion=1, flag=wx.EXPAND)
        self.leftSizer.Add(self.managementSizer, proportion=0,
                           flag=wx.EXPAND | wx.TOP | wx.BOTTOM, border=4)
        self.leftPanel.SetSizerAndFit(self.leftSizer)
        
        ###############################################################
        # End of the left layout
        ###############################################################
        
        ###############################################################
        # In the middle:
        # fileSizer holds the filename text
        # newAttributeSizer holds attribute adder and list
        # Everything else is placed directly into middleSizer
        ###############################################################
        
        # The file name texts
        self.fileSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # The static text
        self.fileLabel = wx.StaticText(self.middlePanel, label='File Name: ')
        
        # The load master file button (changes label based on current file)
        self.fileButton = wx.Button(self.middlePanel,
                                    label='Load Master File...')
        
        # Populate the sizer:
        self.fileSizer.Add(self.fileLabel, flag=wx.CENTER)
        self.fileSizer.Add(self.fileButton, proportion=1,
                           flag=wx.EXPAND | wx.RIGHT, border=20)
        
        # The middle contents
        # The 'More Info:' label
        self.moreLabel = wx.StaticText(self.middlePanel, label='More Info:\n')
        
        # The 'More Info' text box
        self.moreBox = wx.TextCtrl(self.middlePanel, style=wx.TE_MULTILINE | 
                                                          wx.TE_READONLY)
        
        # The move buttons
        self.rightOne = wx.Button(self.middlePanel, label='>')
        self.leftOne = wx.Button(self.middlePanel, label='<')
        self.rightAll = wx.Button(self.middlePanel, label='>>')
        self.leftAll = wx.Button(self.middlePanel, label='<<')
        
        # Tooltips for the move buttons
        self.rightOne.SetToolTipString('Move selected to project')
        self.leftOne.SetToolTipString('Remove selected from project')
        self.rightAll.SetToolTipString('Move all to project')
        self.leftAll.SetToolTipString('Remove all from project')
        
        # The new attribute adder line
        self.newAttributeSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # The drop down list
        self.attrSelect = wx.Choice(self.middlePanel, size=(110, -1))
        self.attrSelect.SetItems(POSSIBLE_ATTRIBUTES)
        
        # The text box
        self.attrSpecify = wx.TextCtrl(self.middlePanel)
        
        # The add button
        self.attrAdd = wx.Button(self.middlePanel, label='+', size=(40, -1))
        
        # Arrange the attribute adder line
        self.newAttributeSizer.Add(self.attrSelect)
        self.newAttributeSizer.Add(self.attrSpecify, proportion=1,
                                   flag=wx.EXPAND)
        self.newAttributeSizer.Add(self.attrAdd)
        #self.newAttributeSizer.AddSpacer(19)
        
        # The added attribute list
        self.attrList = ULC.UltimateListCtrl(self.middlePanel,
                                             agwStyle=wx.LC_REPORT |
                                                wx.LC_VRULES |
                                                wx.LC_HRULES |
                                                ULC.ULC_HAS_VARIABLE_ROW_HEIGHT)
        self.attrList.InsertColumn(0, 'Attribute', width=108)
        self.attrList.InsertColumn(1, 'Value', width=40)
        self.attrList.InsertColumn(2, '', width=36)
        self.attrList.SetColumnWidth(1, ULC.ULC_AUTOSIZE_FILL)
        
        # The load / save attribute buttons file line
        self.loadSaveAttributeSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # The load attribute button
        self.loadAttrButton = wx.Button(self.middlePanel,
                                        label='Load Attribute File...')
        # The save attribute button
        self.saveAttrButton = wx.Button(self.middlePanel,
                                        label='Save Attribute File...')
                                  
        # Arrange the load / save attribute buttons line
        self.loadSaveAttributeSizer.Add(self.loadAttrButton, proportion=1,
                                        flag=wx.EXPAND | wx.LEFT,
                                        border=32)
        self.loadSaveAttributeSizer.AddSpacer(16)
        self.loadSaveAttributeSizer.Add(self.saveAttrButton, proportion=1,
                                        flag=wx.EXPAND | wx.RIGHT,
                                        border=32)
        
        # Arrange the middle panel
        self.middleSizer.Add(self.fileSizer, proportion=0,
                            flag=wx.EXPAND | wx.TOP | wx.LEFT, border=20)
        # Add a border line
        self.middleSizer.Add(wx.StaticLine(self.middlePanel, size=(332, 2)),
                            flag=wx.LEFT | wx.TOP | wx.RIGHT |
                            wx.EXPAND, border=16)
        self.middleSizer.Add(self.moreLabel, proportion=0,
                            flag=wx.TOP | wx.LEFT, border=20)
        self.middleSizer.Add(self.moreBox, proportion=1,
                            flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM,
                            border=24)
        self.middleSizer.Add(self.rightOne, proportion=0,
                             flag=wx.CENTER | wx.TOP | wx.BOTTOM, border=4)
        self.middleSizer.Add(self.leftOne, proportion=0,
                             flag=wx.CENTER | wx.TOP | wx.BOTTOM, border=4)
        self.middleSizer.Add(self.rightAll, proportion=0,
                             flag=wx.CENTER | wx.TOP | wx.BOTTOM, border=4)
        self.middleSizer.Add(self.leftAll, proportion=0,
                             flag=wx.CENTER | wx.TOP | wx.BOTTOM, border=4)
        self.middleSizer.Add(self.newAttributeSizer, proportion=0,
                             flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP,
                             border=24)
        self.middleSizer.Add(self.attrList, proportion=1,
                             flag=wx.EXPAND | wx.LEFT | wx.RIGHT,
                             border=24)
        self.middleSizer.AddSpacer(6)
        self.middleSizer.Add(self.loadSaveAttributeSizer, flag=wx.EXPAND)
        self.middleSizer.AddSpacer(6)

        self.middlePanel.SetSizerAndFit(self.middleSizer)
        
        ###############################################################
        # End of the middle layout
        ###############################################################
        
        ###############################################################
        # On the right:
        # projectNameSizer holds the project name label and box
        # nextStepSizer holds the buttons for moving to the integrator
        # Everything else is placed directly into rightSizer
        ###############################################################
        
        # The project name label and box
        self.projectNameSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # The static label
        self.projectNameText = wx.StaticText(self.rightPanel,
                                             label='Project Name: ')
        
        # The text box
        self.projectNameBox = wx.TextCtrl(self.rightPanel)
        
        # Populate the sizer
        self.projectNameSizer.Add(self.projectNameText,
                                  flag=wx.CENTER | wx.LEFT, border=8)
        self.projectNameSizer.Add(self.projectNameBox, proportion=1,
                                  flag=wx.EXPAND | wx.RIGHT, border=8)
        
        # The new project tree box
        self.newProjectTree = myTreeCtrl(self.rightPanel,
                                         style=wx.TR_MULTIPLE |
                                               wx.TR_EXTENDED | 
                                               wx.TR_DEFAULT_STYLE)
        
        # The integrator buttons
        self.nextStepSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # Create a new integrator window
        self.newIntegrator = wx.Button(self.rightPanel,
                                       label='New Project File')
        
        # Append scans to an existing integrator window
        self.appendIntegrator = wx.Button(self.rightPanel,
                                          label='Append Scans To...')
        
        # Populate the sizer
        self.nextStepSizer.Add(wx.StaticLine(self.rightPanel, size=(2, 24)),
                               flag=wx.LEFT | wx.RIGHT, border=1)
        self.nextStepSizer.Add(self.newIntegrator, proportion=1,
                               flag=wx.EXPAND | wx.LEFT | wx.RIGHT, border=16)
        self.nextStepSizer.Add(self.appendIntegrator, proportion=1,
                               flag=wx.EXPAND | wx.LEFT | wx.RIGHT, border=16)
        
        # Arrange the right panel
        self.rightSizer.Add(self.projectNameSizer, proportion=0,
                            flag=wx.EXPAND | wx.TOP, border=20)
        self.rightSizer.Add(self.newProjectTree, proportion=1,
                             flag=wx.EXPAND | wx.TOP | wx.BOTTOM,
                             border=4)
        self.rightSizer.Add(self.nextStepSizer, proportion=0,
                            flag=wx.EXPAND | wx.TOP | wx.BOTTOM,
                            border=4)
                            
        self.rightPanel.SetSizerAndFit(self.rightSizer)
        
        ###############################################################
        # End of the right layout
        ###############################################################
        
        self.fullSizer.Add(self.leftPanel, proportion=6,
                           flag=wx.EXPAND | wx.LEFT, border=8)
        self.fullSizer.Add(self.middlePanel, proportion=3, flag=wx.EXPAND)
        self.fullSizer.Add(self.rightPanel, proportion=3,
                           flag=wx.EXPAND | wx.RIGHT, border=8)
        
        self.fullWindow.SetSizer(self.fullSizer)
        
        ###############################################################
        # End of window arrangement
        ###############################################################
        
        # Make the menu bar
        self.menuBar = wx.MenuBar()
        
        # The file menu
        self.fileMenu = wx.Menu()
        self.loadFile = self.fileMenu.Append(-1, 'Load HDF file...')
        self.loadAttr = self.fileMenu.Append(-1, 'Load attribute file...')
        self.exitWindow = self.fileMenu.Append(-1, 'Exit')
        
        self.menuBar.Append(self.fileMenu, 'File')
        
        self.SetMenuBar(self.menuBar)
        
        ###############################################################
        # Start bindings
        ###############################################################
        
        # Menu bindings
        self.Bind(wx.EVT_MENU, self.loadHDF, self.loadFile)
        self.Bind(wx.EVT_MENU, self.loadAttributes, self.loadAttr)
        self.Bind(wx.EVT_MENU, self.onClose, self.exitWindow)
        
        # Button bindings
        # The specfile / number button
        self.specButton.Bind(wx.EVT_BUTTON, self.filterSpec)
        # The HK button
        self.hkButton.Bind(wx.EVT_BUTTON, self.filterHK)
        # The L button
        self.LButton.Bind(wx.EVT_BUTTON, self.filterL)
        # The type button
        self.typeButton.Bind(wx.EVT_BUTTON, self.filterType)
        # The date button
        self.dateButton.Bind(wx.EVT_BUTTON, self.filterDate)
        # The keep button
        self.keepButton.Bind(wx.EVT_BUTTON, self.keepSelected)
        # The reset button
        self.resetButton.Bind(wx.EVT_BUTTON, self.readFile)
        # The load project file button
        self.fileButton.Bind(wx.EVT_BUTTON, self.loadHDF)
        # The single right button
        self.rightOne.Bind(wx.EVT_BUTTON, self.moveOneRight)
        # The single left button
        self.leftOne.Bind(wx.EVT_BUTTON, self.moveOneLeft)
        # The all right button
        self.rightAll.Bind(wx.EVT_BUTTON, self.moveAllRight)
        # The all left button
        self.leftAll.Bind(wx.EVT_BUTTON, self.moveAllLeft)
        # The new attribute button
        self.attrAdd.Bind(wx.EVT_BUTTON, self.addAttribute)
        # The load attribute file button
        self.loadAttrButton.Bind(wx.EVT_BUTTON, self.loadAttributes)
        # The save attribute file button
        self.saveAttrButton.Bind(wx.EVT_BUTTON, self.saveAttributes)
        # The new button
        self.newIntegrator.Bind(wx.EVT_BUTTON, self.newProject)
        # The append button
        self.appendIntegrator.Bind(wx.EVT_BUTTON, self.appendTo)
        
        # Lost focus on the project name box (rename project tree)
        self.projectNameBox.Bind(wx.EVT_KILL_FOCUS, self.newName)
        # dataTable click
        self.dataTable.Bind(wx.EVT_LIST_ITEM_SELECTED, self.tableClick)
        
        ###############################################################
        # End bindings
        ###############################################################
        
        # Because the 'Filter Specfile' button is created first,
        # it automatically gets the focus. As a result, the button
        # glows blue when the window first appears. To change this,
        # we put the focus on a static text element before showing
        # the window.
        self.moreLabel.SetFocus()
        self.Show()
        
    # Choose and load an HDF file so it can be filtered
    def loadHDF(self, event):
        '''Open an HDF file and parse its contents into the filter.'''
        
        loadDialog = wx.FileDialog(self, message='Load file...',
                                   defaultDir=os.getcwd(), defaultFile='',
                                   wildcard='Master files (*.mh5)|*.mh5|'+\
                                             'All files (*.*)|*',
                                   style=wx.OPEN)
        if loadDialog.ShowModal() == wx.ID_OK:
            print 'Loading ' + loadDialog.GetPath()
            try:
                self.filterFile.close()
                del self.filterFile
                self.filterFile = None
            except:
                pass
            try:
                self.filterFile = h5py.File(loadDialog.GetPath(), 'r')
                #self.set_data('filter_file', self.filterFile)
            except:
                print 'Error reading file'
                loadDialog.Destroy()
                raise
            
            self.readFile(None)
            self.fileButton.SetLabel(os.path.split(loadDialog.GetPath())[-1])
            if self.projectNameBox.GetValue() == '':
                projName = os.path.basename(loadDialog.GetPath())
                projName = projName.rsplit('.', 1)[0] + '.ph5'
                self.projectNameBox.SetValue(projName)
            self.newProjectTree.DeleteAllItems()
            self.projectDict = {}
        loadDialog.Destroy()
        self.dataTable.SetFocus()
    
    # Reset all filters and read in an HDF file
    def readFile(self, event):
        '''Read the HDF file self.filterFile, building the 
        scan list self.scanItems to place into the filter list.
        
        '''
        
        if self.filterFile is None:
            print 'Error: no file selected'
            return
            
        # Reset all the hdf-related variables
        self.scanItems = []
        self.allSpecs = {}
        self.allHK = {}
        self.LMin, self.LMax = (float('inf'), float('-inf'))
        self.allTypes = {}
        self.possibleYears = []
        self.dateMin, self.dateMax = (float('inf'), float('-inf'))
        self.activeFilters = []
        self.specResult = None
        self.specCases = None
        self.hkResult = None
        self.hkCases = None
        self.LResult = None
        self.LCases = None
        self.typeResult = None
        self.typeCases = None
        self.dateResult = None
        self.dateCases = None
        
        filterItems = self.filterFile.items()
        filterItems.sort()
        for spec, group in filterItems:
            for number, scan in group.items():
                scanAttrs = scan.attrs
                sAbort = scanAttrs.get('aborted', '?')
                if sAbort == 0:
                    sAbort = ''
                elif sAbort == 1:
                    sAbort = 'True'
                self.scanItems.append([scanAttrs.get('spec_name', 'N/A'),
                                       str(scanAttrs.get('index', '0')),
                                       str(scanAttrs.get('h_val', '--')),
                                       str(scanAttrs.get('k_val', '--')),
                                       str(scanAttrs.get('real_L_start', '--')),
                                       str(scanAttrs.get('real_L_stop', '--')),
                                       scanAttrs.get('s_type', 'N/A'),
                                       sAbort,
                                       scanAttrs.get('date', 'N/A'),
                                       scanAttrs.get('hk_dist', '--'),
                                       scan.name])
        # A bit more list formatting, as well as information gathering
        for scan in self.scanItems:
            if scan[6] not in self.allTypes:
                self.allTypes[scan[6]] = 1
            else:
                self.allTypes[scan[6]] += 1
            # Unless it's a rodscan, don't bother cluttering
            # the display with the L start and stop values
            if scan[6] != 'rodscan':
                scan[4] = '--'
                scan[5] = '--'
            else:
                specVal = scan[0]
                hVal = scan[2]
                kVal = scan[3]
                if (hVal, kVal) not in self.allHK:
                    self.allHK[(hVal, kVal)] = {specVal: [scan[9], 1]}
                elif specVal not in self.allHK[(hVal, kVal)]:
                    self.allHK[(hVal, kVal)][specVal] = [scan[9], 1]
                else:
                    self.allHK[(hVal, kVal)][specVal][1] += 1
            # Make sure the specfile names are saved
            if scan[0] not in self.allSpecs:
                self.allSpecs[scan[0]] = [int(scan[1])]
            else:
                self.allSpecs[scan[0]].append(int(scan[1]))
            # Establish the min and max L values:
            try:
                #if float(scan[4]) < self.LMin: self.LMin = float(scan[4])
                self.LMin = min(float(scan[4]), self.LMin)
            except:
                pass
            try:
                #if float(scan[5]) > self.LMax: self.LMax = float(scan[5])
                self.LMax = max(float(scan[5]), self.LMax)
            except:
                pass
            # Make a list of all the years involved in the scans
            if scan[8].split()[-1] not in self.possibleYears:
                self.possibleYears.append(scan[8].split()[-1])
            # Establish the min and max date epochs:
            try:
                self.dateMin = min(time.mktime(time.strptime(scan[8])),
                                   self.dateMin)
            except:
                pass
            try:
                self.dateMax = max(time.mktime(time.strptime(scan[8])),
                                   self.dateMax)
            except:
                pass
        # Sort the scans by index first
        self.scanItems.sort(key=lambda scan : int(scan[1]))
        # Then sort the scans by specfile name, resulting in 
        # scans ordered by specfile first, then scan number
        self.scanItems.sort(key=lambda scan : scan[0])
        # Update the data table with the new scan information
        self.updateTable()
    
    # Delete everything in the table, then add the appropriate scans
    def updateTable(self):
        '''Delete all the scans from the filter list,
        then repopulate it with the scans that pass
        all of the active filters.
        
        '''
        
        self.dataTable.DeleteAllItems()
        self.activeFilters = []
        for filter in [self.specCases, self.hkCases,
                       self.LCases, self.typeCases, self.dateCases]:
            if filter is not None:
                self.activeFilters.append(filter)
        if self.activeFilters != []:
            self.activeFilters = ft.list_intersect(*self.activeFilters)
            for entry in self.scanItems:
                if entry[10] in self.activeFilters:
                    self.dataTable.Append(entry)
                    try:
                        if self.projectDict[entry[0]][entry[1]] is not None:
                            item = self.dataTable.GetItemCount()-1
                            self.dataTable.SetItemTextColour(item,
                                                             wx.Color(128,
                                                                      128,
                                                                      128))
                    except:
                        pass
        else:
            for entry in self.scanItems:
                self.dataTable.Append(entry)
                try:
                    if self.projectDict[entry[0]][entry[1]] is not None:
                        item = self.dataTable.GetItemCount()-1
                        self.dataTable.SetItemTextColour(item,
                                                         wx.Color(128,
                                                                  128,
                                                                  128))
                except:
                    pass
        return
    
    # Update the 'More Info:' panel when a selection is made
    def tableClick(self, event):
        '''When a selection is made in the scan list,
        update the 'More Info:' panel with relevant information.
        
        '''
        
        tableSelection = event.GetItem()
        specName = self.dataTable.GetItem(tableSelection.m_itemId, 0).Text
        scanNumber = self.dataTable.GetItem(tableSelection.m_itemId, 1).Text
        selectionPath = '/'+specName+'/'+scanNumber
        selectionItem = self.filterFile[selectionPath]
        selectionAttrs = selectionItem.attrs
        toShow = 'Command: ' + selectionAttrs.get('cmd', 'N/A') + \
                  '\n\n' + \
                  'Attenuators: ' + selectionAttrs.get('atten', 'N/A') + \
                  '\n\n' + \
                  'Energy: ' + str(selectionAttrs.get('energy', 'N/A')) + \
                  '\n\n' + \
                  'Data points: ' + str(selectionAttrs.get('nl_dat', 'N/A'))
        self.moreBox.SetValue(toShow)
    
    # Add an attribute both to the attribute dictionary and the on-screen list
    def addAttribute(self, event):
        '''When the '+' button is clicked, add the corresponding
        attribute and value to the list (if the attribute isn't
        already present).
        
        '''
        
        attrText = self.attrSelect.GetStringSelection()
        if attrText == '' or attrText in self.attrDict:
            return
        attrValue = str(self.attrSpecify.GetValue())
        if attrValue == '':
            return
        self.attrDict[attrText] = attrValue
        self.attrList.Append([attrText, attrValue, ''])
        self.attrList.SetItemData(self.attrList.GetItemCount()-1, attrText)
        thisButton = wx.Button(self.attrList, label='X', size=(32, 15),
                               name=attrText)
        thisButton.Bind(wx.EVT_BUTTON, self.deleteMe)
        self.attrList.SetItemWindow(self.attrList.GetItemCount()-1, col=2,
                                    wnd=thisButton)
        self.attrList.SortItems(cmp)
    
    # Delete an attribute both from the dictionary and the on-screen list
    def deleteMe(self, event):
            deleteThis = event.GetEventObject().GetName()
            del self.attrDict[deleteThis]
            deleteThis = self.attrList.FindItemData(-1, deleteThis)
            self.attrList.DeleteItem(deleteThis)
    
    # Load a tab-delimited file of attributes
    def loadAttributes(self, event):
        '''Load in a file with attribute / value pairs, separated by
        a tab.
        
        '''
        loadDialog = wx.FileDialog(self, message='Load file...',
                                   defaultDir=os.getcwd(), defaultFile='',
                                   wildcard='txt files (*.txt)|*.txt|'+\
                                            'All files (*.*)|*',
                                   style=wx.OPEN)
        if loadDialog.ShowModal() == wx.ID_OK:
            print 'Loading attribute file ' + loadDialog.GetPath()
            try:
                attributeFile = open(loadDialog.GetPath())
            except:
                print 'Error opening attribute file'
                loadDialog.Destroy()
                raise
            try:
                for line in attributeFile:
                    line = line.strip().split('\t')
                    if len(line) != 2 or \
                            line[0] not in POSSIBLE_ATTRIBUTES or \
                            line[0] in self.attrDict.keys() or \
                            line[1] is '':
                        continue
                    self.attrDict[line[0]] = line[1]
                    self.attrList.Append([line[0], line[1], ''])
                    self.attrList.SetItemData(self.attrList.GetItemCount()-1,
                                              line[0])
                    thisButton = wx.Button(self.attrList, label='X',
                                           size=(32, 15), name=line[0])
                    thisButton.Bind(wx.EVT_BUTTON, self.deleteMe)
                    self.attrList.SetItemWindow(self.attrList.GetItemCount()-1,
                                                col=2, wnd=thisButton)
                attributeFile.close()
                self.attrList.SortItems(cmp)
            except:
                print 'Error reading attribute file'
                loadDialog.Destroy()
                raise
        loadDialog.Destroy()
                    
    # Save the current attribute table into a tab-delimited file
    def saveAttributes(self, event):
        '''Save a file with attribute / value pairs, separated by
        a tab.
        
        '''
        
        saveDialog = wx.FileDialog(self, message='Save file...',
                                   defaultDir=os.getcwd(), defaultFile='',
                                   wildcard='txt files (*.txt)|*.txt|'+\
                                            'All files (*.*)|*',
                                   style=wx.SAVE)
        if saveDialog.ShowModal() == wx.ID_OK:
            print 'Saving attribute file ' + saveDialog.GetPath()
            try:
                attributeFile = open(saveDialog.GetPath(), 'a')
            except:
                print 'Error opening attribute file'
                saveDialog.Destroy()
                raise
            try:
                for key, value in self.attrDict.iteritems():
                    attributeFile.write(key + '\t' + value + '\n')
            except:
                print 'Error writing to file'
                attributeFile.close()
                saveDialog.Destroy()
                raise
            attributeFile.close()
        saveDialog.Destroy()
    
    # Convert a dictionary to a tree
    def dictToTree(self, thisDict, thisTree, thisRoot):
        for key, value in thisDict.items():
            if type(value) == dict:
                parentItem = thisTree.AppendItem(thisRoot, key)
                self.dictToTree(value, thisTree, parentItem)
            else:
                thisTree.AppendItem(thisRoot, str(key) + ': ' + str(value))
        thisTree.SortChildren(thisRoot)
    
    # Update the file name in the project tree
    def newName(self, event):
        if not self.projectNameBox.GetValue().endswith('.ph5'):
            self.projectNameBox.SetValue(self.projectNameBox.GetValue()+'.ph5')
        try:
            self.newProjectTree.SetItemText(self.newProjectTree.GetRootItem(),
                                            self.projectNameBox.GetValue())
        except:
            pass
    
    # Move the selected scans to the project tree
    def moveOneRight(self, event):
        item = self.dataTable.GetFirstSelected()
        if item == -1:
            return
        while item != -1:
            specName = self.dataTable.GetItem(item, 0).Text
            scanNumber = self.dataTable.GetItem(item, 1).Text
            if specName in self.projectDict:
                if scanNumber not in self.projectDict[specName]:
                    self.projectDict[specName][scanNumber] = \
                                                            self.attrDict.copy()
                else:
                    print 'Specfile ' + specName + ', scan ' + \
                          scanNumber + ' is already in the tree.'
            else:
                self.projectDict[specName] = {}
                self.projectDict[specName][scanNumber] = self.attrDict.copy()
            self.dataTable.SetItemTextColour(item, wx.Color(128, 128, 128))
            item = self.dataTable.GetNextSelected(item)
        self.newProjectTree.DeleteAllItems()
        projectRoot = \
                self.newProjectTree.AddRoot(self.projectNameBox.GetValue())
        self.dictToTree(self.projectDict, self.newProjectTree, projectRoot)
        self.newProjectTree.Expand(projectRoot)

    # Move all the scans to the project tree
    def moveAllRight(self, event):
        item = self.dataTable.GetNextItem(-1)
        if item == -1:
            return
        while item != -1:
            specName = self.dataTable.GetItem(item, 0).Text
            scanNumber = self.dataTable.GetItem(item, 1).Text
            if specName in self.projectDict:
                if scanNumber not in self.projectDict[specName]:
                    self.projectDict[specName][scanNumber] = \
                                                            self.attrDict.copy()
                else:
                    print 'Specfile ' + specName + ', scan ' + \
                          scanNumber + ' is already in the tree.'
            else:
                self.projectDict[specName] = {}
                self.projectDict[specName][scanNumber] = self.attrDict.copy()
            self.dataTable.SetItemTextColour(item, wx.Color(128, 128, 128))
            item = self.dataTable.GetNextItem(item)
        self.newProjectTree.DeleteAllItems()
        projectRoot = \
                self.newProjectTree.AddRoot(self.projectNameBox.GetValue())
        self.dictToTree(self.projectDict, self.newProjectTree, projectRoot)
        self.newProjectTree.Expand(projectRoot)
    
    # Move the selected scans out of the project tree
    def moveOneLeft(self, event):
        allSelected = self.newProjectTree.GetSelections()
        allSelected.sort(key = lambda selection: \
                                        self.newProjectTree.getLevel(selection))
        '''for item in allSelected:
            print self.newProjectTree.GetItemText(item)
        '''
        for selection in allSelected:
            try:
                selectionLevel = self.newProjectTree.getLevel(selection)
            except:
                selectionLevel = -1
            if selectionLevel == 0:
                self.moveAllLeft(event)
                return
            elif selectionLevel == 1:
                specName = self.newProjectTree.GetItemText(selection)
                item = self.dataTable.GetNextItem(-1)
                while item != -1:
                    if self.dataTable.GetItem(item, 0).Text == specName:
                        self.dataTable.SetItemTextColour(item, wx.BLACK)
                    item = self.dataTable.GetNextItem(item)
                self.projectDict.pop(specName)
                self.newProjectTree.Delete(selection)
            elif selectionLevel == 2:
                scanParent = self.newProjectTree.GetItemParent(selection)
                scanNumber = self.newProjectTree.GetItemText(selection)
                specName = self.newProjectTree.GetItemText(scanParent)
                item = self.dataTable.GetNextItem(-1)
                while item != -1:
                    if self.dataTable.GetItem(item, 0).Text == specName and \
                            self.dataTable.GetItem(item, 1).Text == scanNumber:
                        self.dataTable.SetItemTextColour(item, wx.BLACK)
                        break
                    item = self.dataTable.GetNextItem(item)
                self.projectDict[specName].pop(scanNumber)
                self.newProjectTree.Delete(selection)
            else:
                pass
        popUs = []
        for key in self.projectDict:
            if self.projectDict[key] == {}:
                popUs.append(key)
        for key in popUs:
            self.projectDict.pop(key)
        projectRoot = self.newProjectTree.GetRootItem()
        if not projectRoot.IsOk():
            return
        if self.newProjectTree.GetChildrenCount(projectRoot) == 0:
            self.newProjectTree.Delete(projectRoot)
        else:
            item, cookie = self.newProjectTree.GetFirstChild(projectRoot)
            while item.IsOk():
                if self.newProjectTree.GetChildrenCount(item) == 0:
                    self.newProjectTree.Delete(item)
                item, cookie = self.newProjectTree.GetNextChild(item, cookie)
    
    # Move all the scans out of the project tree
    def moveAllLeft(self, event):
        self.newProjectTree.DeleteAllItems()
        self.projectDict = {}
        item = self.dataTable.GetNextItem(-1)
        while item != -1:
            self.dataTable.SetItemTextColour(item, wx.BLACK)
            item = self.dataTable.GetNextItem(item)
    
    # Create a new HDF project file, then open it in an integrator window
    def newProject(self, event):
        #print self.projectDict
        if self.projectDict == {}:
            print 'No scans selected'
            return
        self.fileDirectory, holding = os.path.split(self.filterFile.filename)
        file_types = 'Project files (*.ph5)|*.ph5|All files (*.*)|*'
        save_dialog = wx.FileDialog(self, message='Create file...',
                                    defaultDir=self.fileDirectory,
                                    defaultFile=self.projectNameBox.GetValue(),
                                    wildcard=file_types,
                                    style=wx.SAVE | wx.OVERWRITE_PROMPT)
        if save_dialog.ShowModal() == wx.ID_OK:
            print 'Start: ', time.ctime(time.time())
            out_file = save_dialog.GetPath()
            try:
                mtp.master_to_project(self.filterFile.filename,
                                      self.projectDict,
                                      out_file, append=False, gui=True)
            except:
                print 'Error generating project file'
            print 'Finish: ', time.ctime(time.time())
        save_dialog.Destroy()
    
    # Append scans to an existing HDF project file, then open it
    def appendTo(self, event):
        if self.projectDict == {}:
            print 'No scans selected'
            return
        self.fileDirectory, holding = os.path.split(self.filterFile.filename)
        file_types = 'Project files (*.ph5)|*.ph5|All files (*.*)|*'
        save_dialog = wx.FileDialog(self, message='Append to...',
                                    defaultDir=self.fileDirectory,
                                    defaultFile=self.projectNameBox.GetValue(),
                                    wildcard=file_types,
                                    style=wx.SAVE)
        if save_dialog.ShowModal() == wx.ID_OK:
            print 'Start: ', time.ctime(time.time())
            out_file = save_dialog.GetPath()
            try:
                mtp.master_to_project(self.filterFile.filename,
                                      self.projectDict,
                                      out_file, append=True, gui=True)
            except:
                print 'Error saving to project file'
            print 'Finish: ', time.ctime(time.time())
        save_dialog.Destroy()
    
    # Keep only the selected scans by faking a filtered result
    def keepSelected(self, event):
        item = self.dataTable.GetFirstSelected()
        if item == -1:
            return
        self.specResult = {}
        self.specCases = []
        while item != -1:
            specName = self.dataTable.GetItem(item, 0).Text
            scanNumber = self.dataTable.GetItem(item, 1).Text
            if specName in self.specResult:
                self.specResult[specName].append(int(scanNumber))
            else:
                self.specResult[specName] = [int(scanNumber)]
            self.specCases.append('/'+specName+'/'+scanNumber)
            item = self.dataTable.GetNextSelected(item)
        for key in self.allSpecs:
            if key not in self.specResult:
                self.specResult[key] = []
        self.updateTable()
    
    # Filter by specfile and scan number
    def filterSpec(self, event):
        filterWindow = SpecWindow(self,
                                  self.allSpecs,
                                  self.specResult).ShowModal()
        return
    
    # Filter by HK pair
    def filterHK(self, event):
        filterWindow = HKWindow(self, self.allHK, self.hkResult).ShowModal()
        return
    
    # Filter by L range
    def filterL(self, event):
        filterWindow = LWindow(self,
                               (self.LMin, self.LMax),
                               self.LResult).ShowModal()
        return
        
    # Filter by scan type
    def filterType(self, event):
        filterWindow = TypeWindow(self,
                                  self.allTypes,
                                  self.typeResult).ShowModal()
        return
    
    # Filter by date range
    def filterDate(self, event):
        filterWindow = DateWindow(self,
                                  self.possibleYears,
                                  (self.dateMin, self.dateMax),
                                  self.dateResult).ShowModal()
        return
        
    # Close the window
    def onClose(self, event):
        try:
            self.filterFile.close()
        except:
            pass
        self.Destroy()
        del self


class SpecWindow(wx.Dialog):
    '''The GUI for filtering by specfile.'''
    def __init__(self, parent=None, allSpec={}, currentSpec=None):
        
        # Make the window
        wx.Dialog.__init__(self, parent, -1,
                           title="Specfiles", size=(240, 400))
        
        # Make the list
        self.list = treemix.CustomTreeCtrl(self, agwStyle=
                                                  treemix.TR_HAS_BUTTONS |
                                                  treemix.TR_MULTIPLE |
                                                  treemix.TR_EXTENDED | 
                                                  treemix.TR_AUTO_CHECK_CHILD |
                                                  treemix.TR_AUTO_CHECK_PARENT)
        
        # Populate the list
        self.allSpec = allSpec
        if currentSpec is not None:
            self.currentSpec = currentSpec
        else:
            self.currentSpec = allSpec
        self.allRoot = self.list.AddRoot(self.GetParent().fileButton.GetLabel(),
                                         ct_type=1)
        specKeys = self.allSpec.keys()
        specKeys.sort()
        for spec in specKeys:
            specRoot = self.list.AppendItem(self.allRoot, spec, ct_type=1)
            allScans = self.allSpec[spec]
            allScans.sort()
            for scan in allScans:
                scanRoot = self.list.AppendItem(specRoot, str(scan), ct_type=1)
                if scan in self.currentSpec[spec]:
                    self.list.CheckItem(scanRoot)
        self.list.Expand(self.allRoot)
        
        # Create the apply and cancel buttons
        self.apply = wx.Button(self, label="Apply")
        self.cancel = wx.Button(self, label="Cancel")
        
        # Bindings
        self.apply.Bind(wx.EVT_BUTTON, self.onApply)
        self.cancel.Bind(wx.EVT_BUTTON, self.onCancel)
        
        # General layout
        self.windowSizer = wx.BoxSizer(wx.VERTICAL)
        # The apply and cancel buttons
        self.buttonSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # The list
        self.windowSizer.Add(self.list, proportion=1,
                             flag=wx.EXPAND | wx.ALL, border=4)
        
        # The apply and cancel buttons
        self.buttonSizer.Add(self.apply, proportion=1,
                             flag=wx.EXPAND | wx.ALL, border=4)
        self.buttonSizer.Add(self.cancel, proportion=1,
                             flag=wx.EXPAND | wx.ALL, border=4)
        self.windowSizer.Add(self.buttonSizer,
                             flag=wx.EXPAND | wx.ALL, border=4)
        
        # Create the window, center it, and bind the close button
        self.SetSizer(self.windowSizer)
        self.CenterOnScreen()
        self.Bind(wx.EVT_CLOSE, self.onCancel)
    
    # First, release the focus so if something breaks you're not stuck.
    # Then, iterate through the specfiles, filling out the selectedSpec
    # dictionary with selected scan numbers. With this done, filter the
    # hdf file for the matching scans and update the table in the main window.
    def onApply(self, event):
        wx.Window.MakeModal(self, False)
        selectedSpec = {}
        specKids = self.allRoot.GetChildren()
        for spec in specKids:
            scanKids = spec.GetChildren()
            specText = spec.GetText()
            selectedSpec[specText] = []
            for scan in scanKids:
                if scan.IsChecked():
                    selectedSpec[specText].append(int(scan.GetText()))
        self.GetParent().specCases = None
        if self.allSpec == selectedSpec:
            self.GetParent().specResult = None
        else:
            self.GetParent().specResult = selectedSpec
            for spec in selectedSpec.keys():
                bothCases = []
                if selectedSpec[spec] != []:
                    thisCase = ft.cases(self.GetParent().filterFile,
                                        'spec_name',
                                        '.startswith("' + spec + '")')
                    thatCase = ft.cases(self.GetParent().filterFile, 'index',
                                        'in ' + str(selectedSpec[spec]))
                    bothCases = ft.list_intersect(thisCase, thatCase)
                    self.GetParent().specCases = \
                            ft.list_union(bothCases, self.GetParent().specCases)
        self.GetParent().updateTable()
        self.Destroy()
        del self
    
    # Release the focus and close the window, making no changes
    def onCancel(self, event):
        wx.Window.MakeModal(self, False)
        self.Destroy()
        del self


class HKWindow(wx.Dialog):
    '''The GUI for filtering by HK pair.'''
    def __init__(self, parent=None, allHK={}, currentHK=None):
        
        # Make the window
        wx.Dialog.__init__(self, parent, -1,
                           title="HK Pairs", size=(400, 400))
                           
        # Make the list
        self.list = treemix.CustomTreeCtrl(self, agwStyle=
                                                  treemix.TR_HAS_BUTTONS |
                                                  treemix.TR_MULTIPLE |
                                                  treemix.TR_EXTENDED |
                                                  treemix.TR_AUTO_CHECK_CHILD |
                                                  treemix.TR_AUTO_CHECK_PARENT)
        
        # Populate the list
        self.allHK = allHK
        if currentHK is not None:
            self.currentHK = currentHK
        else:
            self.currentHK = allHK
        self.allRoot = self.list.AddRoot('HK Pairs', ct_type=1)
        self.allHKKeys = self.allHK.keys()
        self.allHKKeys = sorted(self.allHKKeys, key=lambda key: eval(key[1]))
        self.allHKKeys = sorted(self.allHKKeys, key=lambda key: eval(key[0]))
        for hk in self.allHKKeys:
            hkRoot = self.list.AppendItem(self.allRoot,
                                          str(tuple(map(float, hk))), ct_type=1)
            allSpecs = self.allHK[hk].keys()
            allSpecs.sort()
            for spec in allSpecs:
                specRoot = self.list.AppendItem(hkRoot,
                                                spec + ': ' + 
                                                str(self.allHK[hk][spec][1]) + \
                                                ' at a distance of ' + \
                                                str(self.allHK[hk][spec][0]),
                                                ct_type=1)
                if spec in self.currentHK[hk]:
                    self.list.CheckItem(specRoot)
        self.list.Expand(self.allRoot)
        
        # Create the apply and cancel buttons
        self.apply = wx.Button(self, label="Apply")
        self.cancel = wx.Button(self, label="Cancel")
        
        # Bindings
        self.apply.Bind(wx.EVT_BUTTON, self.onApply)
        self.cancel.Bind(wx.EVT_BUTTON, self.onCancel)
    
        # General layout
        self.windowSizer = wx.BoxSizer(wx.VERTICAL)
        # The apply and cancel buttons
        self.buttonSizer = wx.BoxSizer(wx.HORIZONTAL)

        # The list
        self.windowSizer.Add(self.list, proportion=1,
                             flag=wx.EXPAND | wx.ALL, border=4)
    
        # The apply and cancel buttons
        self.buttonSizer.Add(self.apply, proportion=1,
                             flag=wx.EXPAND | wx.ALL, border=4)
        self.buttonSizer.Add(self.cancel, proportion=1,
                             flag=wx.EXPAND | wx.ALL, border=4)
        self.windowSizer.Add(self.buttonSizer,
                             flag=wx.EXPAND | wx.ALL, border=4)
        
        # Create the window, center it, and bind the close button
        self.SetSizer(self.windowSizer)
        self.CenterOnScreen()
        self.Bind(wx.EVT_CLOSE, self.onCancel)

    # First, release the focus so if something breaks you're not stuck.
    # Then, fill out selectedHK to match the format of allHKs from the 
    # parent class. Once this is done, filter the hdf file to return the
    # desired scans and update the table in the main window.
    def onApply(self, event):
        wx.Window.MakeModal(self, False)
        selectedHK = {}
        hkKids = self.allRoot.GetChildren()
        for hk in hkKids:
            specKids = hk.GetChildren()
            hkText = eval(hk.GetText())
            selectedHK[hkText] = {}
            for spec in specKids:
                if spec.IsChecked():
                    specName = spec.GetText().split(':')[0]
                    allHKText = tuple(map(str, hkText))
                    selectedHK[hkText][specName] = \
                                    self.allHK[allHKText][specName]
        self.GetParent().hkCases = None
        if self.allHK == selectedHK:
            self.GetParent().hkResult = None
        else:
            self.GetParent().hkResult = selectedHK
            for hk in selectedHK.keys():
                bothCases = []
                if selectedHK[hk] != {}:
                    thisCase = ft.cases(self.GetParent().filterFile,
                                        'h_val',
                                        '== ' + str(hk[0]))
                    thatCase = ft.cases(self.GetParent().filterFile,
                                        'k_val',
                                        '== ' + str(hk[1]))
                    otherCase = ft.cases(self.GetParent().filterFile,
                                         'spec_name',
                                         'in ' + str(selectedHK[hk].keys()))
                    bothCases = ft.list_intersect(thisCase, thatCase, otherCase)
                    self.GetParent().hkCases = \
                            ft.list_union(bothCases, self.GetParent().hkCases)
        self.GetParent().updateTable()
        self.Destroy()
        del self
    
    # Release the focus and close the window, making no changes
    def onCancel(self, event):
        wx.Window.MakeModal(self, False)
        self.Destroy()
        del self


class LWindow(wx.Dialog):
    '''The GUI for filtering by L value.'''
    def __init__(self, parent=None, allL=(-100, 100), currentL=None):
    
        # Make the window
        wx.Dialog.__init__(self, parent, -1,
                           title="L Range", size=(232, 100))
        
        self.LMin, self.LMax = allL
        if currentL is not None:
            self.currentLMin, self.currentLMax = currentL
        else:
            self.currentLMin, self.currentLMax = allL
        
        # 'From' value components
        self.fromText = wx.StaticText(self, label="From:")
        self.fromValue = wx.TextCtrl(self, size=(60, -1))
        self.fromValue.SetValue(str(self.currentLMin))
        
        # 'To' value components
        self.toText = wx.StaticText(self, label="To:")
        self.toValue = wx.TextCtrl(self, size=(60, -1))
        self.toValue.SetValue(str(self.currentLMax))
        
        # Buttons (apply and cancel)
        self.apply = wx.Button(self, label="Apply")
        self.cancel = wx.Button(self, label="Cancel")
        
        # Bindings
        self.apply.Bind(wx.EVT_BUTTON, self.onApply)
        self.cancel.Bind(wx.EVT_BUTTON, self.onCancel)
        
        # General layout
        self.windowSizer = wx.BoxSizer(wx.VERTICAL)
        # 'From' and 'To' layout
        self.fromToSizer = wx.BoxSizer(wx.HORIZONTAL)
        # Button layout
        self.buttonSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # 'From' input followed by 'To' input
        self.fromToSizer.Add(self.fromText, flag=wx.EXPAND | wx.TOP |
                                                 wx.BOTTOM | wx.LEFT, border=8)
        self.fromToSizer.Add(self.fromValue, flag=wx.EXPAND | wx.ALL, border=4)
        self.fromToSizer.Add(self.toText, flag=wx.EXPAND | wx.TOP |
                                               wx.BOTTOM | wx.LEFT, border=8)
        self.fromToSizer.Add(self.toValue, flag=wx.EXPAND | wx.ALL, border=4)
        self.windowSizer.Add(self.fromToSizer,
                             flag=wx.EXPAND | wx.ALL, border=4)
        
        # Button positions
        self.buttonSizer.Add(self.apply, proportion=1,
                             flag=wx.EXPAND | wx.RIGHT | wx.BOTTOM | wx.LEFT,
                             border=4)
        self.buttonSizer.Add(self.cancel, proportion=1,
                             flag=wx.EXPAND | wx.RIGHT | wx.BOTTOM | wx.LEFT,
                             border=4)
        self.windowSizer.Add(self.buttonSizer,
                             flag=wx.EXPAND | wx.ALL, border=4)
        
        # Create the window, center it, and bind the close button
        self.SetSizer(self.windowSizer)
        self.CenterOnScreen()
        self.Bind(wx.EVT_CLOSE, self.onCancel)

    # First, release the focus so if something breaks you're not stuck.
    # Then, capture the 'from' and 'to' values and filter the hdf file
    # to pick out scans within the specified range. Also, filter out 
    # anything that's not a rodscan to reduce clutter.
    def onApply(self, event):
        wx.Window.MakeModal(self, False)
        try:
            fromL = float(self.fromValue.GetValue())
        except:
            print 'Error: Invalid minimum L; setting to -100'
            fromL = -100.0
        try:
            toL = float(self.toValue.GetValue())
        except:
            print 'Error: Invalid maximum L; setting to 100'
            toL = 100.0
        self.GetParent().LCases = None
        if fromL <= self.LMin and toL >= self.LMax:
            self.GetParent().LResult = None
        else:
            self.GetParent().LResult = (fromL, toL)
            thisCase = ft.cases(self.GetParent().filterFile,
                                'real_L_start',
                                '>= ' + str(fromL))
            thatCase = ft.cases(self.GetParent().filterFile,
                                'real_L_stop',
                                '<= ' + str(toL))
            otherCase = ft.cases(self.GetParent().filterFile,
                                 's_type',
                                 '.startswith("rodscan")')
            bothCases = ft.list_intersect(thisCase, thatCase, otherCase)
            self.GetParent().LCases = bothCases
        self.GetParent().updateTable()
        self.Destroy()
        del self
    
    #Release the focus and close the window, making no changes
    def onCancel(self, event):
        wx.Window.MakeModal(self, False)
        self.Destroy()
        del self


class TypeWindow(wx.Dialog):
    '''The GUI for filtering by scan type.'''
    def __init__(self, parent=None, allTypes={}, currentTypes=None):
    
        # Make the window
        wx.Dialog.__init__(self, parent, -1,
                           title="Scan Types", size=(200, 400))

        # Make the list
        self.list = NumberListCtrl(self, style=wx.LC_REPORT | wx.LC_NO_HEADER |
                                               wx.LC_HRULES | wx.LC_VRULES)
        self.list.InsertColumn(0, "", width=24)
        self.list.InsertColumn(1, "Scan Type", width=66)
        self.list.InsertColumn(2, "Num. Scans")
        
        # Populate the list with the different scan types; check
        # the appropriate boxes to show currently selected types
        self.allTypes = allTypes
        if currentTypes is not None:
            self.currentTypes = currentTypes
        else:
            self.currentTypes = allTypes
        self.allTypeKeys = self.allTypes.keys()
        self.allTypeKeys.sort()
        for key in self.allTypeKeys:
            newType = self.list.Append(["", key, self.allTypes[key]])
            if key in self.currentTypes:
                self.list.CheckItem(newType)
        
        # Create the 'check all' check box (selector), the apply and
        # cancel buttons, and the input field for scan ranges
        self.selector = wx.CheckBox(self, style=wx.CHK_3STATE)
        self.apply = wx.Button(self, label="Apply")
        self.cancel = wx.Button(self, label="Cancel")
        
        # Set the 'check all' check box to the
        # appropriate state (all, some, or none)
        self.setCheck(None)
        
        # Bindings
        self.apply.Bind(wx.EVT_BUTTON, self.onApply)
        self.cancel.Bind(wx.EVT_BUTTON, self.onCancel)
        self.selector.Bind(wx.EVT_CHECKBOX, self.onCheckAll)
        
        # If checking an item takes too long, try commenting this out
        self.Bind(wx.EVT_CHECKBOX, self.setCheck)
        
        # General layout
        self.windowSizer = wx.BoxSizer(wx.VERTICAL)
        # The column headers
        self.columnSizer = wx.BoxSizer(wx.HORIZONTAL)
        # The apply and cancel buttons
        self.buttonSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # The column headers
        self.columnSizer.Add(self.selector,
                             flag=wx.EXPAND | wx.LEFT | wx.RIGHT, border=4)
        self.columnSizer.Add(wx.StaticText(self, label="Scan Type"),
                             flag=wx.EXPAND | wx.LEFT, border=12)
        self.columnSizer.Add(wx.StaticText(self, label="Num. Scans"),
                             flag=wx.EXPAND | wx.LEFT, border=12)
        self.windowSizer.Add(self.columnSizer,
                             flag=wx.EXPAND | wx.ALL, border=4)
        
        # The list
        self.windowSizer.Add(self.list, proportion=1,
                             flag=wx.EXPAND | wx.ALL, border=4)
        
        # The apply and cancel buttons
        self.buttonSizer.Add(self.apply, proportion=1,
                             flag=wx.EXPAND | wx.ALL, border=4)
        self.buttonSizer.Add(self.cancel, proportion=1,
                             flag=wx.EXPAND | wx.ALL, border=4)
        self.windowSizer.Add(self.buttonSizer,
                             flag=wx.EXPAND | wx.ALL, border=4)
        
        # Create the window, center it, and bind the close button
        self.SetSizer(self.windowSizer)
        self.CenterOnScreen()
        self.Bind(wx.EVT_CLOSE, self.onCancel)
    
    # Examine the state of all the check boxes until enough
    # is known to set the state of the 'check all' box
    def setCheck(self, event):
        num = self.list.GetItemCount()
        if num == 0: return
        thirdState = self.list.IsChecked(0)
        for i in range(1, num):
            if self.list.IsChecked(i) != thirdState:
                self.selector.Set3StateValue(wx.CHK_UNDETERMINED)
                return
        self.selector.SetValue(thirdState)
    
    # Select all or deselect all; if the 'check all' check box
    # is in the mixed state, deselects all
    def onCheckAll(self, event):
        if self.selector.Get3StateValue():
            self.onSelectAll(None)
        else:
            self.onDeselectAll(None)
    
    # Check all check boxes
    def onSelectAll(self, event):
        num = self.list.GetItemCount()
        for i in range(num):
            self.list.CheckItem(i)

    # Uncheck all check boxes
    def onDeselectAll(self, event):
        num = self.list.GetItemCount()
        for i in range(num):
            self.list.CheckItem(i, False)

    # First, release the focus so if something breaks you're not stuck.
    # Then, filter the hdf file to select the checked scans and update
    # the table in the main window.
    def onApply(self, event):
        wx.Window.MakeModal(self, False)
        selectedTypes = []
        num = self.list.GetItemCount()
        for i in range(num):
            if self.list.IsChecked(i):
                selectedTypes.append(self.list.GetItem(i, 1).GetText())
        self.GetParent().typeCases = None
        if set(self.allTypes.keys()) == set(selectedTypes):
            self.GetParent().typeResult = None
        else:
            self.GetParent().typeResult = selectedTypes
            thisCase = ft.cases(self.GetParent().filterFile,
                                's_type',
                                'in ' + str(selectedTypes))
            self.GetParent().typeCases = thisCase
        self.GetParent().updateTable()
        self.Destroy()
        del self
    
    # Release the focus and close the window, making no changes
    def onCancel(self, event):
        wx.Window.MakeModal(self, False)
        self.Destroy()
        del self


class DateWindow(wx.Dialog):
    '''The GUI for filtering by date.'''
    def __init__(self, parent=None, possibleYears=[],
                                    allDates=(float('-inf'), float('inf')),
                                    currentDates=None):
        
        # Make the window
        wx.Dialog.__init__(self, parent, -1,
                           title="Date Range", size=(288, 117))
        
        self.dateMin, self.dateMax = allDates
        if currentDates is not None:
            self.currentDateMin, self.currentDateMax = currentDates
        else:
            self.currentDateMin, self.currentDateMax = allDates
        # Set the default choices for year, month, day, hour, and minute.
        # The year choices are passed as an argument, sorted early to late.
        # If the user selects February 31, it automaticalls sets it to 
        # February 28/29 at 23:59:59 to prevent an error when converting
        # to epoch.
        self.yearChoices = possibleYears
        if self.yearChoices == []:
            self.yearChoices = ['N/A']
        self.monthChoices = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                             'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        self.dayChoices = ['01', '02', '03', '04', '05', '06', '07', '08',
                           '09', '10', '11', '12', '13', '14', '15', '16',
                           '17', '18', '19', '20', '21', '22', '23', '24',
                           '25', '26', '27', '28', '29', '30', '31']
        self.hourChoices = ['00', '01', '02', '03', '04', '05', '06', '07',
                            '08', '09', '10', '11', '12', '13', '14', '15',
                            '16', '17', '18', '19', '20', '21', '22', '23']
        self.minuteChoices = ['00', '01', '02', '03', '04', '05', '06', '07',
                              '08', '09', '10', '11', '12', '13', '14', '15',
                              '16', '17', '18', '19', '20', '21', '22', '23',
                              '24', '25', '26', '27', '28', '29', '30', '31',
                              '32', '33', '34', '35', '36', '37', '38', '39',
                              '40', '41', '42', '43', '44', '45', '46', '47',
                              '48', '49', '50', '51', '52', '53', '54', '55',
                              '56', '57', '58', '59']
        
        # Create the colons for between the hour and minute choices,
        # set to a larger font so they're easier to see
        self.colonText = wx.StaticText(self, label=":")
        self.colonText.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        self.colonText2 = wx.StaticText(self, label=":")
        self.colonText2.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        
        # Buttons
        self.apply = wx.Button(self, label="Apply")
        self.cancel = wx.Button(self, label="Cancel")
        
        # Bindings
        self.apply.Bind(wx.EVT_BUTTON, self.onApply)
        self.cancel.Bind(wx.EVT_BUTTON, self.onCancel)

        # Create the 'From' fields
        self.fromText = wx.StaticText(self, label="From:")
        self.fromMonth = wx.Choice(self, choices=self.monthChoices,
                                   size=(46, 21))
        self.fromDay = wx.Choice(self, choices=self.dayChoices,
                                 size=(36, 21))
        self.fromYear = wx.Choice(self, choices=self.yearChoices,
                                  size=(54, 21))
        self.fromHour = wx.Choice(self, choices=self.hourChoices,
                                  size=(36, 21))
        self.fromMinute = wx.Choice(self, choices=self.minuteChoices,
                                    size=(36, 21))
        
        # Set the initial selections for the 'From' fields
        fromDate = time.ctime(self.currentDateMin).split()
        self.fromMonth.SetStringSelection(fromDate[1])
        self.fromDay.SetStringSelection(fromDate[2])
        self.fromYear.SetStringSelection(fromDate[4])
        self.fromHour.SetStringSelection(fromDate[3].split(':')[0])
        self.fromMinute.SetStringSelection(fromDate[3].split(':')[1])
        
        #Create the 'To' fields
        self.toText = wx.StaticText(self, label="To:")
        self.toMonth = wx.Choice(self, choices=self.monthChoices,
                                 size=(46, 21))
        self.toDay = wx.Choice(self, choices=self.dayChoices,
                               size=(36, 21))
        self.toYear = wx.Choice(self, choices=self.yearChoices,
                                size=(54, 21))
        self.toHour = wx.Choice(self, choices=self.hourChoices,
                                size=(36, 21))
        self.toMinute = wx.Choice(self, choices=self.minuteChoices,
                                  size=(36, 21))
        
        # Set the initial selections for the 'To' fields
        toDate = time.ctime(self.currentDateMax).split()
        self.toMonth.SetStringSelection(toDate[1])
        self.toDay.SetStringSelection(toDate[2])
        self.toYear.SetStringSelection(toDate[4])
        self.toHour.SetStringSelection(toDate[3].split(':')[0])
        minuteBump = str(int(toDate[3].split(':')[1]) + 1)
        self.toMinute.SetStringSelection(minuteBump)
        
        # General layout
        self.windowSizer = wx.BoxSizer(wx.VERTICAL)
        # 'From' values
        self.fromSizer = wx.BoxSizer(wx.HORIZONTAL)
        # 'To' values
        self.toSizer = wx.BoxSizer(wx.HORIZONTAL)
        # Apply and cancel buttons
        self.buttonSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # 'From' sizer
        self.fromSizer.Add(self.fromText, flag=wx.EXPAND | wx.TOP |
                                               wx.BOTTOM | wx.LEFT, border=4)
        self.fromSizer.Add(self.fromMonth, flag=wx.EXPAND | wx.LEFT, border=2)
        self.fromSizer.Add(self.fromDay, flag=wx.EXPAND | wx.LEFT, border=2)
        self.fromSizer.Add(self.fromYear, flag=wx.EXPAND | wx.LEFT, border=2)
        self.fromSizer.Add(self.fromHour, flag=wx.EXPAND | wx.LEFT, border=8)
        self.fromSizer.Add(self.colonText, flag=wx.EXPAND | wx.LEFT, border=0)
        self.fromSizer.Add(self.fromMinute, flag=wx.EXPAND | wx.LEFT, border=0)
        self.windowSizer.Add(self.fromSizer, flag=wx.EXPAND | wx.ALL, border=4)
        
        # 'To' sizer
        self.toSizer.Add(self.toText, flag=wx.EXPAND | wx.TOP |
                                           wx.BOTTOM | wx.LEFT, border=4)
        self.toSizer.Add(self.toMonth, flag=wx.EXPAND | wx.LEFT, border=14)
        self.toSizer.Add(self.toDay, flag=wx.EXPAND | wx.LEFT, border=2)
        self.toSizer.Add(self.toYear, flag=wx.EXPAND | wx.LEFT, border=2)
        self.toSizer.Add(self.toHour, flag=wx.EXPAND | wx.LEFT, border=8)
        self.toSizer.Add(self.colonText2, flag=wx.EXPAND | wx.LEFT, border=0)
        self.toSizer.Add(self.toMinute, flag=wx.EXPAND | wx.LEFT, border=0)
        self.windowSizer.Add(self.toSizer, flag=wx.EXPAND | wx.ALL, border=4)
        
        # Button sizer
        self.buttonSizer.Add(self.apply, proportion=1,
                             flag=wx.EXPAND | wx.ALL, border=4)
        self.buttonSizer.Add(self.cancel, proportion=1,
                             flag=wx.EXPAND | wx.ALL, border=4)
        self.windowSizer.Add(self.buttonSizer,
                             flag=wx.EXPAND | wx.LEFT |
                                  wx.RIGHT | wx.BOTTOM, border=4)
        
        # Create the window, center it, and bind the close button
        self.SetSizer(self.windowSizer)
        self.CenterOnScreen()
        self.Bind(wx.EVT_CLOSE, self.onCancel)

    # First, release the focus so if something breaks you're not stuck.
    # Then, capture the selected 'from' and 'to' dates and filter the 
    # hdf file to select the scans that were initiated in the given range.
    def onApply(self, event):
        wx.Window.MakeModal(self, False)
        # Check the 'From' date to make sure it exists
        if self.fromMonth.GetStringSelection() in ['Sep', 'Apr', 'Jun', 'Nov'] \
                        and self.fromDay.GetStringSelection() == '31':
            self.fromDay.SetStringSelection('30')
            self.fromHour.SetStringSelection('23')
            self.fromMinute.SetStringSelection('59')
        elif self.fromMonth.GetStringSelection() == 'Feb':
            if int(self.fromYear.GetStringSelection()) % 4 == 0:
                if self.fromDay.GetStringSelection() in ['30', '31']:
                    self.fromDay.SetStringSelection('29')
                    self.fromHour.SetStringSelection('23')
                    self.fromMinute.SetStringSelection('59')
            else:
                if self.fromDay.GetStringSelection() in ['29', '30', '31']:
                    self.fromDay.SetStringSelection('28')
                    self.fromHour.SetStringSelection('23')
                    self.fromMinute.SetStringSelection('59')
        # Concatenate the date string and convert it to epoch
        fromDate = 'Mon ' + self.fromMonth.GetStringSelection() + ' ' + \
                    self.fromDay.GetStringSelection() + ' ' + \
                    self.fromHour.GetStringSelection() + ':' + \
                    self.fromMinute.GetStringSelection() + ':00 ' + \
                    self.fromYear.GetStringSelection()
        fromDate = time.mktime(time.strptime(fromDate))
        # Check the 'To' date to make sure it exists
        if self.toMonth.GetStringSelection() in ['Sep', 'Apr', 'Jun', 'Nov'] \
                        and self.toDay.GetStringSelection() == '31':
            self.toDay.SetStringSelection('30')
            self.toHour.SetStringSelection('23')
            self.toMinute.SetStringSelection('59')
        elif self.toMonth.GetStringSelection() == 'Feb':
            if int(self.toYear.GetStringSelection()) % 4 == 0:
                if self.toDay.GetStringSelection() in ['30', '31']:
                    self.toDay.SetStringSelection('29')
                    self.toHour.SetStringSelection('23')
                    self.toMinute.SetStringSelection('59')
            else:
                if self.toDay.GetStringSelection() in ['29', '30', '31']:
                    self.toDay.SetStringSelection('28')
                    self.toHour.SetStringSelection('23')
                    self.toMinute.SetStringSelection('59')
        # Concatenate the date string and convert it to epoch
        toDate = 'Mon ' + self.toMonth.GetStringSelection() + ' ' + \
                    self.toDay.GetStringSelection() + ' ' + \
                    self.toHour.GetStringSelection() + ':' + \
                    self.toMinute.GetStringSelection() + ':00 ' + \
                    self.toYear.GetStringSelection()
        toDate = time.mktime(time.strptime(toDate))
        if fromDate > toDate:
            print "Error: 'To' date precedes 'From' date"
            self.Destroy()
            del self
            return
        self.GetParent().dateCases = None
        if fromDate <= self.dateMin and toDate >= self.dateMax:
            self.GetParent().dateResult = None
        else:
            self.GetParent().dateResult = (fromDate, toDate)
            thisCase = ft.cases(self.GetParent().filterFile,
                                'epoch',
                                '>= ' + str(fromDate))
            thatCase = ft.cases(self.GetParent().filterFile,
                                'epoch',
                                '<= ' + str(toDate))
            bothCases = ft.list_intersect(thisCase, thatCase)
            self.GetParent().dateCases = bothCases
        self.GetParent().updateTable()
        self.Destroy()
        del self
    
    # Release the focus and close the window, making no changes
    def onCancel(self, event):
        wx.Window.MakeModal(self, False)
        self.Destroy()
        del self


class TableDataCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin):
    '''This is a simple wrapper class to combine a standard
        list control with one that autosizes the final column.
    
    '''
    def __init__(self, parent, ID=-1, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)


class NumberListCtrl(wx.ListCtrl, listmix.CheckListCtrlMixin,
                     listmix.ListCtrlAutoWidthMixin):
    '''This class is used to make a list with check boxes'''
    def __init__(self, *args, **kwargs):
        wx.ListCtrl.__init__(self, *args, **kwargs)
        listmix.CheckListCtrlMixin.__init__(self)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
    
    # This function can be used to update the 'check all' check box,
    # but for large numbers ofitems it can be noticeably slow
    def OnCheckItem(self, index, flag):
        wx.PostEvent(self.GetEventHandler(),
                     wx.PyCommandEvent(wx.EVT_CHECKBOX.typeId,
                                       self.GetId()))
        

class myTreeCtrl(wx.TreeCtrl):
    '''Tree class identical to a regular wx.TreeCtrl except for
        an overridden sort function and a getLevel function that
        returns the number of layers down from the root the item
        is (root is 0, root's child is 1, etc).
    
    '''
    def __init__(self, *args, **kwargs):
        wx.TreeCtrl.__init__(self, *args, **kwargs)
        
    def getLevel(self, item):
        level = 0
        rootItem = self.GetRootItem()
        while item != rootItem:
            level += 1
            item = self.GetItemParent(item)
        return level
    
    def OnCompareItems(self, item1, item2):
        try:
            return cmp(int(self.GetItemText(item1)),
                       int(self.GetItemText(item2)))
        except:
            return cmp(self.GetItemText(item1), self.GetItemText(item2))

        
if __name__ == '__main__':
    app = wx.App()
    myFilter = filterGUI(None)
    testData = [
                ['uo2-30a_nov11a.spc', 1, 2, 3, 4, 5, 'ascan', '', 'Sun Jan 11 08:42:32 2011'],
                ['AllTheSamplesEverRunForever.spc', 2, 2, 3, 4, 5, 'ascan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 3, 2, 3, 4, 5, 'ascan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 4, 2, 3, 4, 5, 'ascan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 5, 2, 3, 4, 5, 'ascan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 6, 2, 3, 4, 5, 'ascan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 7, 2, 3, 4, 5, 'ascan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 8, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 9, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 10, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 11, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 12, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 13, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 14, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 15, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 16, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 17, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 18, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 19, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 20, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 21, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 22, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 23, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 24, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 25, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 26, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 27, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 28, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 29, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 30, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 31, 2, 3, 4, 5, 'rodscan', '', 'Today'],
                ['AllTheSamplesEverRunForever.spc', 32, 2, 3, 4, 5, 'rodscan', '', 'Today']
               ]
    for entry in testData:
        myFilter.dataTable.Append(entry)
    #myFilter.dataTable.SetHighlightFocusColour('Red')
    #myFilter.dataTable.SetItemBackgroundColour(2, wx.Colour(0, 200, 255))
    app.MainLoop()