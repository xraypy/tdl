"""
Scan Filter GUI
Author: Craig Biwer (cbiwer@uchicago.edu)
Last modified: 11.14.2011
"""
############################################################################

#Imports
from PythonCard import model
import wx
import os
import math
from wxUtil import wxUtil
from   specfile import SpecFile
import wx.lib.mixins.listctrl as listmix
from pds.lib.shellutil import mod_import

#Global variables, used to easily move parameters
# to and from the child filter windows

#All the scans read in from the specfile
allScans = []
#All sets of [H, K] values found in the file
# Form is: {[HK pair], [visible (T/F), [scan indexes], [distances]]}
allHKPairs = {}
#All the different scan types read in from the specfile
# Form is: {scan type, [visible (T/F), [scan indexes]]}
allTypes = {}
#List of all scan numbers to show, used when filtering by number
expandedRange = []
#List of all scans presently displayed; subset of allScans
outScans = []
#List of months present in the scan file, used when filtering by date
allMonths = []
#List of years present in the scan file, used when filtering by date
allYears = []
#Smallest L found in the specfile
fromL = 0.0
#Largest L found in the specfile
toL = 0.0
#Formatted strings corresponding to the date filter
fromDate = ""
toDate = ""

############################################################################

class wxScanFilter(model.Background, wxUtil):

    ###########################################################
    # Init and util methods
    ###########################################################
    def on_initialize(self, event):
        # Initialization
        
        #Set the GUI to display the current filename
        self.components.FileName.text = 'File name: ' + str(os.path.split(self.filename)[1])
        
        #Create a list of child integrator windows
        self.childCount = 0
        self.integratorChildren = {}
        
        #Parse the file and update the GUI
        self.parseFile(self.filename)
    
    ##############################################################
    
    #Using specfile.py (tdl\modules\xray\ana), parses the given
    # specfile, reading the header info for each scan and pulling
    # out the the relevant information. Updates the GUI at the end.
    # This only needs to be run once, though it is run again if the
    # user chooses to reset the filters and reread the file
    def parseFile(self, filename):
        
        #Declare global variables
        global allScans
        global allHKPairs
        global allTypes
        global outScans
        global expandedRange
        global allMonths
        global allYears
        global fromDate
        global toDate
        global fromL
        global toL
        
        #This reads the file and returns a list of dictionaries,
        # one per scan, with the relevant information
        # Information returned per scan:
        #   cmnd: Actual command issued
        #   index: Scan index number
        #   nl_start: Line number for the start of the scan
        #   date: The date and time of the scan
        #   time: How long the scan took
        #   G: All the G values
        #   Q: The Q value
        #   mot_names: The motor names
        #   P: The P values
        #   ncols: The number of columns in the data
        #   labels: The column headers for the scan data
        #   atten: The states of the attenuators
        #   energy: The energy used for the scan
        #   lineno: Line number for the end of the scan
        #   lStart: The actual starting L value
        #   lStop: The actual stopping L value
        #   aborted: Whether or not the scan was aborted
        #   dataStart: Line number for the start of the data
        #   dataStop: Line number for the end of the data
        specFile = SpecFile(filename)
        
        #Get set to parse the list of dictionaries, populating
        # allScansLocal while running then copying it to allScans
        # when complete. allYears keeps track of the years present
        # in the file, allMonths is a list of all the months in a
        # year, and missingMonths is used to keep track of which
        # months aren't present in the file (this is then used to
        # complement allMonths, which assures the months appear in
        # sorted order).
        allScansLocal = []
        #allDates = []
        allYears = []
        allMonths = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        missingMonths = allMonths[:]
        
        #Loops through the list of dictionaries, grabbing the desired
        # data and storing it for filtering
        for i in specFile._summary:
            scanIndex = str(i['index'])
            #hValD and kValD are the actual distances of H and K, respectively
            hValD = float(i['G'].split()[28])
            kValD = float(i['G'].split()[29])
            thetaValD = math.radians(180-float(i['G'].split()[33]))
            fullLine = i['cmd'].split()
            scanType = fullLine[0]
            #If it's a rodscan or hklscan, grab the H, K, and L values
            if scanType == 'rodscan':
                #hVal and kVal are the points in reciprocal space
                hVal = fullLine[1]
                kVal = fullLine[2]
                Li = i['lStart']
                Lf = i['lStop']
                #To get the distance to the point, we have to multiply the
                # reciprocal point by the distances of H and K
                # For example: if [H, K] is [2, -1] but H is 1.8 and K is 2,
                # we would need to use [3.6, -2] as our distances
                HKDist = round(((float(hVal)*hValD)**2 + (float(kVal)*kValD)**2 - 2*float(hVal)*hValD*float(kVal)*kValD*math.cos(thetaValD))**0.5, 8)
                if str([int(hVal), int(kVal)]) not in allHKPairs:
                    allHKPairs[str([int(hVal), int(kVal)])] = [True, [scanIndex], [HKDist]]
                else:
                    allHKPairs[str([int(hVal), int(kVal)])][1].append(scanIndex)
                    if HKDist not in allHKPairs[str([int(hVal), int(kVal)])][2]:
                        allHKPairs[str([int(hVal), int(kVal)])][2].append(HKDist)
                        allHKPairs[str([int(hVal), int(kVal)])][2].sort()
            elif scanType == 'hklscan':
                hVal = fullLine[1]
                kVal = fullLine[3]
                Li = fullLine[5]
                Lf = fullLine[6]
                HKDist = ((float(hVal)*hValD)**2 + (float(kVal)*kValD)**2 - 2*float(hVal)*hValD*float(kVal)*kValD*math.cos(thetaValD))**0.5
                if '[' + hVal + ', ' + kVal + ']' not in allHKPairs:
                    allHKPairs['[' + hVal + ', ' + kVal + ']'] = [True, [scanIndex], [HKDist]]
                else:
                    allHKPairs['[' + hVal + ', ' + kVal + ']'][1].append(scanIndex)
                    if HKDist not in allHKPairs['[' + hVal + ', ' + kVal + ']'][2]:
                        allHKPairs['[' + hVal + ', ' + kVal + ']'][2].append(HKDist)
                        allHKPairs['[' + hVal + ', ' + kVal + ']'][2].sort()
            else:
                hVal, kVal, HKDist, Li, Lf = ("--", "--", "--", "--", "--")
            if i['aborted']:
                aborted = 'Yes'
            else:
                aborted = ''
            scanDate = i['date']
            scanData = (i['cmd'], i['time'], i['atten'], i['energy'], i['nl_start'], i['lineno'], i['dataStart'], i['dataStop'], i['mot_names'])
            #Each loop, check that the fromL and toL
            # are still the min and max, respectively
            fromL = min(fromL, toL, Li, Lf)
            toL = max(fromL, toL, Li, Lf)
            
            allScansLocal.append([scanIndex, hVal, kVal, Li, Lf, scanType, aborted, scanDate, 'True', scanData, HKDist, i['G'], i['labels']])
            
            #Keep track of all the scan types
            if scanType not in allTypes:
                allTypes[scanType] = [True, [scanIndex]]
            else:
                allTypes[scanType][1].append(scanIndex)
        
            #Keep track of the months
            #allDates.append(scanDate.split())
            if scanDate.split()[1] in missingMonths:
                missingMonths.remove(scanDate.split()[1])
            if scanDate.split()[-1] not in allYears:
                allYears.append(scanDate.split()[-1])
            
        #With all scans parsed, sets the GUI to display the scans
        # and sets the column widths appropriately
        allScans = allScansLocal[:]
        outScans = allScansLocal[:]
        expandedRange = range(specFile.scan_min(), specFile.scan_max()+1)
        for thisMonth in missingMonths:
            allMonths.remove(thisMonth)
        allYears.sort()
        fromDate = 'Mon ' + allMonths[0] + ' 01 00:00:00 ' + allYears[0]
        toDate = 'Mon ' + allMonths[-1] + ' 31 23:59:00 ' + allYears[-1]
        self.components.unfilteredList.items = allScansLocal[:]
        self.components.unfilteredList.SetColumnWidth(0, 40) #Scan number
        self.components.unfilteredList.SetColumnWidth(1, 40) #H Val
        self.components.unfilteredList.SetColumnWidth(2, 38) #K Val
        self.components.unfilteredList.SetColumnWidth(3, 60) #L Start
        self.components.unfilteredList.SetColumnWidth(4, 60) #L Stop
        self.components.unfilteredList.SetColumnWidth(5, 66) #Scan Type
        self.components.unfilteredList.SetColumnWidth(6, 56) #Aborted
        #The 'Date' column is left to fill the remaining space
            
    #This function applies each filter to allScans, saving the results
    # to outScans. It then updates the GUI with the new output scans.
    def updateScans(self, event):
        global allScans
        global expandedRange
        global outScans
        #iterScans = allScans[:]
        outScans = allScans[:]
        
        for scan in allScans:
            #Determine whether the scan lies outside the L range
            #Set to False if not applicable
            try:
                self.lTest = float(scan[3]) <= fromL or float(scan[4]) >= toL
            except:
                self.lTest = False
            #Scan number
            if int(scan[0]) not in expandedRange:
                outScans.remove(scan)
            #HK Value
            elif scan[8] != 'True':
                outScans.remove(scan)
            #L Range
            elif self.lTest:
                outScans.remove(scan)
            #Scan Type
            elif not allTypes[scan[5]][0]:
                outScans.remove(scan)
            #Date Range
            elif self.isAfter(scan[7], str(toDate)) or self.isBefore(scan[7], str(fromDate)):
                outScans.remove(scan)
        #Update the GUI and reset the column widths (they resize on update)
        self.components.unfilteredList.items = outScans
        self.components.unfilteredList.SetColumnWidth(0, 40)
        self.components.unfilteredList.SetColumnWidth(1, 40)
        self.components.unfilteredList.SetColumnWidth(2, 38)
        self.components.unfilteredList.SetColumnWidth(3, 60)
        self.components.unfilteredList.SetColumnWidth(4, 60)
        self.components.unfilteredList.SetColumnWidth(5, 66)
        self.components.unfilteredList.SetColumnWidth(6, 56)
    
    #Given two dates, checks that the first comes before the second
    # Format: 'Day Month Date Time Year', i.e. 'Mon Jan 12 08:42:32 2011'
    def isBefore(self, amI, beforeMe):
        allMonths = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        (testDay, testMonth, testDate, testTime, testYear) = amI.split()
        (otherDay, otherMonth, otherDate, otherTime, otherYear) = beforeMe.split()
        if int(testYear) < int(otherYear):
            return True
        elif int(testYear) == int(otherYear):
            if allMonths.index(testMonth) < allMonths.index(otherMonth):
                return True
            elif allMonths.index(testMonth) == allMonths.index(otherMonth):
                if int(testDate) < int(otherDate):
                    return True
                elif int(testDate) == int(otherDate):
                    if testTime < otherTime:
                        return True
        return False
        
    #Given two dates, checks that the first comes after the second
    # Format: 'Day Month Date Time Year', i.e. 'Mon Jan 12 08:42:32 2011'
    def isAfter(self, amI, afterMe):
        return amI == afterMe or not self.isBefore(amI, afterMe)
    
    #Filter by date
    def on_DateFilter_mouseClick(self, event):
        dateFilters = self.DateWindow(self)
        self.Bind(wx.EVT_BUTTON, self.updateScans, id=dateFilters.apply.GetId())
    
    #Filter by [H, K] value
    def on_HKFilter_mouseClick(self, event):
        hkFilters = self.HKWindow(self)
        self.Bind(wx.EVT_BUTTON, self.updateScans, id=hkFilters.apply.GetId())
    
    #Filter by scan type
    def on_TypeFilter_mouseClick(self, event):
        typeFilters = self.TypeWindow(self)
        self.Bind(wx.EVT_BUTTON, self.updateScans, id=typeFilters.apply.GetId())
    
    #Keep the selected items, discard the rest
    def on_Update_mouseClick(self, event):
        numberOne = self.components.unfilteredList.GetFirstSelected()
        if numberOne == -1: return
        global outScans
        global expandedRange
        expandedRange = []
        while numberOne != -1:
            expandedRange.append(int(outScans[numberOne][0]))
            numberOne = self.components.unfilteredList.GetNextSelected(numberOne)
        self.updateScans(None)
    
    #Show more scan information when selected in the list
    def on_unfilteredList_select(self, event):
        scanNum = int(self.components.unfilteredList.getStringSelection()[0][0])
        fullScanData = allScans[scanNum - 1][9]
        #This is what's shown; more can be easily added
        toShow = 'Command: ' + fullScanData[0].strip() + '\n\n' + \
                    'Time: ' + fullScanData[1].strip() + '\n\n' + \
                    'Attenuators: ' + fullScanData[2].strip() + '\n\n' + \
                    'Energy: ' + fullScanData[3].strip() + '\n\n' + \
                    'Start line: ' + str(fullScanData[4]) + '\n\n' + \
                    'End line: ' + str(fullScanData[5]) + '\n\n' + \
                    'Data points: ' + str(int(fullScanData[7]) - int(fullScanData[6])) + '\n\n' + \
                    'HK Distance(s): ' + str(allScans[scanNum - 1][10]) + '\n\n'
        
        self.components.MoreInfoHere.text = toShow
    
    #Filter by scan number
    def on_NumberFilter_mouseClick(self, event):
        numberFilters = self.NumberWindow(self)
        self.Bind(wx.EVT_BUTTON, self.updateScans, id=numberFilters.apply.GetId())
    
    #Filter by L-range
    def on_FilterLs_mouseClick(self, event):
        lFilters = self.LWindow(self)
        self.Bind(wx.EVT_BUTTON, self.updateScans, id=lFilters.apply.GetId())
    
    #Reset all the global variables, then reread and reparse the file
    def on_ResetButton_mouseClick(self, event):
        global allScans
        global expandedRange
        global allHKPairs
        global allTypes
        global outScans
        global allMonths
        global allYears
        global fromL
        global toL
        global fromDate
        global toDate
        allScans = []
        expandedRange = []
        allHKPairs = {}
        allTypes = {}
        outScans = []
        allMonths = []
        allYears = []
        fromL = 0.0
        toL = 0.0
        fromDate = ""
        toDate = ""
        self.parseFile(self.filename)
    
    #Currently just prints out the filename, followed by the number, start
    # line, and ending line of each selected scan
    def on_SaveFiltered_mouseClick(self, event):
        global outScans
        #outStrings = []
        #for scan in outScans:
        #    outStrings.append("Scan " + scan[0] + " starts " + str(scan[9][4]) + \
        #                        " ends " + str(scan[9][5]))
        #toSave = str(os.path.split(self.filename)[1]) + '\n' + '\n'.join(outStrings)
        #print toSave
        
        import wxIntegrator
        self.childCount += 1
        integrator = mod_import(wxIntegrator)
        integrator = wxIntegrator.Integrator(self, self.filename, outScans)
        self.integratorChildren[self.childCount] = integrator
    
    #Format the integratorChildren keys to place in the AppendFiltered dialog
    def formatChoices(self, formatUs):
        formatedChoices = []
        for item in formatUs:
            formatedChoices.append("Specfile Integrator " + str(item))
        return formatedChoices
    
    #Pick an existing integrator child to which to append the currently selected scans
    def on_AppendFiltered_mouseClick(self,event):
        if len(self.integratorChildren) == 0:
            print 'No existing children.'
            return
        #for integrator, number
        appendOptions = self.formatChoices(self.integratorChildren.keys())
        chooseChild = wx.SingleChoiceDialog(self, 'Choose an existing integrator session:',\
                                                'Append to...',\
                                                appendOptions)
        if chooseChild.ShowModal() == wx.ID_OK:
            appendTo = chooseChild.GetStringSelection()[20:]
            appendTo = self.integratorChildren[int(appendTo)]
            appendTo.appendScans(outScans)
        chooseChild.Destroy()
    
    #When the filter is closed, close all its children as well
    def on_close(self,event):
        for key, child in self.integratorChildren.items():
            child.onClose(None)
        self.Destroy()

        
############################################################################

    class NumberWindow(wx.Frame):
        #This class is for filtering by scan number
        
        #Declare the global variables to be used
        global allScans
        global expandedRange
        global outScans
        def __init__(self, *args, **kwargs):
            #Make the window
            wx.Frame.__init__(self, args[0], -1, title="#s", size=(160, 800))
            self.panel = wx.Panel(self)
            self.list = self.NumberListCtrl(self.panel, style=wx.LC_REPORT | wx.LC_NO_HEADER)
            self.list.InsertColumn(0, "", width=20)
            self.list.InsertColumn(1, "Scan Number", width=80)
            
            #Populate the check box list with all the scan numbers,
            # checking the box if it's currently displayed
            for i in allScans:
                self.list.Append(["", i[0]])
                if i in outScans:
                    self.list.CheckItem(int(i[0]) - 1)
                    
            #Create the 'check all' check box (selector), the apply and cancel
            # buttons, and the input field for scan ranges
            self.selector = wx.CheckBox(self.panel, name="Selector", style=wx.CHK_3STATE, pos = (8, 8))
            self.apply = wx.Button(self.panel, label="Apply", name="Apply")
            self.cancel = wx.Button(self.panel, label="Cancel", name="Cancel")
            self.scanRange = wx.TextCtrl(self.panel, name="ScanNumbers", style=wx.TE_PROCESS_ENTER)
            self.scanRange.SetValue('Scan numbers')
            self.scanRange.SetToolTip(wx.ToolTip("Type scan numbers and ranges separated\nby commas, i.e. \'2, 4, 6-8\'"))
            
            #Set the 'check all' check box to the appropriate state (all, some, or none)
            self.setCheck(None)
            
            #Bindings (tie events to functions)
            self.Bind(wx.EVT_BUTTON, self.onApply, id=self.apply.GetId())
            self.Bind(wx.EVT_BUTTON, self.onCancel, id=self.cancel.GetId())
            self.selector.Bind(wx.EVT_CHECKBOX, self.onCheckAll)
            self.scanRange.Bind(wx.EVT_KEY_DOWN, self.onApplyRange)
            
            #Sizers (create the layout of the window)
            self.sizer = wx.BoxSizer(wx.VERTICAL)       #General layout
            self.sizer2 = wx.BoxSizer(wx.HORIZONTAL)    #The apply and cancel buttons
            self.sizer3 = wx.BoxSizer(wx.HORIZONTAL)    #The column headers
            
            #Top text box
            self.sizer.Add(self.scanRange, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #Column headers
            self.sizer3.Add(self.selector, flag=wx.EXPAND | wx.LEFT | wx.RIGHT, border = 4)
            self.sizer3.Add(wx.StaticText(self.panel, label="Scan Number"), flag=wx.EXPAND | wx.LEFT, border = 4)
            self.sizer.Add(self.sizer3, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #The list
            self.sizer.Add(self.list, proportion=1, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #The apply and cancel buttons
            self.sizer2.Add(self.apply, proportion=1, flag=wx.EXPAND | wx.ALL, border = 4)
            self.sizer2.Add(self.cancel, proportion=1, flag=wx.EXPAND | wx.ALL, border = 4)
            self.sizer.Add(self.sizer2, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #Create the window, center it, show it, bind the close button,
            # and lock the focus on it
            self.panel.SetSizerAndFit(self.sizer)
            self.CenterOnScreen()
            self.Show()
            self.Bind(wx.EVT_CLOSE, self.on_close)
            wx.Window.MakeModal(self)
            
        #Examine the state of all the check boxes until enough
        # is known to set the state of the 'check all' box
        def setCheck(self, event):
            num = self.list.GetItemCount()
            if num == 0: return
            thirdState = self.list.IsChecked(0)
            for i in range(1, num):
                if (thirdState == False and self.list.IsChecked(i)) or \
                    (thirdState == True and not self.list.IsChecked(i)):
                    self.selector.Set3StateValue(wx.CHK_UNDETERMINED)
                    return
                if self.list.IsChecked(i):
                    thirdState = True
                else:
                    thirdState = False
            self.selector.SetValue(thirdState)
        
        #When a user hits the enter (return) key from the text field, this
        # updates the check boxes in the list to copy the desired range
        def onApplyRange(self, event):
            global expandedRange
            keyCode = event.GetKeyCode()
            if keyCode == wx.WXK_RETURN or keyCode == 372:
                self.holdRange = self.readRange(self.scanRange.GetValue())
                if self.holdRange == None:
                    return
                expandedRange = self.holdRange
                for i in range(self.list.GetItemCount()):
                    self.list.CheckItem(i, i+1 in expandedRange)
                self.setCheck(None)
            else:
                event.Skip()
            return
        
        #Given a string of comma-delimited values and ranges, expands
        # the set to a list of single integer values
        # i.e. '2, 4, 6-8' becomes [2, 4, 6, 7, 8] 
        def readRange(self, inRange=''):
            if inRange == '' or inRange == None:
                return None
            rangeList = inRange.strip().replace(' ', '').split(',')
            fullList = []
            for thisPart in rangeList:
                try:
                    int(thisPart)
                    fullList.append(int(thisPart))
                except:
                    try:
                        miniRange = thisPart.split('-')
                        fullList.extend(range(int(miniRange[0]), int(miniRange[-1])+1))
                    except:
                        fullList = None
            return fullList
        
        #Select all or deselect all; if the 'check all' check box
        # is in the mixed state, deselects all
        def onCheckAll(self, event):
            if self.selector.Get3StateValue():
                self.onSelectAll(None)
            else:
                self.onDeselectAll(None)
        
        #Check all check boxes
        def onSelectAll(self, event):
            num = self.list.GetItemCount()
            for i in range(num):
                self.list.CheckItem(i)
    
        #Uncheck all check boxes
        def onDeselectAll(self, event):
            num = self.list.GetItemCount()
            for i in range(num):
                self.list.CheckItem(i, False)
        
        #First, release the focus so if something breaks you're
        # not stuck.
        # Then, get the indexes of the checked boxes, and keep them
        # in the  main window; hide the rest
        # Finally, skip the event so it can be caught by the main
        # window (allowing it to know when to update the list)
        def onApply(self, event):
            global expandedRange
            wx.Window.MakeModal(self, False)
            expandedRange = []
            num = self.list.GetItemCount()
            for i in range(num):
                if self.list.IsChecked(i):
                    expandedRange.append(i+1)
            self.Destroy()
            event.Skip()
        
        #Release the focus and close the window, making no changes
        def onCancel(self, event):
            wx.Window.MakeModal(self, False)
            self.Destroy()
        
        #Release the focus and close the window, making no changes
        def on_close(self,event):
            self.onCancel(event)
        
        #This class is used simply to make the check box list
        class NumberListCtrl(wx.ListCtrl, listmix.CheckListCtrlMixin, listmix.ListCtrlAutoWidthMixin):
            def __init__(self, *args, **kwargs):
                wx.ListCtrl.__init__(self, *args, **kwargs)
                listmix.CheckListCtrlMixin.__init__(self)
                listmix.ListCtrlAutoWidthMixin.__init__(self)
            
            #This function can be used to update the 'check all' check box,
            # but for specfiles with large numbers of scans it is noticeably slow
            #def OnCheckItem(self, index, flag):
            #    wx.PostEvent(self.GetEventHandler(), wx.PyCommandEvent(wx.EVT_CHECKBOX.typeId, self.GetId()))
                
#######################################################################

    class HKWindow(wx.Frame):
        #This class is for filtering by HK pair
        global allHKPairs
        def __init__(self, *args, **kwargs):
            #Make the window
            wx.Frame.__init__(self, args[0], -1, title="HK Pairs", size=(400, 400))
            self.panel = wx.Panel(self)
            self.list = self.HKListCtrl(self.panel, style=wx.LC_REPORT | wx.LC_NO_HEADER | wx.LC_HRULES | wx.LC_VRULES)
            self.list.InsertColumn(0, "", width=24)
            self.list.InsertColumn(1, "HK Pair", width=64)
            self.list.InsertColumn(2, "Num. Scans", width = 64)
            self.list.InsertColumn(3, "HK Distance(s)")
            self.list.Arrange()
            self.counter = 0
            self.keyToIndex = {}
            #Sort the HK pairs first by H and K, then by distance,
            # resulting in distances grouped together and sorted
            # internally by H and K
            allHKKeys = allHKPairs.keys()
            allHKKeys = sorted(allHKKeys, key = lambda key: eval(key))
            allHKKeys = sorted(allHKKeys, key = lambda key: allHKPairs[key][2])
            for key in allHKKeys:
                self.list.Append(["", key, len(allHKPairs[key][1]), allHKPairs[key][2]])
                self.keyToIndex[key] = self.counter
                self.list.CheckItem(self.keyToIndex[key], allHKPairs[key][0])
                self.counter += 1
            
            #Creates the 'check all' check box ('Selector'), as well as the
            # apply and cancel buttons
            self.selector = wx.CheckBox(self.panel, name="Selector", style=wx.CHK_3STATE, pos = (8, 8))
            self.apply = wx.Button(self.panel, label="Apply", name="Apply")
            self.cancel = wx.Button(self.panel, label="Cancel", name="Cancel")
            #Set the appropriate state for the 'check all' box
            self.setCheck(None)
            
            #Bindings
            self.Bind(wx.EVT_BUTTON, self.onApply, id=self.apply.GetId())
            self.Bind(wx.EVT_BUTTON, self.onCancel, id=self.cancel.GetId())
            self.Bind(wx.EVT_CHECKBOX, self.setCheck)
            self.selector.Bind(wx.EVT_CHECKBOX, self.onCheckAll)
        
            #Sizers
            self.sizer = wx.BoxSizer(wx.VERTICAL)       #General layout
            self.sizer2 = wx.BoxSizer(wx.HORIZONTAL)    #The apply and cancel buttons
            self.sizer3 = wx.BoxSizer(wx.HORIZONTAL)    #The column headings
            
            #The column headings (including the 'check all' box)
            self.sizer3.Add(self.selector, flag=wx.EXPAND | wx.LEFT | wx.RIGHT, border = 4)
            self.sizer3.Add(wx.StaticText(self.panel, label="HK Pair"), flag=wx.EXPAND | wx.LEFT, border = 12)
            self.sizer3.Add(wx.StaticText(self.panel, label="Num. Scans"), flag=wx.EXPAND | wx.LEFT, border = 24)
            self.sizer3.Add(wx.StaticText(self.panel, label="HK Distance(s)"), flag=wx.EXPAND | wx.LEFT, border = 16)
            self.sizer.Add(self.sizer3, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #The list
            self.sizer.Add(self.list, proportion=1, flag=wx.EXPAND | wx.ALL, border=4)
        
            #The apply and cancel buttons
            self.sizer2.Add(self.apply, proportion=1, flag=wx.EXPAND | wx.ALL, border = 4)
            self.sizer2.Add(self.cancel, proportion=1, flag=wx.EXPAND | wx.ALL, border = 4)
            self.sizer.Add(self.sizer2, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #Create the window, center it, show it, bind the close button,
            # and lock the focus on it
            self.panel.SetSizerAndFit(self.sizer)
            self.CenterOnScreen()
            self.Show()
            self.Bind(wx.EVT_CLOSE, self.on_close)
            wx.Window.MakeModal(self)
        
        #Examine the state of all the check boxes until enough
        # is known to set the state of the 'check all' box
        def setCheck(self, event):
            num = self.list.GetItemCount()
            if num == 0: return
            thirdState = self.list.IsChecked(0)
            for i in range(1, num):
                if (thirdState == False and self.list.IsChecked(i)) or \
                    (thirdState == True and not self.list.IsChecked(i)):
                    self.selector.Set3StateValue(wx.CHK_UNDETERMINED)
                    return
                if self.list.IsChecked(i):
                    thirdState = True
                else:
                    thirdState = False
            self.selector.SetValue(thirdState)
        
        #Select all or deselect all; if the 'check all' check box
        # is in the mixed state, deselects all
        def onCheckAll(self, event):
            if self.selector.Get3StateValue():
                self.onSelectAll(None)
            else:
                self.onDeselectAll(None)
            
        #Check all check boxes
        def onSelectAll(self, event):
            num = self.list.GetItemCount()
            for i in range(num):
                self.list.CheckItem(i)
        
        #Uncheck all check boxes
        def onDeselectAll(self, event):
            num = self.list.GetItemCount()
            for i in range(num):
                self.list.CheckItem(i, False)
    
        #First, release the focus so if something breaks you're
        # not stuck.
        # Then, for each key in the list, examine the state of its check
        # box and apply it to the corresponding element in allHKPairs
        # Finally, skip the event so it can be caught by the main
        # window (allowing it to know when to update the list)
        def onApply(self, event):
            wx.Window.MakeModal(self, False)
            for key in allHKPairs:
                allHKPairs[key][0] = self.list.IsChecked(self.keyToIndex[key])
                trueFalse, iterList, distList = allHKPairs[key]
                for i in iterList:
                    intI = int(i) - 1
                    allScans[intI][8] = str(trueFalse)
            self.Destroy()
            event.Skip()
        
        #Release the focus and close the window, making no changes
        def onCancel(self, event):
            wx.Window.MakeModal(self, False)
            self.Destroy()
        
        #Release the focus and close the window, making no changes
        def on_close(self,event):
            self.onCancel(event)
        
        #This class is used simply to make the check box list
        class HKListCtrl(wx.ListCtrl, listmix.CheckListCtrlMixin, listmix.ListCtrlAutoWidthMixin):
            def __init__(self, *args, **kwargs):
                wx.ListCtrl.__init__(self, *args, **kwargs)
                listmix.CheckListCtrlMixin.__init__(self)
                listmix.ListCtrlAutoWidthMixin.__init__(self)
                
            #This function can be used to update the 'check all' check box,
            # but for specfiles with large numbers of different HK pairs
            # it is noticeably slow
            def OnCheckItem(self, index, flag):
                wx.PostEvent(self.GetEventHandler(), wx.PyCommandEvent(wx.EVT_CHECKBOX.typeId, self.GetId()))
            
###########################################################################
            
    class LWindow(wx.Frame):
        #This class is for filtering by L value
        global fromL
        global toL
        def __init__(self, *args, **kwargs):
            #Make the window
            wx.Frame.__init__(self, args[0], -1, title="Scan Types", size=(232, 108))
            self.panel = wx.Panel(self)
            
            #'From' value components
            self.fromText = wx.StaticText(self.panel, label="From:")
            self.fromValue = wx.TextCtrl(self.panel, name="fromValue", size=(60, -1))
            self.fromValue.SetValue(str(fromL))
            
            #'To' value components
            self.toText = wx.StaticText(self.panel, label="To:")
            self.toValue = wx.TextCtrl(self.panel, name="toValue", size=(60, -1))
            self.toValue.SetValue(str(toL))
            
            #Buttons (apply and cancel)
            self.apply = wx.Button(self.panel, label="Apply", name="Apply")
            self.cancel = wx.Button(self.panel, label="Cancel", name="Cancel")
            
            #Bindings
            self.Bind(wx.EVT_BUTTON, self.onApply, id=self.apply.GetId())
            self.Bind(wx.EVT_BUTTON, self.onCancel, id=self.cancel.GetId())
            
            #Sizers
            self.sizer = wx.BoxSizer(wx.VERTICAL)       #General layout
            self.sizer2 = wx.BoxSizer(wx.HORIZONTAL)    #'From' and 'To' layout
            self.sizer4 = wx.BoxSizer(wx.HORIZONTAL)    #Button layout
            
            #'From' input followed by 'To' input
            self.sizer2.Add(self.fromText, flag=wx.EXPAND | wx.TOP | wx.BOTTOM | wx.LEFT, border = 8)
            self.sizer2.Add(self.fromValue, flag=wx.EXPAND | wx.ALL, border = 4)
            self.sizer2.Add(self.toText, flag=wx.EXPAND | wx.TOP | wx.BOTTOM | wx.LEFT, border = 8)
            self.sizer2.Add(self.toValue, flag=wx.EXPAND | wx.ALL, border = 4)
            self.sizer.Add(self.sizer2, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #Button positions
            self.sizer4.Add(self.apply, proportion=1, flag=wx.EXPAND | wx.RIGHT | wx.BOTTOM | wx.LEFT, border = 4)
            self.sizer4.Add(self.cancel, proportion=1, flag=wx.EXPAND | wx.RIGHT | wx.BOTTOM | wx.LEFT, border = 4)
            self.sizer.Add(self.sizer4, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #Create the window, center it, show it, bind the close button,
            # and lock the focus on it
            self.panel.SetSizerAndFit(self.sizer)
            self.CenterOnScreen()
            self.Show()
            self.Bind(wx.EVT_CLOSE, self.on_close)
            wx.Window.MakeModal(self)
    
        #First, release the focus so if something breaks you're
        # not stuck.
        # Then, update the global fromL and toL values so the
        # parent class can access the changes
        # Finally, skip the event so it can be caught by the main
        # window (allowing it to know when to update the list)
        def onApply(self, event):
            wx.Window.MakeModal(self, False)
            global fromL
            global toL
            fromL = float(self.fromValue.GetValue())
            toL = float(self.toValue.GetValue())
            self.Destroy()
            event.Skip()
        
        #Release the focus and close the window, making no changes
        def onCancel(self, event):
            wx.Window.MakeModal(self, False)
            self.Destroy()
        
        #Release the focus and close the window, making no changes
        def on_close(self,event):
            self.onCancel(event)
            
###############################################################################
        
    class TypeWindow(wx.Frame):
        #This class is for filtering by scan type
        global allTypes
        def __init__(self, *args, **kwargs):
            #Make the window
            wx.Frame.__init__(self, args[0], -1, title="Scan Types", size=(200, 400))
            self.panel = wx.Panel(self)
            self.list = self.TypeListCtrl(self.panel, style=wx.LC_REPORT | wx.LC_NO_HEADER | wx.LC_HRULES | wx.LC_VRULES)
            self.list.InsertColumn(0, "", width=24)
            self.list.InsertColumn(1, "Scan Type", width=66)
            self.list.InsertColumn(2, "Num. Scans")
            
            #Populate the list with the different scan types; check
            # the appropriate boxes to show currently selected types
            self.counter = 0
            self.typeToIndex = {}
            for type in allTypes:
                self.list.Append(["", type, len(allTypes[type][1])])
                self.typeToIndex[type] = self.counter
                self.list.CheckItem(self.typeToIndex[type], allTypes[type][0])
                self.counter += 1
            
            #Create the 'check all' check box (selector), the apply and cancel
            # buttons, and the input field for scan ranges
            self.selector = wx.CheckBox(self.panel, name="Selector", style=wx.CHK_3STATE, pos = (8, 8))
            self.apply = wx.Button(self.panel, label="Apply", name="Apply")
            self.cancel = wx.Button(self.panel, label="Cancel", name="Cancel")
            
            #Set the 'check all' check box to the appropriate state (all, some, or none)
            self.setCheck(None)
            
            #Bindings
            self.Bind(wx.EVT_BUTTON, self.onApply, id=self.apply.GetId())
            self.Bind(wx.EVT_BUTTON, self.onCancel, id=self.cancel.GetId())
            self.Bind(wx.EVT_CHECKBOX, self.setCheck)
            self.selector.Bind(wx.EVT_CHECKBOX, self.onCheckAll)
            
            #Sizers
            self.sizer = wx.BoxSizer(wx.VERTICAL)       #General layout
            self.sizer2 = wx.BoxSizer(wx.HORIZONTAL)    #The apply and cancel buttons
            self.sizer3 = wx.BoxSizer(wx.HORIZONTAL)    #The column headers
            
            #The column headers
            self.sizer3.Add(self.selector, flag=wx.EXPAND | wx.LEFT | wx.RIGHT, border = 4)
            self.sizer3.Add(wx.StaticText(self.panel, label="Scan Type"), flag=wx.EXPAND | wx.LEFT, border = 12)
            self.sizer3.Add(wx.StaticText(self.panel, label="Num. Scans"), flag=wx.EXPAND | wx.LEFT, border = 12)
            self.sizer.Add(self.sizer3, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #The list
            self.sizer.Add(self.list, proportion=1, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #The apply and cancel  buttons
            self.sizer2.Add(self.apply, proportion=1, flag=wx.EXPAND | wx.ALL, border = 4)
            self.sizer2.Add(self.cancel, proportion=1, flag=wx.EXPAND | wx.ALL, border = 4)
            self.sizer.Add(self.sizer2, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #Create the window, center it, show it, bind the close button,
            # and lock the focus on it
            self.panel.SetSizerAndFit(self.sizer)
            self.CenterOnScreen()
            self.Show()
            self.Bind(wx.EVT_CLOSE, self.on_close)
            wx.Window.MakeModal(self)
        
        #Examine the state of all the check boxes until enough
        # is known to set the state of the 'check all' box
        def setCheck(self, event):
            num = self.list.GetItemCount()
            if num == 0: return
            thirdState = self.list.IsChecked(0)
            for i in range(1, num):
                if (thirdState == False and self.list.IsChecked(i)) or \
                    (thirdState == True and not self.list.IsChecked(i)):
                    self.selector.Set3StateValue(wx.CHK_UNDETERMINED)
                    return
                if self.list.IsChecked(i):
                    thirdState = True
                else:
                    thirdState = False
            self.selector.SetValue(thirdState)
        
        #Select all or deselect all; if the 'check all' check box
        # is in the mixed state, deselects all
        def onCheckAll(self, event):
            if self.selector.Get3StateValue():
                self.onSelectAll(None)
            else:
                self.onDeselectAll(None)
        
        #Check all check boxes
        def onSelectAll(self, event):
            num = self.list.GetItemCount()
            for i in range(num):
                self.list.CheckItem(i)
    
        #Uncheck all check boxes
        def onDeselectAll(self, event):
            num = self.list.GetItemCount()
            for i in range(num):
                self.list.CheckItem(i, False)
    
        #First, release the focus so if something breaks you're
        # not stuck.
        # Then, update the global allTypes so the parent class can
        # tell which scans to show and which to hide.
        # Finally, skip the event so it can be caught by the main
        # window (allowing it to know when to update the list)
        def onApply(self, event):
            wx.Window.MakeModal(self, False)
            for type in allTypes:
                allTypes[type][0] = self.list.IsChecked(self.typeToIndex[type])
            self.Destroy()
            event.Skip()
        
        #Release the focus and close the window, making no changes
        def onCancel(self, event):
            wx.Window.MakeModal(self, False)
            self.Destroy()
        
        #Release the focus and close the window, making no changes
        def on_close(self,event):
            self.onCancel(event)
        
        #This class is used simply to make the check box list
        class TypeListCtrl(wx.ListCtrl, listmix.CheckListCtrlMixin, listmix.ListCtrlAutoWidthMixin):
            def __init__(self, *args, **kwargs):
                wx.ListCtrl.__init__(self, *args, **kwargs)
                listmix.CheckListCtrlMixin.__init__(self)
                listmix.ListCtrlAutoWidthMixin.__init__(self)
            
            #This function can be used to update the 'check all' check box,
            # but for specfiles with large numbers of scans it is noticeably slow
            def OnCheckItem(self, index, flag):
                wx.PostEvent(self.GetEventHandler(), wx.PyCommandEvent(wx.EVT_CHECKBOX.typeId, self.GetId()))
                
######################################################################

    class DateWindow(wx.Frame):
        #This class is for filtering by date
        global allMonths
        global allYears
        global fromDate
        global toDate
        def __init__(self, *args, **kwargs):
            #Make the window
            wx.Frame.__init__(self, args[0], -1, title="Scan Types", size=(288, 132))
            self.panel = wx.Panel(self)
            
            #Set the default choices for day, hour, and minute
            # Note the user would be allowed to select Febuary 31 since the
            # choices are static, but it won't break the program
            # i.e. February 31 comes after February 30 but before March 1
            self.dayChoices = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', \
                            '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', \
                            '25', '26', '27', '28', '29', '30', '31']
            self.hourChoices = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', \
                            '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23']
            self.minuteChoices = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', \
                            '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', \
                            '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', \
                            '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', \
                            '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59']
            
            #Create the colons for between the hour and minute choices, set to a
            # larger font so they're easier to see
            self.colonText = wx.StaticText(self.panel, label=":")
            self.colonText.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
            self.colonText2 = wx.StaticText(self.panel, label=":")
            self.colonText2.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
            
            #Buttons
            self.apply = wx.Button(self.panel, label="Apply", name="Apply")
            self.cancel = wx.Button(self.panel, label="Cancel", name="Cancel")
            
            #Bindings
            self.Bind(wx.EVT_BUTTON, self.onApply, id=self.apply.GetId())
            self.Bind(wx.EVT_BUTTON, self.onCancel, id=self.cancel.GetId())

            #Create the 'From' fields
            self.fromText = wx.StaticText(self.panel, label="From:")
            self.fromMonth = wx.Choice(self.panel, name="fromMonth", choices = allMonths, size=(46, 21))
            self.fromDay = wx.Choice(self.panel, name="fromDay", choices = self.dayChoices, size=(36, 21))
            self.fromYear = wx.Choice(self.panel, name="fromYear", choices = allYears, size = (54, 21))
            self.fromHour = wx.Choice(self.panel, name="fromHour", choices = self.hourChoices, size=(36, 21))
            self.fromMinute = wx.Choice(self.panel, name="fromMinute", choices = self.minuteChoices, size=(36, 21))
            
            #Set the initial selections for the 'From' fields
            self.fromMonth.SetStringSelection(fromDate.split()[1])
            self.fromDay.SetStringSelection(fromDate.split()[2])
            self.fromYear.SetStringSelection(fromDate.split()[4])
            self.fromHour.SetStringSelection(fromDate.split()[3].split(':')[0])
            self.fromMinute.SetStringSelection(fromDate.split()[3].split(':')[1])
            
            #Create the 'To' fields
            self.toText = wx.StaticText(self.panel, label="To:")
            self.toMonth = wx.Choice(self.panel, name="toMonth", choices = allMonths, size=(46, 21))
            self.toDay = wx.Choice(self.panel, name="toDay", choices = self.dayChoices, size=(36, 21))
            self.toYear = wx.Choice(self.panel, name="toYear", choices = allYears, size = (54, 21))
            self.toHour = wx.Choice(self.panel, name="toHour", choices = self.hourChoices, size=(36, 21))
            self.toMinute = wx.Choice(self.panel, name="toMinute", choices = self.minuteChoices, size=(36, 21))
            
            #Set the initial selections for the 'To' fields
            self.toMonth.SetStringSelection(toDate.split()[1])
            self.toDay.SetStringSelection(toDate.split()[2])
            self.toYear.SetStringSelection(toDate.split()[4])
            self.toHour.SetStringSelection(toDate.split()[3].split(':')[0])
            self.toMinute.SetStringSelection(toDate.split()[3].split(':')[1])
            
            #Sizers
            self.sizer = wx.BoxSizer(wx.VERTICAL)       #General layout
            self.sizer2 = wx.BoxSizer(wx.HORIZONTAL)    #'From' values
            self.sizer3 = wx.BoxSizer(wx.HORIZONTAL)    #'To' values
            self.sizer4 = wx.BoxSizer(wx.HORIZONTAL)    #Apply and cancel buttons
            
            #'From' sizer
            self.sizer2.Add(self.fromText, flag=wx.EXPAND | wx.TOP | wx.BOTTOM | wx.LEFT, border = 4)
            self.sizer2.Add(self.fromMonth, flag=wx.EXPAND | wx.LEFT, border = 2)
            self.sizer2.Add(self.fromDay, flag=wx.EXPAND | wx.LEFT, border = 2)
            self.sizer2.Add(self.fromYear, flag=wx.EXPAND | wx.LEFT, border = 2)
            self.sizer2.Add(self.fromHour, flag=wx.EXPAND | wx.LEFT, border = 8)
            self.sizer2.Add(self.colonText, flag=wx.EXPAND | wx.LEFT, border = 0)
            self.sizer2.Add(self.fromMinute, flag=wx.EXPAND | wx.LEFT, border = 0)
            self.sizer.Add(self.sizer2, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #'To' sizer
            self.sizer3.Add(self.toText, flag=wx.EXPAND | wx.TOP | wx.BOTTOM | wx.LEFT, border = 4)
            self.sizer3.Add(self.toMonth, flag=wx.EXPAND | wx.LEFT, border = 14)
            self.sizer3.Add(self.toDay, flag=wx.EXPAND | wx.LEFT, border = 2)
            self.sizer3.Add(self.toYear, flag=wx.EXPAND | wx.LEFT, border = 2)
            self.sizer3.Add(self.toHour, flag=wx.EXPAND | wx.LEFT, border = 8)
            self.sizer3.Add(self.colonText2, flag=wx.EXPAND | wx.LEFT, border = 0)
            self.sizer3.Add(self.toMinute, flag=wx.EXPAND | wx.LEFT, border = 0)
            self.sizer.Add(self.sizer3, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #Apply and cancel buttons
            self.sizer4.Add(self.apply, proportion=1, flag=wx.EXPAND | wx.ALL, border = 4)
            self.sizer4.Add(self.cancel, proportion=1, flag=wx.EXPAND | wx.ALL, border = 4)
            self.sizer.Add(self.sizer4, flag=wx.EXPAND | wx.ALL, border = 4)
            
            #Create the window, center it, show it, bind the close button,
            # and lock the focus on it
            self.panel.SetSizerAndFit(self.sizer)
            self.CenterOnScreen()
            self.Show()
            self.Bind(wx.EVT_CLOSE, self.on_close)
            wx.Window.MakeModal(self)
    
        #First, release the focus so if something breaks you're
        # not stuck.
        # Then, update the global fromDate and toDate values so the
        # parent class can access the changes
        # Finally, skip the event so it can be caught by the main
        # window (allowing it to know when to update the list)
        def onApply(self, event):
            global fromDate
            global toDate
            wx.Window.MakeModal(self, False)
            fromDate = 'Mon ' + self.fromMonth.GetStringSelection() + ' ' + \
                        self.fromDay.GetStringSelection() + ' ' + self.fromHour.GetStringSelection() + ':' + \
                        self.fromMinute.GetStringSelection() + ':00 ' + self.fromYear.GetStringSelection()
            toDate = 'Mon ' + self.toMonth.GetStringSelection() + ' ' + \
                        self.toDay.GetStringSelection() + ' ' + self.toHour.GetStringSelection() + ':' + \
                        self.toMinute.GetStringSelection() + ':00 ' + self.toYear.GetStringSelection()
            self.Destroy()
            event.Skip()
        
        #Release the focus and close the window, making no changes
        def onCancel(self, event):
            wx.Window.MakeModal(self, False)
            self.Destroy()
        
        #Release the focus and close the window, making no changes
        def on_close(self,event):
            self.onCancel(event)
            
    