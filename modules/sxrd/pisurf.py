"""
Functions and classes used to build the pi-surf GUI

Authors/modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""
###############################################################################

import wx
import os
import time
import numpy as Num

from tdl.modules.sxrd.ctrfitcalcs import *
from tdl.modules.sxrd.simplex import simplex
from tdl.modules.sxrd.genetic_algorithm import genetic
from tdl.modules.sxrd.pisurf_resonant import ResonantDataPanel
from tdl.modules.sxrd.pisurf_resonant import read_RSD, read_f1f2 
from tdl.modules.sxrd.pisurf_resonant import RASD_Fourier, write_rids
from tdl.modules.sxrd.fourierframe import createfourierframe

from tdl.modules.xtab.atomic import f0data as database
###############################################################################
class wxCtrFitFrame(wx.Frame):
    def __init__(self, parent, title, size):
        wx.Frame.__init__(self, parent, title= title, size=size)

        self.filename = ['','','','','','']
        self.dirname = ['','','','','','']
        self.filename7 = ''
        # A status bar
        self.statusbar = self.CreateStatusBar(2)

        # Setting up the ReadFilesmenu.
        filemenu= wx.Menu()
        menuReadData = filemenu.Append(wx.ID_ANY,"&Read Datafile"," Read in a Data file")
        menuReadBulk = filemenu.Append(wx.ID_ANY,"&Read Bulkfile"," Read in a Bulk file")
        menuReadSurface = filemenu.Append(wx.ID_ANY,"&Read Surfacefile"," Read in a Surface file")
        menuReadParameter = filemenu.Append(wx.ID_ANY,"&Read Parameterfile"," Read in a Parameter file")
        menuReadRigidbody = filemenu.Append(wx.ID_ANY,"&Read Rigidbodyfile"," Read in a Rigid Body file")
        menuReadBV = filemenu.Append(wx.ID_ANY,"&Read Bondvalencefile"," Read in a Bond Valence file")
        filemenu.AppendSeparator()
        menuLoadSession = filemenu.Append(wx.ID_ANY,"&Load Session"," Read in Files previously saved in a Session file") 
        menuExit = filemenu.Append(wx.ID_EXIT,"E&xit"," Terminate the program")

        #Setting up the WriteFiles
        writemenu = wx.Menu()
        menuWriteCif = writemenu.Append(wx.ID_ANY, "&Write .cif file", " Write surface structure to a .cif file")
        menuWritePar = writemenu.Append(wx.ID_ANY, "&Write Parameter file", " Write parameters to a .par file")
        menuWriteData = writemenu.Append(wx.ID_ANY, "&Write Data file", " Write data and fit results to a .dat file")
        menuWriteEdens = writemenu.Append(wx.ID_ANY, "&Write Edensity file", " Write E-density plot curves to a file")
        menuWriteSurf = writemenu.Append(wx.ID_ANY, "&Write Surface file", " Write surface in fractional coordinates to a .sur file")
        writemenu.AppendSeparator()
        menuSaveSession = writemenu.Append(wx.ID_ANY,"&Save Session"," Save paths to currently loaded files to file for fast restart")
        #menuWriteBulk = writemenu.Append(wx.ID_ANY, "&Write Bulk file", " Write bulk structure to a .bul file")

        #Setting up the resonant files menu
        ridsmenu = wx.Menu()
        menuReadf1f2 = ridsmenu.Append(wx.ID_ANY, "&Read .f1f2 file", " Read in f1f2 Data ")
        menuReadridsdata = ridsmenu.Append(wx.ID_ANY, "&Read resonant data files", " Specify the file containing the names of all the .rsd files to be read ")
        ridsmenu.AppendSeparator()
        menuWriteRids = ridsmenu.Append(wx.ID_ANY, "&Write RIDS data file", " write RIDS data and fit results to a file ")
        
        # Creating the menubar.
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&Read Files") # Adding the "filemenu" to the MenuBar
        menuBar.Append(writemenu,"&Write Files")
        menuBar.Append(ridsmenu,"&RIDS Files")
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        
        # define events
        self.Bind(wx.EVT_MENU, self.OnReadData, menuReadData)
        self.Bind(wx.EVT_MENU, self.OnReadBulk, menuReadBulk)
        self.Bind(wx.EVT_MENU, self.OnReadSurface, menuReadSurface)
        self.Bind(wx.EVT_MENU, self.OnReadParameter, menuReadParameter)
        self.Bind(wx.EVT_MENU, self.OnReadRigidbody, menuReadRigidbody)
        self.Bind(wx.EVT_MENU, self.OnReadBV, menuReadBV)
        self.Bind(wx.EVT_MENU, self.LoadSession, menuLoadSession)
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)

        self.Bind(wx.EVT_MENU, self.OnWriteCif, menuWriteCif)
        self.Bind(wx.EVT_MENU, self.OnWritePar, menuWritePar)
        self.Bind(wx.EVT_MENU, self.OnWriteData, menuWriteData)
        self.Bind(wx.EVT_MENU, self.OnWriteEdens, menuWriteEdens)
        self.Bind(wx.EVT_MENU, self.OnWriteSurf, menuWriteSurf)
        self.Bind(wx.EVT_MENU, self.SaveSession, menuSaveSession)

        self.Bind(wx.EVT_MENU, self.OnReadf1f2, menuReadf1f2)
        self.Bind(wx.EVT_MENU, self.OnReadRids, menuReadridsdata)
        self.Bind(wx.EVT_MENU, self.OnWriteRids, menuWriteRids)

        self.nb = CtrNotebook(self)

        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.Show(True)
        

    def OnReadData(self,e):
        dlg = wx.FileDialog(self, "Choose a data file ", self.dirname[0], ".dat", "*.dat", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename[0] = dlg.GetFilename()
            self.dirname[0] = dlg.GetDirectory()
            self.nb.data = read_data(self.dirname[0]+'/'+self.filename[0])
            self.nb.MainControlPage.Rdata = calc_Rdata(self.nb.data)
            self.nb.MainControlPage.datafile.SetValue(self.filename[0])
            self.nb.MainControlPage.Rod_weight = []
            self.nb.MainControlPage.rodweight = []
            self.nb.MainControlPage.rmsflaglog.SetValue(self.nb.MainControlPage.RMS_flag_log[self.nb.MainControlPage.RMS_flag]+self.nb.MainControlPage.Rdata[self.nb.MainControlPage.RMS_flag])
            if self.nb.bulk != []:
                for x in self.nb.data:
                    x.calcFbulk(self.nb.cell, self.nb.bulk, self.nb.g_inv, database)
   
            for i in range(len(self.nb.data)):
                wx.StaticText(self.nb.MainControlPage, label = (str(int(self.nb.data[i].H))+' '+str(int(self.nb.data[i].K))+' L'), pos=(350,25*i+67), size=(40,20))
                self.nb.MainControlPage.rodweight.append(wx.TextCtrl(self.nb.MainControlPage,1000+i, pos=(400,25*i+65), size=(30,20)))
                self.nb.MainControlPage.rodweight[i].SetValue('1')
                self.nb.MainControlPage.Rod_weight.append(1)
                self.Bind(wx.EVT_TEXT, self.nb.MainControlPage.setrodweight, self.nb.MainControlPage.rodweight[i])
            self.nb.SetSelection(0)
        dlg.Destroy()

    def OnReadBulk(self,e):
        """ Read in bulk file"""
        dlg = wx.FileDialog(self, "Choose a bulk file", self.dirname[1], ".bul", "*.bul", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename[1] = dlg.GetFilename()
            self.dirname[1] = dlg.GetDirectory()
            self.nb.bulk, self.nb.cell, self.nb.NLayers = read_bulk(self.dirname[1]+'/'+self.filename[1])
            self.nb.MainControlPage.bulkfile.SetValue(self.filename[1])
            self.nb.MainControlPage.bulkfile.AppendText('\na = '+str(round(self.nb.cell[0],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\nb = '+str(round(self.nb.cell[1],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\nc = '+str(round(self.nb.cell[2],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\nalpha = '+str(round(self.nb.cell[3],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\nbeta = '+str(round(self.nb.cell[4],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\ngamma = '+str(round(self.nb.cell[5],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\ndelta1 = '+str(round(self.nb.cell[6],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\ndelta2 = '+str(round(self.nb.cell[7],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\nNLayers = '+str(self.nb.NLayers))
            self.nb.g_inv = calc_g_inv(self.nb.cell)
            self.nb.ResonantDataPage.allrasd.cell = self.nb.cell
            self.nb.ResonantDataPage.allrasd.g_inv = self.nb.g_inv
            if self.nb.data != []:
                for x in self.nb.data:
                    x.calcFbulk(self.nb.cell, self.nb.bulk, self.nb.g_inv, database)                
            self.nb.SetSelection(0)
        dlg.Destroy()

    def OnReadSurface(self,e):
        """ Read in surface file"""
        dlg = wx.FileDialog(self, "Choose a surface file", self.dirname[2], ".sur", "*.sur", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename[2] = dlg.GetFilename()
            self.dirname[2] = dlg.GetDirectory()
            self.nb.surface, self.nb.parameter_usage, self.nb.runningDB = read_surface(self.dirname[2]+'/'+self.filename[2], database)
            self.nb.MainControlPage.surfacefile.SetValue(self.filename[2])
            self.nb.SetSelection(0)
        dlg.Destroy()

    def OnReadParameter(self,e):
        """ Read in parameter file"""
        dlg = wx.FileDialog(self, "Choose a parameter file", self.dirname[3], ".par", "*.par", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename[3] = dlg.GetFilename()
            self.dirname[3] = dlg.GetDirectory()
            self.nb.parameter, self.nb.param_labels = read_parameters(self.dirname[3]+'/'+self.filename[3])
            self.nb.MainControlPage.parameterfile.SetValue(self.filename[3])
           
            self.nb.DeletePage(2)
            self.nb.ParameterPage = ParameterPanel(self.nb)
            self.nb.AddPage(self.nb.ParameterPage, " Parameters ")
            
                        
            for i in range(len(self.nb.param_labels)):

                control0_tmp = wx.StaticText(self.nb.ParameterPage, label = self.nb.param_labels[i], pos=(20, 23*(i+1)+25), size=(120, 20))
                self.nb.ParameterPage.control0.append(control0_tmp)
                
                control1_tmp = wx.TextCtrl(self.nb.ParameterPage,2000+i, pos=(150, 23*(i+1)+20), size=(120,20))
                self.nb.ParameterPage.control1.append(control1_tmp)
                self.nb.ParameterPage.control1[i].SetValue(str(round(self.nb.parameter[self.nb.param_labels[i]][0], 12)))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparvalue, self.nb.ParameterPage.control1[i])
                
                control2_tmp = (wx.TextCtrl(self.nb.ParameterPage,2000+i+len(self.nb.param_labels), pos=(280, 23*(i+1)+20), size=(80,20)))
                self.nb.ParameterPage.control2.append(control2_tmp)
                self.nb.ParameterPage.control2[i].SetValue(str(self.nb.parameter[self.nb.param_labels[i]][1]))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparmin, self.nb.ParameterPage.control2[i])
                
                control3_tmp = (wx.TextCtrl(self.nb.ParameterPage,2000+i+2*len(self.nb.param_labels), pos=(370, 23*(i+1)+20), size=(80,20)))
                self.nb.ParameterPage.control3.append(control3_tmp)
                self.nb.ParameterPage.control3[i].SetValue(str(self.nb.parameter[self.nb.param_labels[i]][2]))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparmax, self.nb.ParameterPage.control3[i])
                                
                control4_tmp = (wx.CheckBox(self.nb.ParameterPage,2000+i+3*len(self.nb.param_labels), label = '', pos = (460, 23*(i+1)+25)))
                self.nb.ParameterPage.control4.append(control4_tmp)
                self.nb.ParameterPage.control4[i].SetValue(self.nb.parameter[self.nb.param_labels[i]][3])
                wx.EVT_CHECKBOX(self.nb.ParameterPage, self.nb.ParameterPage.control4[i].GetId(), self.nb.ParameterPage.editparstate)

                control5_tmp = wx.Button(self.nb.ParameterPage, 2000+i+4*len(self.nb.param_labels), label = '<', pos = (520, 23*(i+1)+20), size = (20,20))
                self.nb.ParameterPage.control5.append(control5_tmp)
                self.Bind(wx.EVT_BUTTON, self.nb.ParameterPage.toggleminus , self.nb.ParameterPage.control5[i])

                control6_tmp = (wx.TextCtrl(self.nb.ParameterPage,2000+i+5*len(self.nb.param_labels), pos=(550, 23*(i+1)+20), size=(50,20)))
                self.nb.ParameterPage.control6.append(control6_tmp)
                self.nb.ParameterPage.control6[i].SetValue('0')
                self.nb.ParameterPage.togglesteps.append(0)
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.togglestep, self.nb.ParameterPage.control6[i])

                control7_tmp = wx.Button(self.nb.ParameterPage,2000+i+6*len(self.nb.param_labels), label = '>', pos = (610, 23*(i+1)+20), size = (20,20))
                self.nb.ParameterPage.control7.append(control7_tmp)
                self.Bind(wx.EVT_BUTTON, self.nb.ParameterPage.toggleplus , self.nb.ParameterPage.control7[i])

            self.nb.ParameterPage.SetScrollbars(0, 10, 0, int((len(self.nb.param_labels)+4)*2.3)+1)
            self.nb.SetSelection(0)
        dlg.Destroy()
        

    def OnReadRigidbody(self,e):
        """ Read in rigid body file"""
        dlg = wx.FileDialog(self, "Choose a rigid body file", self.dirname[4], ".rbf", "*.rbf", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename[4] = dlg.GetFilename()
            self.dirname[4] = dlg.GetDirectory()
            self.nb.rigid_bodies = read_rigid_bodies(self.dirname[4]+'/'+self.filename[4])
            self.nb.MainControlPage.rigidbodyfile.SetValue(self.filename[4])
        dlg.Destroy()

    def OnReadBV(self, e):
        """ Read in a Bond Valence File """
        dlg = wx.FileDialog(self, "Choose a bond valence file", self.dirname[5], ".bvf", "*.bvf", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename[5] = dlg.GetFilename()
            self.dirname[5] = dlg.GetDirectory()
            self.nb.MainControlPage.BVclusters = []
            self.nb.MainControlPage.BVclusters = read_BV((self.dirname[5]+'/'+self.filename[5]), self.nb.cell)
            self.nb.MainControlPage.bondvalencefile.SetValue(self.filename[5])
        dlg.Destroy()

    def LoadSession(self,e):
        self.filename = ['','','','','','']
        self.dirname = ['','','','','','']
        dlg = wx.FileDialog(self, "Load files saved in a session", '', ".ssn", "*.ssn", wx.OPEN)
        ok = True
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            f = file(dirname+'/'+filename, 'r')
            paths = f.readlines()
            f.close()
            for i in range(len(paths)):
                path = str.rsplit(paths[i])
                if len(path) == 2:
                    self.dirname[i] = path[0]
                    self.filename[i] = path[1]
                    if i < 4:
                        try:
                            open(self.dirname[i]+'/'+self.filename[i])
                        except IOError:
                            print 'Cannot open file '+ self.filename[i]
                            ok = False
        dlg.Destroy()
        if ok:
            #read data
            self.nb.data = read_data(self.dirname[0]+'/'+self.filename[0])
            self.nb.MainControlPage.Rdata = calc_Rdata(self.nb.data)
            self.nb.MainControlPage.datafile.SetValue(self.filename[0])
            self.nb.MainControlPage.Rod_weight = []
            self.nb.MainControlPage.rodweight = []
            self.nb.MainControlPage.rmsflaglog.SetValue(self.nb.MainControlPage.RMS_flag_log[self.nb.MainControlPage.RMS_flag]+self.nb.MainControlPage.Rdata[self.nb.MainControlPage.RMS_flag])
            for i in range(len(self.nb.data)):
                wx.StaticText(self.nb.MainControlPage, label = (str(int(self.nb.data[i].H))+' '+str(int(self.nb.data[i].K))+' L'), pos=(350,25*i+67), size=(40,20))
                self.nb.MainControlPage.rodweight.append(wx.TextCtrl(self.nb.MainControlPage,1000+i, pos=(400,25*i+65), size=(30,20)))
                self.nb.MainControlPage.rodweight[i].SetValue('1')
                self.nb.MainControlPage.Rod_weight.append(1)
                self.Bind(wx.EVT_TEXT, self.nb.MainControlPage.setrodweight, self.nb.MainControlPage.rodweight[i])
            #read bulk
            self.nb.bulk, self.nb.cell, self.nb.NLayers = read_bulk(self.dirname[1]+'/'+self.filename[1])
            self.nb.MainControlPage.bulkfile.SetValue(self.filename[1])
            self.nb.MainControlPage.bulkfile.AppendText('\na = '+str(round(self.nb.cell[0],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\nb = '+str(round(self.nb.cell[1],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\nc = '+str(round(self.nb.cell[2],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\nalpha = '+str(round(self.nb.cell[3],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\nbeta = '+str(round(self.nb.cell[4],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\ngamma = '+str(round(self.nb.cell[5],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\ndelta1 = '+str(round(self.nb.cell[6],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\ndelta2 = '+str(round(self.nb.cell[7],5)))
            self.nb.MainControlPage.bulkfile.AppendText('\nNLayers = '+str(self.nb.NLayers))
            self.nb.g_inv = calc_g_inv(self.nb.cell)
            self.nb.ResonantDataPage.allrasd.cell = self.nb.cell
            self.nb.ResonantDataPage.allrasd.g_inv = self.nb.g_inv
            if self.nb.data != []:
                for x in self.nb.data:
                    x.calcFbulk(self.nb.cell, self.nb.bulk, self.nb.g_inv, database)                
            # read surface
            self.nb.surface, self.nb.parameter_usage, self.nb.runningDB = read_surface(self.dirname[2]+'/'+self.filename[2], database)
            self.nb.MainControlPage.surfacefile.SetValue(self.filename[2])
            # read parameters
            self.nb.parameter, self.nb.param_labels = read_parameters(self.dirname[3]+'/'+self.filename[3])
            self.nb.MainControlPage.parameterfile.SetValue(self.filename[3])
            self.nb.DeletePage(2)
            self.nb.ParameterPage = ParameterPanel(self.nb)
            self.nb.AddPage(self.nb.ParameterPage, " Parameters " )            
            for i in range(len(self.nb.param_labels)):

                control0_tmp = wx.StaticText(self.nb.ParameterPage, label = self.nb.param_labels[i], pos=(20, 23*(i+1)+25), size=(120, 20))
                self.nb.ParameterPage.control0.append(control0_tmp)
                
                control1_tmp = wx.TextCtrl(self.nb.ParameterPage,2000+i, pos=(150, 23*(i+1)+20), size=(120,20))
                self.nb.ParameterPage.control1.append(control1_tmp)
                self.nb.ParameterPage.control1[i].SetValue(str(round(self.nb.parameter[self.nb.param_labels[i]][0], 12)))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparvalue, self.nb.ParameterPage.control1[i])
                
                control2_tmp = (wx.TextCtrl(self.nb.ParameterPage,2000+i+len(self.nb.param_labels), pos=(280, 23*(i+1)+20), size=(80,20)))
                self.nb.ParameterPage.control2.append(control2_tmp)
                self.nb.ParameterPage.control2[i].SetValue(str(self.nb.parameter[self.nb.param_labels[i]][1]))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparmin, self.nb.ParameterPage.control2[i])
                
                control3_tmp = (wx.TextCtrl(self.nb.ParameterPage,2000+i+2*len(self.nb.param_labels), pos=(370, 23*(i+1)+20), size=(80,20)))
                self.nb.ParameterPage.control3.append(control3_tmp)
                self.nb.ParameterPage.control3[i].SetValue(str(self.nb.parameter[self.nb.param_labels[i]][2]))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparmax, self.nb.ParameterPage.control3[i])
                                
                control4_tmp = (wx.CheckBox(self.nb.ParameterPage,2000+i+3*len(self.nb.param_labels), label = '', pos = (460, 23*(i+1)+25)))
                self.nb.ParameterPage.control4.append(control4_tmp)
                self.nb.ParameterPage.control4[i].SetValue(self.nb.parameter[self.nb.param_labels[i]][3])
                wx.EVT_CHECKBOX(self.nb.ParameterPage, self.nb.ParameterPage.control4[i].GetId(), self.nb.ParameterPage.editparstate)

                control5_tmp = wx.Button(self.nb.ParameterPage, 2000+i+4*len(self.nb.param_labels), label = '<', pos = (520, 23*(i+1)+20), size = (20,20))
                self.nb.ParameterPage.control5.append(control5_tmp)
                self.Bind(wx.EVT_BUTTON, self.nb.ParameterPage.toggleminus , self.nb.ParameterPage.control5[i])

                control6_tmp = (wx.TextCtrl(self.nb.ParameterPage,2000+i+5*len(self.nb.param_labels), pos=(550, 23*(i+1)+20), size=(50,20)))
                self.nb.ParameterPage.control6.append(control6_tmp)
                self.nb.ParameterPage.control6[i].SetValue('0')
                self.nb.ParameterPage.togglesteps.append(0)
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.togglestep, self.nb.ParameterPage.control6[i])

                control7_tmp = wx.Button(self.nb.ParameterPage,2000+i+6*len(self.nb.param_labels), label = '>', pos = (610, 23*(i+1)+20), size = (20,20))
                self.nb.ParameterPage.control7.append(control7_tmp)
                self.Bind(wx.EVT_BUTTON, self.nb.ParameterPage.toggleplus , self.nb.ParameterPage.control7[i])

            self.nb.ParameterPage.SetScrollbars(0, 10, 0, int((len(self.nb.param_labels)+4)*2.3)+1)
            #read rigid bodies
            if self.filename[4] != '':
                self.nb.rigid_bodies = read_rigid_bodies(self.dirname[4]+'/'+self.filename[4])
                self.nb.MainControlPage.rigidbodyfile.SetValue(self.filename[4])
            # read BV clusters
            if self.filename[5] != '':
                self.nb.MainControlPage.BVclusters = []
                self.nb.MainControlPage.BVclusters = read_BV((self.dirname[5]+'/'+self.filename[5]), self.nb.cell)
                self.nb.MainControlPage.bondvalencefile.SetValue(self.filename[5])
            self.nb.SetSelection(0)
        pass
            
    def OnExit(self,e):
        self.Close(True)  # Close the frame.

    def OnWriteCif(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Write surface structure to .cif file", dirname, ".cif", "*.cif", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            write_cif(self.nb.cell, self.nb.surface, self.nb.parameter,self.nb.parameter_usage, self.nb.rigid_bodies, self.nb.MainControlPage.UBW_flag, self.nb.MainControlPage.use_lay_el, filename)
        dlg.Destroy()

    def OnWritePar(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Write parameters to .par file", dirname, ".par", "*.par", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            write_par(self.nb.parameter, self.nb.param_labels, filename)
            self.filename[3] = filename
            self.dirname[3] = dirname
            self.nb.MainControlPage.parameterfile.SetValue(self.filename[3])
        dlg.Destroy()

    def OnWriteData(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Write Data and Fit Results to .dat file", dirname, ".dat", "*.dat", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            write_data(self.nb.data, filename)
        dlg.Destroy()

    def OnWriteEdens(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Write E-density plot curves to a file", dirname, "", "*.*", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            self.Figure2 = plot_edensity(self.nb.MainControlPage.Figure2, self.nb.surface, self.nb.parameter,\
                                         self.nb.parameter_usage, self.nb.cell, database, self.nb.rigid_bodies,\
                                         self.nb.MainControlPage.UBW_flag, self.nb.MainControlPage.use_lay_el,\
                                         self.nb.ResonantDataPage.resel, filename)
            self.Figure2.canvas.draw()
        dlg.Destroy()
        
    def OnWriteSurf(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Write surface structure to .sur file", dirname, ".sur", "*.sur", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            write_surface(self.nb.cell, self.nb.surface, self.nb.parameter, self.nb.parameter_usage, self.nb.rigid_bodies, self.nb.MainControlPage.UBW_flag, self.nb.MainControlPage.use_lay_el, filename)
        dlg.Destroy()

    def SaveSession(self,e):
        dlg = wx.FileDialog(self,"Save paths to currently loaded files for fast restart",'',".ssn","*.ssn", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            f = file(dirname+'/'+filename, 'w')
            for i in range(len(self.filename)):
                f.write(str(self.dirname[i])+'  '+str(self.filename[i])+'\n')
            f.close()
        dlg.Destroy()
        pass

    def OnReadf1f2(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Read f1f2 data from .f1f2 file", dirname, ".f1f2", "*.f1f2", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename7 = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            self.nb.ResonantDataPage.allrasd = read_f1f2(self.nb.ResonantDataPage.allrasd,dirname+'/'+self.filename7)
            self.nb.ResonantDataPage.f1f2file.SetValue(self.filename7)
            self.nb.SetSelection(1)
        dlg.Destroy()
        
    def OnReadRids(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Specify the file with a list of *.rsd filenames", dirname, "","*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            f = file(dirname+'/'+filename,'r')
            rasddata = f.readlines()
            f.close()
            self.nb.ResonantDataPage.allrasd = read_RSD(self.nb.ResonantDataPage.allrasd, self.nb.bulk, self.nb.surface, self.nb.parameter,\
                                                        self.nb.parameter_usage, self.nb.rigid_bodies, database, rasddata,\
                                                        self.nb.MainControlPage.UBW_flag, self.nb.MainControlPage.use_lay_el,\
                                                        self.nb.MainControlPage.el, dirname, self)

            n = len(self.nb.ResonantDataPage.allrasd.list)
            for i in range(n):
                self.nb.ResonantDataPage.allrasd.list[i] = RASD_Fourier(self.nb.ResonantDataPage.allrasd, i)

                AR = self.nb.ResonantDataPage.allrasd.list[i].AR
                PR = self.nb.ResonantDataPage.allrasd.list[i].PR

                if AR < 0 and PR < 0 and PR > -0.5:
                    PR = PR + 0.5
                    AR = -AR
                elif AR > 0 and PR < 0:
                    PR = PR + 1
                elif AR < 0 and PR > 0.5:
                    AR = -AR
                    PR = PR -0.5
                elif AR < 0 and PR < 0.5 and PR > 0:
                    AR = -AR
                    PR = PR +0.5
                if PR > 1: PR = PR -1
                self.nb.ResonantDataPage.allrasd.list[i].AR = AR
                self.nb.ResonantDataPage.allrasd.list[i].PR = PR

                self.nb.ResonantDataPage.allrasd.RMS = self.nb.ResonantDataPage.allrasd.RMS + Num.sum(self.nb.ResonantDataPage.allrasd.list[i].delta_F)
                 
                
                control0_tmp = wx.Button(self.nb.ResonantDataPage, 3000+i, label = self.nb.ResonantDataPage.allrasd.list[i].file, pos = (290, 23*(i+1)+20), size = (120,20))
                self.nb.ResonantDataPage.datacontrol0.append(control0_tmp)
                self.Bind(wx.EVT_BUTTON, self.nb.ResonantDataPage.ClickDataButton , self.nb.ResonantDataPage.datacontrol0[i])

                control1_tmp = wx.TextCtrl(self.nb.ResonantDataPage,3000+i+n, pos=(420, 23*(i+1)+20), size=(60,20), style = wx.TE_READONLY)
                self.nb.ResonantDataPage.datacontrol1.append(control1_tmp)
                self.nb.ResonantDataPage.datacontrol1[i].SetValue(str(round(AR, 3)))

                control2_tmp = wx.TextCtrl(self.nb.ResonantDataPage,3000+i+2*n, pos=(490, 23*(i+1)+20), size=(60,20))
                self.nb.ResonantDataPage.datacontrol2.append(control2_tmp)
                self.nb.ResonantDataPage.datacontrol2[i].SetValue(str(round(PR, 3)))
                self.Bind(wx.EVT_TEXT, self.nb.ResonantDataPage.editPR, self.nb.ResonantDataPage.datacontrol2[i])

                control4_tmp = wx.CheckBox(self.nb.ResonantDataPage,3000+i+4*n, label = '', pos = (560, 23*(i+1)+25))
                self.nb.ResonantDataPage.datacontrol4.append(control4_tmp)
                self.nb.ResonantDataPage.datacontrol4[i].SetValue(False)
                wx.EVT_CHECKBOX(self.nb.ResonantDataPage, self.nb.ResonantDataPage.datacontrol4[i].GetId(), self.nb.ResonantDataPage.editUIF)

                control5_tmp = wx.CheckBox(self.nb.ResonantDataPage,3000+i+5*n, label = '', pos = (590, 23*(i+1)+25))
                self.nb.ResonantDataPage.datacontrol5.append(control5_tmp)
                self.nb.ResonantDataPage.datacontrol5[i].SetValue(False)
                wx.EVT_CHECKBOX(self.nb.ResonantDataPage, self.nb.ResonantDataPage.datacontrol5[i].GetId(), self.nb.ResonantDataPage.editUIR)

                control6_tmp = wx.TextCtrl(self.nb.ResonantDataPage,3000+i+6*n, pos=(620, 23*(i+1)+20), size=(60,20), style = wx.TE_READONLY)
                self.nb.ResonantDataPage.datacontrol6.append(control6_tmp)
                self.nb.ResonantDataPage.datacontrol6[i].SetValue('')

                control7_tmp = wx.TextCtrl(self.nb.ResonantDataPage,3000+i+7*n, pos=(690, 23*(i+1)+20), size=(60,20), style = wx.TE_READONLY)
                self.nb.ResonantDataPage.datacontrol7.append(control7_tmp)
                self.nb.ResonantDataPage.datacontrol7[i].SetValue('')

            self.nb.ResonantDataPage.SetScrollbars(0, 10, 0, int((n+4)*2.3)+1)
            self.nb.SetSelection(1)
            self.nb.ResonantDataPage.allrasd.RMS = self.nb.ResonantDataPage.allrasd.RMS / self.nb.ResonantDataPage.allrasd.ndata
            message = 'best obtainable chi**2 for a fit on all resonant data is: '+str(self.nb.ResonantDataPage.allrasd.RMS)
            print message
        dlg.Destroy()

    def OnWriteRids(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Write RIDS data and fit results to a file", dirname, "", "*.*", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            write_rids(self.nb.ResonantDataPage.allrasd, filename)
        dlg.Destroy()
        
        
############################################################################################################
class CtrNotebook(wx.Notebook):
    def __init__(self, parent):
        wx.Notebook.__init__(self,parent)

        self.parameter_usage = []
        self.surface = []
        self.runningDB = {}
        self.bulk = []
        self.rigid_bodies = []
        self.data= []
        self.NLayers = 1
        self.cell = []
        self.g_inv = []

        self.parameter = {}
        self.param_labels = []
        self.frame = self.GetParent()

        self.MainControlPage = MainControlPanel(self)
        self.ParameterPage = ParameterPanel(self)
        self.ResonantDataPage = ResonantDataPanel(self)

        self.AddPage(self.MainControlPage, " Main Controls " )
        self.AddPage(self.ResonantDataPage, "Resonant Data ")
        self.AddPage(self.ParameterPage, " Parameters ")
        

        
############################################################################################################
class MainControlPanel(wx.Panel):
    def __init__(self,parent):
        wx.Panel.__init__(self,parent)
        self.doplotbulk = False
        self.doplotsurf = False
        self.doplotrough = False
        self.doplotwater = False
        self.plotdims = [2,3]

        self.BVclusters = []
        self.use_BVC = False

        self.Rod_weight = [] #List of Rod weights
        self.param_best = {}
        self.RMS = -1
        self.Rdata = ['','','','','']
        self.RMS_flag = 1
        self.RMS_options = ['1','2','3','4']
        self.RMS_flag_log = [' invalid R flag ',' R = norm. average of |log(F) - log(Fcalc)|, Rdata = ',' R = norm. average |F - Fcalc|, Rdata = ',' R = norm. Ferr weighted average of |F - Fcalc|, Rdata = ', ' R = norm. chi**2, Rdata = ']
        self.UBW_flag = False
        self.use_lay_el = False
        self.el = 'h'
        
        self.sim_an_params = [50, 20, 0.7, 10, 0.01, 500000, False]
        self.simplex_params = [1.0,0.5,2.0,1.0,0.0001,0.001,10000, False]
        self.genetic_params = [0.7,0.5,3.0,1000,0.5,10,0.001,0.0001, False]

        self.opt_options = ['Simulated Annealing', 'Genetic Algorithm', 'Downhill Simplex']
        self.opt_flag = 2
        self.nb = self.GetParent()

        self.xf = 1.
        self.yf = 1.
        self.zf = 3.
        self.an = 10
        self.bn = 10
        self.cn = 60
        self.zmin = 0.
        self.flag = 0
        self.DFMoptions = ['F_obs', 'F model', 'F_bulk', 'F obs - F model', 'F obs - F bulk']
        self.fframe = None
        
        self.Figure1 = None
        self.Figure2 = None
        self.Figure3 = None
        self.StopFit = False

        #Panel Headings
        wx.StaticText(self, label = 'File Information: ', pos=(20, 12), size=(100, 20))
        wx.StaticText(self, label = 'Data Information: ', pos=(290, 12), size=(100, 20))
        wx.StaticText(self, label = 'Fit options: ', pos=(520, 12), size=(100, 20))
        wx.StaticLine(self, pos = (270,0), size = (5,650), style = wx.LI_VERTICAL)
        wx.StaticLine(self, pos = (500,0), size = (5,650), style = wx.LI_VERTICAL)
        wx.StaticLine(self, pos = (0,265), size = (270,5), style = wx.LI_HORIZONTAL)
        wx.StaticLine(self, pos = (0,515), size = (270,5), style = wx.LI_HORIZONTAL)
        wx.StaticLine(self, pos = (270,425), size = (230,5), style = wx.LI_HORIZONTAL)
        wx.StaticLine(self, pos = (500,165), size = (285,5), style = wx.LI_HORIZONTAL)
        wx.StaticLine(self, pos = (500,280), size = (285,5), style = wx.LI_HORIZONTAL)
        wx.StaticLine(self, pos = (500,420), size = (285,5), style = wx.LI_HORIZONTAL)
        wx.StaticLine(self, pos = (500,560), size = (285,5), style = wx.LI_HORIZONTAL)

        wx.StaticText(self, label = 'Rods:', pos=(350,42), size=(40,20))
        wx.StaticText(self, label = 'weight', pos=(400,42), size=(30,20))
        #rod weights
        self.rodweight = []# List of TextControls to read rod weights

        #display uploaded files
        wx.StaticText(self, label = 'data file = ', pos=(20, 42), size=(100, 20))
        self.datafile = wx.TextCtrl(self, pos=(130,40), size=(120,20), style= wx.TE_READONLY)

        wx.StaticText(self, label = 'bulk file = ', pos=(20, 67), size=(100, 20))
        self.bulkfile = wx.TextCtrl(self, pos=(130,65), size=(120,90), style=wx.TE_MULTILINE | wx.TE_READONLY)

        wx.StaticText(self, label = 'surface file = ', pos=(20,162), size=(100, 20))
        self.surfacefile = wx.TextCtrl(self, pos=(130,160), size=(120,20), style= wx.TE_READONLY)

        wx.StaticText(self, label = 'parameter file = ', pos=(20, 187), size=(100, 20))
        self.parameterfile = wx.TextCtrl(self, pos=(130,185), size=(120,20), style= wx.TE_READONLY)

        wx.StaticText(self, label = 'rigid body file = ', pos=(20, 212), size=(100, 20))
        self.rigidbodyfile = wx.TextCtrl(self, pos=(130,210), size=(120,20), style= wx.TE_READONLY)

        wx.StaticText(self, label = 'bond valence file = ', pos=(20, 237), size=(100, 20))
        self.bondvalencefile = wx.TextCtrl(self, pos=(130,235), size=(120,20), style= wx.TE_READONLY)

        #CTR plotting
        wx.StaticText(self, label = 'Plotting Options: ', pos=(20, 282), size=(100, 20))
        wx.StaticText(self, label = 'CTR plot options (Fig. 1): ', pos=(20, 307), size=(200, 20))

        wx.StaticText(self, label = 'CTR plot dimensions', pos=(40, 332), size=(110, 20))
        self.getplotdims = wx.TextCtrl(self, pos=(160, 330), size=(40,-1))
        self.getplotdims.SetValue(str(self.plotdims[0])+' '+str(self.plotdims[1]))
        self.Bind(wx.EVT_TEXT, self.setplotdims, self.getplotdims)

        self.plot_bulk = wx.CheckBox(self, label =  '  plot bulk (green)', pos = (40, 357))
        self.plot_bulk.SetValue(False)
        wx.EVT_CHECKBOX(self, self.plot_bulk.GetId(), self.setplotbulk)
        
        self.plot_surf = wx.CheckBox(self, label =  '  plot surface (red)', pos = (40, 382))
        self.plot_surf.SetValue(False)
        wx.EVT_CHECKBOX(self, self.plot_surf.GetId(), self.setplotsurf)
        
        self.plot_rough = wx.CheckBox(self, label = '  plot roughness (cyan)', pos = (40, 407))
        self.plot_rough.SetValue(False)
        wx.EVT_CHECKBOX(self, self.plot_rough.GetId(), self.setplotrough)
        
        self.plot_water = wx.CheckBox(self, label = '  plot water (magenta)', pos = (40, 432))
        self.plot_water.SetValue(False)
        wx.EVT_CHECKBOX(self, self.plot_water.GetId(), self.setplotwater)

        self.button = wx.Button(self, label = 'Plot CTRs', pos =(170,425))
        self.Bind(wx.EVT_BUTTON, self.OnClick, self.button)

        #edensity plotting
        wx.StaticText(self, label = 'plot e-density', pos=(20, 455), size=(140, 20))
        wx.StaticText(self, label = '(Fig. 2, current surface)', pos=(20, 474), size=(140, 20))
        self.plotedens = wx.Button(self, label = 'Plot edens', pos =(170,465))
        self.Bind(wx.EVT_BUTTON, self.OnPlotedens, self.plotedens)

        #Bond Valence constraints
        wx.StaticText(self, label = 'Bond Valence: ', pos=(20, 542), size=(120, 20))
        self.set_use_BVC = wx.CheckBox(self, label = ' use Bond Valence Constraints in Fit', pos = (25, 570))
        self.set_use_BVC.SetValue(False)
        wx.EVT_CHECKBOX(self, self.set_use_BVC.GetId(), self.setusebvc)
        self.calc_BVS = wx.Button(self, label = 'calculate Bond Valence Sums', pos =(20,600), size =(230, 35))
        self.Bind(wx.EVT_BUTTON, self.on_calc_BVS, self.calc_BVS)

        #Fourier Mapping
        wx.StaticText(self, label = 'Fourier Maps:', pos=(290, 440), size=(160, 20))
        wx.StaticText(self, label = 'cells in x', pos=(290, 467), size=(40, 20))
        self.getxf = wx.TextCtrl(self, pos=(340,465), size=(40,20))
        self.getxf.SetValue(str(self.xf))
        self.Bind(wx.EVT_TEXT, self.setxf, self.getxf)

        wx.StaticText(self, label = '         y', pos=(290, 492), size=(40, 20))
        self.getyf = wx.TextCtrl(self, pos=(340,490), size=(40,20))
        self.getyf.SetValue(str(self.yf))
        self.Bind(wx.EVT_TEXT, self.setyf, self.getyf)

        wx.StaticText(self, label = '         z', pos=(290, 517), size=(40, 20))
        self.getzf = wx.TextCtrl(self, pos=(340,515), size=(40,20))
        self.getzf.SetValue(str(self.zf))
        self.Bind(wx.EVT_TEXT, self.setzf, self.getzf)

        wx.StaticText(self, label = '      zmin', pos=(290, 542), size=(40, 20))
        self.getzmin = wx.TextCtrl(self, pos=(340,540), size=(40,20))
        self.getzmin.SetValue(str(self.zmin))
        self.Bind(wx.EVT_TEXT, self.setzmin, self.getzmin)

        wx.StaticText(self, label = 'pixels in x', pos=(390, 467), size=(50, 20))
        self.getan = wx.TextCtrl(self, pos=(450,465), size=(40,20))
        self.getan.SetValue(str(self.an))
        self.Bind(wx.EVT_TEXT, self.setan, self.getan)

        wx.StaticText(self, label = '          y', pos=(390, 492), size=(50, 20))
        self.getbn = wx.TextCtrl(self, pos=(450,490), size=(40,20))
        self.getbn.SetValue(str(self.bn))
        self.Bind(wx.EVT_TEXT, self.setbn, self.getbn)

        wx.StaticText(self, label = '          z', pos=(390, 517), size=(50, 20))
        self.getcn = wx.TextCtrl(self, pos=(450,515), size=(40,20))
        self.getcn.SetValue(str(self.cn))
        self.Bind(wx.EVT_TEXT, self.setcn, self.getcn)

        self.DFMchoice = wx.ComboBox(self,-1, value=self.DFMoptions[self.flag],\
                                     pos=(340,570),size=(150,20),\
                                     choices=self.DFMoptions,style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.SelectDFM,self.DFMchoice)

        self.DFMButton1 = wx.Button(self, label = 'Fourier Synthesis', pos =(290,600), size=(200,35))
        self.Bind(wx.EVT_BUTTON, self.OnClickDFM1, self.DFMButton1)

        # Fitting options, choose optimization method
        self.opt_choice = wx.ComboBox(self,-1, value=self.opt_options[self.opt_flag],\
                                     pos=(520,35),size=(240,20),\
                                     choices=self.opt_options,style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.Select_opt,self.opt_choice)
        
        #Fitting option use bulk water
        self.douse_bulk_water = wx.CheckBox(self, label = '  use bulk water', pos = (520, 65))
        self.douse_bulk_water.SetValue(False)
        wx.EVT_CHECKBOX(self, self.douse_bulk_water.GetId(), self.set_bulk_water)

        #Fitting option use random start parameters 
        self.random_pars = wx.CheckBox(self, label = '  start Fit with random parameters', pos = (520, 90))
        self.random_pars.SetValue(False)
        wx.EVT_CHECKBOX(self, self.random_pars.GetId(), self.setrandom_pars)
        
        # R-flags
        wx.StaticText(self, label = 'R flag', pos=(520, 117), size=(140, 20))
        self.getrmsflag = wx.ComboBox(self,-1, value=self.RMS_options[self.RMS_flag-1],\
                                     pos=(700, 115), size=(60,20),\
                                     choices=self.RMS_options,style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.setrmsflag, self.getrmsflag)
        self.rmsflaglog = wx.TextCtrl(self, pos=(520, 140), size=(240,20), style = wx.TE_READONLY)
        self.rmsflaglog.SetValue(self.RMS_flag_log[self.RMS_flag]+self.Rdata[self.RMS_flag])

        #Simulated annealing parameters and options
        wx.StaticText(self, label = 'Simulated Annealing Parameters:  ', pos=(520, 182), size=(240, 20))

        wx.StaticText(self, label = 'T start ', pos=(520, 207), size=(40, 20))
        self.Tstart = wx.TextCtrl(self, pos=(570,205), size=(60,20))
        self.Tstart.SetValue(str(self.sim_an_params[0]))
        self.Bind(wx.EVT_TEXT, self.setTstart, self.Tstart)

        wx.StaticText(self, label = 'T end ', pos=(640, 207), size=(50, 20))
        self.Tend = wx.TextCtrl(self, pos=(700,205), size=(60,20))
        self.Tend.SetValue(str(self.sim_an_params[1]))
        self.Bind(wx.EVT_TEXT, self.setTend, self.Tend)

        wx.StaticText(self, label = 'cooling ', pos=(520, 232), size=(40, 20))
        self.cool = wx.TextCtrl(self, pos=(570,230), size=(60,20))
        self.cool.SetValue(str(self.sim_an_params[2]))
        self.Bind(wx.EVT_TEXT, self.setcool, self.cool)

        wx.StaticText(self, label = 'iter/T ', pos=(640, 232), size=(50, 20))
        self.iter = wx.TextCtrl(self, pos=(700,230), size=(60,20))
        self.iter.SetValue(str(self.sim_an_params[3]))
        self.Bind(wx.EVT_TEXT, self.setiter, self.iter)

        wx.StaticText(self, label = 'Boltzmann ', pos=(640, 257), size=(50, 20))
        self.Boltz = wx.TextCtrl(self, pos=(700,255), size=(60,20))
        self.Boltz.SetValue(str(self.sim_an_params[5]))
        self.Bind(wx.EVT_TEXT, self.setBoltz, self.Boltz)
        
        #### Genetic algorithm parameters and options #################################### 
        wx.StaticText(self, label = 'Genetic Algorithm Parameters:  ', pos=(520, 295), size=(200, 20))

        wx.StaticText(self, label = 'mutate', pos=(520, 322), size=(40, 20))
        self.km = wx.TextCtrl(self, pos=(570,320), size=(60,20))
        self.km.SetValue(str(self.genetic_params[0]))
        self.Bind(wx.EVT_TEXT, self.setkm, self.km)

        wx.StaticText(self, label = 'crossover', pos=(640, 322), size=(50, 20))
        self.kc = wx.TextCtrl(self, pos=(700,320), size=(60,20))
        self.kc.SetValue(str(self.genetic_params[1]))
        self.Bind(wx.EVT_TEXT, self.setkc, self.kc)

        wx.StaticText(self, label = 'jump', pos=(520, 347), size=(40, 20))
        self.jump = wx.TextCtrl(self, pos=(570,345), size=(60,20))
        self.jump.SetValue(str(self.genetic_params[4]))
        self.Bind(wx.EVT_TEXT, self.setjump, self.jump)

        wx.StaticText(self, label = 'popmult', pos=(640, 347), size=(50, 20))
        self.popmult = wx.TextCtrl(self, pos=(700,345), size=(60,20))
        self.popmult.SetValue(str(self.genetic_params[2]))
        self.Bind(wx.EVT_TEXT, self.setpopmult, self.popmult)

        wx.StaticText(self, label = 'Ftol', pos=(520, 372), size=(40, 20))
        self.genftol = wx.TextCtrl(self, pos=(570,370), size=(60,20))
        self.genftol.SetValue(str(self.genetic_params[7]))
        self.Bind(wx.EVT_TEXT, self.setgenftol, self.genftol)

        wx.StaticText(self, label = 'Xtol', pos=(640, 372), size=(50, 20))
        self.genxtol = wx.TextCtrl(self, pos=(700,370), size=(60,20))
        self.genxtol.SetValue(str(self.genetic_params[6]))
        self.Bind(wx.EVT_TEXT, self.setgenxtol, self.genxtol)

        wx.StaticText(self, label = 'stagn.', pos=(520, 397), size=(40, 20))
        self.stagnate = wx.TextCtrl(self, pos=(570,395), size=(60,20))
        self.stagnate.SetValue(str(self.genetic_params[5]))
        self.Bind(wx.EVT_TEXT, self.setstagnate, self.stagnate)

        wx.StaticText(self, label = 'maxgen', pos=(640, 397), size=(50, 20))
        self.maxgen = wx.TextCtrl(self, pos=(700,395), size=(60,20))
        self.maxgen.SetValue(str(self.genetic_params[3]))
        self.Bind(wx.EVT_TEXT, self.setmaxgen, self.maxgen)
        
        #### Downhill Simplex parameters and options #################################### 
        wx.StaticText(self, label = 'Downhill Simplex Parameters:  ', pos=(520, 435), size=(200, 20))

        wx.StaticText(self, label = 'alpha:    ', pos=(520, 462), size=(40, 20))
        self.alpha = wx.TextCtrl(self, pos=(570,460), size=(60,20))
        self.alpha.SetValue(str(self.simplex_params[0]))
        self.Bind(wx.EVT_TEXT, self.setalpha, self.alpha)

        wx.StaticText(self, label = 'beta:    ', pos=(640, 462), size=(50, 20))
        self.beta = wx.TextCtrl(self, pos=(700,460), size=(60,20))
        self.beta.SetValue(str(self.simplex_params[1]))
        self.Bind(wx.EVT_TEXT, self.setbeta, self.beta)

        wx.StaticText(self, label = 'gamma:    ', pos=(520, 487), size=(40, 20))
        self.gamma = wx.TextCtrl(self, pos=(570,485), size=(60,20))
        self.gamma.SetValue(str(self.simplex_params[2]))
        self.Bind(wx.EVT_TEXT, self.setgamma, self.gamma)

        wx.StaticText(self, label = 'delta:    ', pos=(640, 487), size=(50, 20))
        self.delta = wx.TextCtrl(self, pos=(700,485), size=(60,20))
        self.delta.SetValue(str(self.simplex_params[3]))
        self.Bind(wx.EVT_TEXT, self.setdelta, self.delta)

        wx.StaticText(self, label = 'Ftol:    ', pos=(520, 512), size=(40, 20))
        self.ftol = wx.TextCtrl(self, pos=(570,510), size=(60,20))
        self.ftol.SetValue(str(self.simplex_params[4]))
        self.Bind(wx.EVT_TEXT, self.setftol, self.ftol)

        wx.StaticText(self, label = 'Xtol:    ', pos=(640, 512), size=(50, 20))
        self.xtol = wx.TextCtrl(self, pos=(700,510), size=(60,20))
        self.xtol.SetValue(str(self.simplex_params[5]))
        self.Bind(wx.EVT_TEXT, self.setxtol, self.xtol)

        wx.StaticText(self, label = 'maxiter:    ', pos=(640, 537), size=(50, 20))
        self.maxiter = wx.TextCtrl(self, pos=(700,535), size=(60,20))
        self.maxiter.SetValue(str(self.simplex_params[6]))
        self.Bind(wx.EVT_TEXT, self.setmaxiter, self.maxiter)

        
        # Start and Stop Fit ###################################################################
        self.Startfitbutton = wx.Button(self, label = 'Start Fit', pos =(520,575), size=(170,60))
        self.Bind(wx.EVT_BUTTON, self.OnClickStartFit, self.Startfitbutton)
        self.Stopfitbutton = wx.Button(self, label = 'Stop Fit', pos =(700,575), size=(60,60))
        self.Bind(wx.EVT_BUTTON, self.OnClickStopFit, self.Stopfitbutton)
                                       
   ################################# Plotting options event functions #########################################
    def setplotdims(self, event):
        if event.GetString() == '':
            None
        else:
            dims = str.rsplit(str(event.GetString()))
            if len(dims) != 2:
                print 'need two integer numbers'
            else:
                self.plotdims = []
                for i in dims:
                    self.plotdims.append(int(i))
                
    def setplotbulk(self,e): self.doplotbulk = self.plot_bulk.GetValue()
    def setplotsurf(self,e): self.doplotsurf = self.plot_surf.GetValue()
    def setplotrough(self,e): self.doplotrough = self.plot_rough.GetValue()
    def setplotwater(self,e): self.doplotwater = self.plot_water.GetValue()

    def setrodweight(self,event):
        rod = event.GetId()-1000
        if (event.GetString() == '') or (event.GetString() == '-'):
            None
        else:
            a = float(event.GetString())
            if a < 0:
                print 'rod weight must be >= 0'
            else:
                self.Rod_weight[rod] = a

    def OnClick(self,e):
        t0 = time.clock()
        self.nb.data, self.RMS = calc_CTRs(self.nb.parameter,self.nb.parameter_usage, self.nb.data, self.nb.cell,\
                               self.nb.surface, self.nb.NLayers, database, self.nb.g_inv, self.Rod_weight, self.nb.rigid_bodies, self.UBW_flag, self.use_BVC,\
                               self.BVclusters, self.RMS_flag, self.use_lay_el, self.el)
        print 'time in [s] needed to process calc_CTRs: '+str(time.clock() -t0)
        self.Figure1 = plot_rods(self.Figure1, self.nb.data, self.plotdims, self.doplotbulk, self.doplotsurf, self.doplotrough,\
                                                 self.doplotwater, self.RMS)
        self.Figure1.canvas.draw()
        check_vibes(self.nb.surface,self.nb.parameter, self.nb.parameter_usage)
    def OnPlotedens(self, e):
        self.Figure2 = plot_edensity(self.Figure2, self.nb.surface, self.nb.parameter, self.nb.parameter_usage, self.nb.cell, database, self.nb.rigid_bodies, self.UBW_flag, self.use_lay_el, self.nb.ResonantDataPage.resel)
        self.Figure2.canvas.draw()
        check_vibes(self.nb.surface,self.nb.parameter, self.nb.parameter_usage)
    ################################# Bond Valence constraints event functions #########################################
    def setusebvc(self,e): self.use_BVC = self.set_use_BVC.GetValue()
    def on_calc_BVS(self, e):
        for i in range(len(self.BVclusters)):
            global_parms, surface = param_unfold(self.nb.parameter,self.nb.parameter_usage, self.nb.surface, self.UBW_flag, self.use_lay_el)
            surface = RB_update(self.nb.rigid_bodies, surface, self.nb.parameter, self.nb.cell)
            BVS, dist = self.BVclusters[i].calc_BVS(surface)
            print 'Bond Valence sum in BV-cluster '+str(i+1)+' = '+str(round(BVS, 4))+'; distances are: \n'
            for j in range(len(dist)):
                print str(round(dist[j],5))
            print '\n'
    ###########################Differenz Fourier Map param setting and execution##########################################
    def setxf(self,e):
        a = e.GetString()
        try:
            value = float(a)
            if value > 0 and value <=1:
                self.xf = value
            else:
                print 'the fractional cell size in x direction should be 1/dH_min (>0 and <1)'
        except ValueError:
            pass
    def setyf(self,e):
        a = e.GetString()
        try:
            value = float(a)
            if value > 0 and value <=1:
                self.yf = value
            else:
                print 'the fractional cell size in y direction should be 1/dK_min (>0 and <1)'
        except ValueError:
            pass
    def setzf(self,e):
        a = e.GetString()
        try:
            value = float(a)
            if value > 0:
                self.zf = value
            else:
                print 'the fractional cell size in z direction should be > 0 (~1/dL)'
        except ValueError:
            pass
    def setzmin(self,e):
        a = e.GetString()
        try:
            self.zmin = float(a)
        except ValueError:
            pass
    def setan(self,e):
        a = e.GetString()
        try:
            value = int(a)
            if value > 0:
                self.an = value
            else:
                print 'an must be >=1 '
        except ValueError:
            pass
    def setbn(self,e):
        a = e.GetString()
        try:
            value = int(a)
            if value > 0:
                self.bn = value
            else:
                print 'bn must be >=1 '
        except ValueError:
            pass
    def setcn(self,e):
        a = e.GetString()
        try:
            value = int(a)
            if value > 0:
                self.cn = value
            else:
                print 'cn must be >=1 '
        except ValueError:
            pass
    def SelectDFM(self,e):
        self.flag = e.GetSelection()
    def OnClickDFM1(self,e):
        self.nb.frame.SetStatusText('Preparing Fourier Components')
        if self.flag==0:
            F_u = None
            F_k = create_F_obs(self.nb.data,self.nb.surface,self.nb.g_inv,\
                                 self.nb.cell,self.nb.runningDB,\
                                 self.nb.parameter,self.nb.parameter_usage,\
                                 self.nb.rigid_bodies,self.UBW_flag,\
                                 self.nb.MainControlPage.use_lay_el,\
                                 self.nb.NLayers)
        if self.flag==1:
            F_u = None
            F_k = create_F_model(self.nb.data,self.nb.surface,self.nb.g_inv,\
                                 self.nb.cell,self.nb.runningDB,\
                                 self.nb.parameter,self.nb.parameter_usage,\
                                 self.nb.rigid_bodies,self.UBW_flag,\
                                 self.nb.MainControlPage.use_lay_el)
        elif self.flag==2:
            F_u = None
            F_k= create_F_bulk(self.nb.data)
            
        elif self.flag==3:
            F_u = create_F_data(self.nb.data)
            F_k = create_F_model(self.nb.data,self.nb.surface,self.nb.g_inv,\
                                 self.nb.cell,self.nb.runningDB,\
                                 self.nb.parameter,self.nb.parameter_usage,\
                                 self.nb.rigid_bodies,self.UBW_flag,\
                                 self.nb.MainControlPage.use_lay_el)
        elif self.flag==4:
            F_u = create_F_data(self.nb.data)
            F_k= create_F_bulk(self.nb.data)
        
        self.nb.frame.SetStatusText('Performing Fourier Synthesis')
        rho, Fcell = fdiff2rho(F_k,F_u,self.nb.cell,self.xf,self.yf,\
                               self.zf,self.an,self.bn,self.cn,self.zmin,self.flag,\
                               self)        
        if self.fframe != None:
            self.fframe.Close()
            self.fframe = None
        self.fframe = createfourierframe(self,rho,Fcell,\
                                         self.zmin*self.nb.cell[2])
        self.fframe.Show(True)
        self.nb.frame.SetStatusText('Done, have a nice day!')
        
    ################################# Fitting options event functions #########################################
    def Select_opt(self, e):
        self.opt_flag = e.GetSelection()
        
    def setrandom_pars(self,e):
        self.sim_an_params[6] = self.random_pars.GetValue()
        self.simplex_params[7] = self.random_pars.GetValue()
        self.genetic_params[8] = self.random_pars.GetValue()      
        
    def set_bulk_water(self,event):
        self.UBW_flag = self.douse_bulk_water.GetValue()

    def setrmsflag(self,event):
        self.RMS_flag = event.GetSelection()+1
        self.rmsflaglog.SetValue(self.RMS_flag_log[self.RMS_flag]+self.Rdata[self.RMS_flag])
                                     
    ################################# Simulated Annealing options event functions #########################################
    def setTstart(self,event):
        a = event.GetString()
        try:
            value = float(a)
            if value > self.sim_an_params[1]:
                self.sim_an_params[0] = value
            else:
                print 'T start should be > T end'
        except ValueError:
            pass    
    def setTend(self,event):
        a = event.GetString()
        try:
            value = float(a)
            if value < self.sim_an_params[0]:
                self.sim_an_params[1] = value
            else:
                print 'T end should be < T start'
        except ValueError:
            pass        
    def setcool(self,event):
        a = event.GetString()
        try:
            value = float(a)
            if value > 0 and value < 1:
                self.sim_an_params[2] = value
            else:
                print 'the cooling factor should be between 0 and 1'
        except ValueError:
            pass   
    def setiter(self,event):
        a = event.GetString()
        try:
            value = int(a)
            if value > 0:
                self.sim_an_params[3] = value
            else:
                print 'you must allow at least one iteration at each Temperature'
        except ValueError:
            pass   
    def setBoltz(self,event):
        a = event.GetString()
        try:
            value = float(a)
            if value > 1000:
                self.sim_an_params[5] = value
            else:
                print 'the Boltzmann weight should be > 1000'
        except ValueError:
            pass  
    ###########################################################################################
    def setkc(self, e):
        a = e.GetString()
        try:
            value = float(a)
            if value > 0 and value <=1:
                self.genetic_params[1] = value
            else:
                print 'the crossover probability should be between 0 and 1'
        except ValueError:
            pass
    def setkm(self, e):
        a = e.GetString()
        try:
            value = float(a)
            if value > 0 and value <=1:
                self.genetic_params[0] = value
            else:
                print 'the muataion probability should be between 0 and 1'
        except ValueError:
            pass
    def setjump(self, e):
        a = e.GetString()
        try:
            value = float(a)
            if value > 0 and value <=1:
                self.genetic_params[4] = value
            else:
                print 'the maximum fractional mutation JUMP should be between 0 and 1'
        except ValueError:
            pass
    def setpopmult(self, e):
        a = e.GetString()
        try:
            value = float(a)
            if value >= 1:
                self.genetic_params[2] = value
            else:
                print 'popmult should be > 1'
        except ValueError:
            pass
    def setgenftol(self, e):
        a = e.GetString()
        try:
            value = float(a)
            if value > 0:
                self.genetic_params[7] = value
            else:
                print 'Ftol must be > 0'
        except ValueError:
            pass
    def setgenxtol(self, e):
        a = e.GetString()
        try:
            value = float(a)
            if value > 0:
                self.genetic_params[6] = value
            else:
                print 'Xtol must be > 0'
        except ValueError:
            pass
    def setstagnate(self, e):
        a = e.GetString()
        try:
            value = int(a)
            if value >= 2:
                self.genetic_params[5] = value
            else:
                print 'stagnate must be >= 2'
        except ValueError:
            pass
    def setmaxgen(self, e):
        a = e.GetString()
        try:
            value = int(a)
            if value >= 1:
                self.genetic_params[3] = value
            else:
                print 'maxgen must be >= 1'
        except ValueError:
            pass
    ################################# Simplex options event functions #########################################
    def setalpha(self,event):
        try:
            a = float(event.GetString())
            if a <= 0 or a > 1:
                print 'alpha must be > 0 and <= 1'
            else:
                self.simplex_params[0] = a
        except ValueError:
            pass
    def setbeta(self,event):
        try:
            a = float(event.GetString())
            if a <= 0 or a > 1:
                print 'beta must be > 0 and <= 1'
            else:
                self.simplex_params[1] = a
        except ValueError:
            pass
    def setgamma(self,event):
        try:
            a = float(event.GetString())
            if a <= 0:
                print 'gamma must be > 0'
            else:
                self.simplex_params[2] = a
        except ValueError:
            pass
    def setdelta(self,event):
        try:
            a = float(event.GetString())
            if a <= 0 or a > 1:
                print 'delta must be > 0 and <= 1'
            else:
                self.simplex_params[3] = a
        except ValueError:
            pass
    def setftol(self,event):
        try:
            a = float(event.GetString())
            if a <= 0:
                print 'ftol must be > 0'
            else:
                self.simplex_params[4] = a
        except ValueError:
            pass
    def setxtol(self,event):
        try:
            a = float(event.GetString())
            if a <= 0:
                print 'xtol must be > 0'
            else:
                self.simplex_params[5] = a
        except ValueError:
            pass
    def setmaxiter(self,event):
        try:
            a = int(event.GetString())
            if a <= 0:
                print 'maxiter must be > 0'
            else:
                self.simplex_params[6] = a
        except ValueError:
            pass
    ################################################################################################################################
    def OnClickStartFit(self,e):
        flag = check_model_consistency(self.nb.param_labels, self.nb.parameter, self.nb.parameter_usage, self.nb.rigid_bodies, self.UBW_flag)
        if flag:
            self.param_best = {}
            self.RMS = -1
            self.StopFit = False
            if self.opt_flag == 0:
                self.nb.data, self.param_best, self.RMS = simulated_annealing(self.nb.data, self.nb.cell, self.nb.NLayers, self.nb.surface,\
                                                                          self.nb.runningDB, self.nb.parameter, self.nb.parameter_usage, self.nb.rigid_bodies, self)
                while wx.GetApp().Pending():
                    wx.GetApp().Dispatch()
                    wx.GetApp().Yield(True)
            elif self.opt_flag == 1:
                self.nb.data, self.param_best, self.RMS = genetic(self, self.nb.data, self.nb.cell, self.nb.surface, self.nb.NLayers, self.nb.runningDB,\
                                                                  self.nb.rigid_bodies, self.nb.parameter, self.nb.parameter_usage)
                while wx.GetApp().Pending():
                    wx.GetApp().Dispatch()
                    wx.GetApp().Yield(True)
            elif self.opt_flag == 2:
                self.nb.data, self.param_best, self.RMS = simplex(self.nb.frame,self.nb.parameter, self.nb.parameter_usage, self.nb.data, self.nb.cell, self.nb.surface, \
                                        self.nb.NLayers, self.nb.runningDB, self.nb.rigid_bodies, self)
                while wx.GetApp().Pending():
                    wx.GetApp().Dispatch()
                    wx.GetApp().Yield(True)
            self.Figure1 = plot_rods(self.Figure1,self.nb.data, self.plotdims, self.doplotbulk, self.doplotsurf, self.doplotrough, self.doplotwater, self.RMS)
            self.Figure1.canvas.draw()
            check_vibes(self.nb.surface,self.nb.parameter, self.nb.parameter_usage)
            dlg = wx.MessageDialog(self, "Keep refined parameters ?","", wx.YES_NO | wx.STAY_ON_TOP)
            if dlg.ShowModal() == wx.ID_YES:
                for i in range(len(self.nb.param_labels)):
                    self.nb.parameter[self.nb.param_labels[i]][0] = self.param_best[self.nb.param_labels[i]][0]
                    self.nb.ParameterPage.control1[i].SetValue(str(round(self.param_best[self.nb.param_labels[i]][0], 12)))
            dlg.Destroy()
            self.nb.SetSelection(2)
    def OnClickStopFit(self,e):
        self.StopFit = True

##########################################################################################################################
class ParameterPanel(wx.ScrolledWindow):
    def __init__(self,parent):
        wx.ScrolledWindow.__init__(self,parent)

        self.control0 = []
        self.control1 = []
        self.control2 = []
        self.control3 = []
        self.control4 = []
        self.control5 = []
        self.control6 = []
        self.control7 = []
        self.togglesteps = []

        wx.StaticText(self, label = 'value', pos=(170, 10), size=(100, 20))
        wx.StaticText(self, label = 'min', pos=(300, 10), size=(70, 20))
        wx.StaticText(self, label = 'max', pos=(390, 10), size=(70, 20))
        wx.StaticText(self, label = 'refine', pos=(455, 10), size=(60, 15))
        wx.StaticText(self, label = ' - ', pos=(525, 10), size=(20, 20))
        wx.StaticText(self, label = ' step ', pos=(560, 10), size=(40, 20))
        wx.StaticText(self, label = ' + ', pos=(615, 10), size=(20, 20))

        self.uncheckbutton = wx.Button(self, label = 'u', pos =(459,25), size=(15,15))
        self.Bind(wx.EVT_BUTTON, self.uncheck, self.uncheckbutton)
        
        self.nb = self.GetParent()

    def uncheck(self,e):
        for i in range(len(self.control4)):
            self.nb.parameter[self.nb.param_labels[i]][3] = False
            self.control4[i].SetValue(False)

    def editparvalue(self,event):
        item = event.GetId()-2000
        try:
            self.nb.parameter[self.nb.param_labels[item]][0] = float(event.GetString())
        except ValueError:
            pass        
    def editparmin(self,event):
        item = event.GetId()-len(self.nb.param_labels)-2000
        try:
            self.nb.parameter[self.nb.param_labels[item]][1] = float(event.GetString())
        except ValueError:
            pass        
    def editparmax(self,event):
        item = event.GetId()-2*len(self.nb.param_labels)-2000
        try:
            self.nb.parameter[self.nb.param_labels[item]][2] = float(event.GetString())
        except ValueError:
            pass        
    def editparstate(self,event):
        item = event.GetId()- 3*len(self.nb.param_labels)-2000
        self.nb.parameter[self.nb.param_labels[item]][3] = self.control4[item].GetValue()
    def toggleminus(self, event):
        item = event.GetId()-4*len(self.nb.param_labels)-2000
        step = self.togglesteps[item]
        self.nb.parameter[self.nb.param_labels[item]][0] = self.nb.parameter[self.nb.param_labels[item]][0] - step
        self.control1[item].SetValue(str(self.nb.parameter[self.nb.param_labels[item]][0]))
        self.nb.data, self.nb.MainControlPage.RMS = calc_CTRs(self.nb.parameter,self.nb.parameter_usage, self.nb.data, self.nb.cell,\
                               self.nb.surface, self.nb.NLayers, database, self.nb.g_inv, self.nb.MainControlPage.Rod_weight, self.nb.rigid_bodies, \
                                                              self.nb.MainControlPage.UBW_flag, self.nb.MainControlPage.use_BVC, self.nb.MainControlPage.BVclusters,\
                                                              self.nb.MainControlPage.RMS_flag, self.nb.MainControlPage.use_lay_el, self.nb.MainControlPage.el)
        self.nb.MainControlPage.Figure1 = plot_rods(self.nb.MainControlPage.Figure1,self.nb.data, self.nb.MainControlPage.plotdims, self.nb.MainControlPage.doplotbulk, self.nb.MainControlPage.doplotsurf, self.nb.MainControlPage.doplotrough,\
                                                 self.nb.MainControlPage.doplotwater, self.nb.MainControlPage.RMS)
        self.nb.MainControlPage.Figure1.canvas.draw()
        self.nb.MainControlPage.Figure2 = plot_edensity(self.nb.MainControlPage.Figure2, self.nb.surface, self.nb.parameter, self.nb.parameter_usage, self.nb.cell, database, self.nb.rigid_bodies, self.nb.MainControlPage.UBW_flag, self.nb.MainControlPage.use_lay_el, self.nb.ResonantDataPage.resel)
        self.nb.MainControlPage.Figure2.canvas.draw()
        check_vibes(self.nb.surface,self.nb.parameter, self.nb.parameter_usage)
    def togglestep(self,event):
        item = event.GetId()-5*len(self.nb.param_labels)-2000
        try:
            self.togglesteps[item] = float(event.GetString())
        except ValueError:
            pass
    def toggleplus(self, event):
        item = event.GetId()-6*len(self.nb.param_labels)-2000
        step = self.togglesteps[item]
        self.nb.parameter[self.nb.param_labels[item]][0] = self.nb.parameter[self.nb.param_labels[item]][0] + step
        self.control1[item].SetValue(str(self.nb.parameter[self.nb.param_labels[item]][0]))
        self.nb.data, self.nb.MainControlPage.RMS = calc_CTRs(self.nb.parameter,self.nb.parameter_usage, self.nb.data, self.nb.cell,\
                               self.nb.surface, self.nb.NLayers, database, self.nb.g_inv, self.nb.MainControlPage.Rod_weight, self.nb.rigid_bodies, \
                                                              self.nb.MainControlPage.UBW_flag, self.nb.MainControlPage.use_BVC, self.nb.MainControlPage.BVclusters,\
                                                              self.nb.MainControlPage.RMS_flag, self.nb.MainControlPage.use_lay_el, self.nb.MainControlPage.el)
        self.nb.MainControlPage.Figure1 = plot_rods(self.nb.MainControlPage.Figure1,self.nb.data, self.nb.MainControlPage.plotdims, self.nb.MainControlPage.doplotbulk, self.nb.MainControlPage.doplotsurf, self.nb.MainControlPage.doplotrough,\
                                                 self.nb.MainControlPage.doplotwater, self.nb.MainControlPage.RMS)
        self.nb.MainControlPage.Figure1.canvas.draw()
        self.nb.MainControlPage.Figure2 = plot_edensity(self.nb.MainControlPage.Figure2,self.nb.surface, self.nb.parameter, self.nb.parameter_usage, self.nb.cell, database, self.nb.rigid_bodies, self.nb.MainControlPage.UBW_flag, self.nb.MainControlPage.use_lay_el,self.nb.ResonantDataPage.resel)
        self.nb.MainControlPage.Figure2.canvas.draw()
        check_vibes(self.nb.surface,self.nb.parameter, self.nb.parameter_usage)
############################################################################################################
############################################################################################################
def start():
    frame = wxCtrFitFrame(parent = None, title = u"\u03c0-surf,  Python Interface StrUcture ReFinement", size = (800,750))
    return frame


