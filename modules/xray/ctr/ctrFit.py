from ctrFitcalcs import *
import wx
import os
"""
Functions and classes used to build the CtrFit GUI
by Frank Heberling (Frank.Heberling@kit.edu)
"""
class wxCtrFitFrame(wx.Frame):
    def __init__(self, parent, title, size):
        wx.Frame.__init__(self, parent, title= title, size=(800,800))

        self.filename1 = ''
        self.filename2 = ''
        self.filename3 = ''
        self.filename4 = ''
        self.filename5 = ''
        
        # A status bar
        self.CreateStatusBar()

        # Setting up the ReadFilesmenu.
        filemenu= wx.Menu()
        menuReadData = filemenu.Append(wx.ID_ANY,"&Read Datafile"," Read in a Data file")
        menuReadBulk = filemenu.Append(wx.ID_ANY,"&Read Bulkfile"," Read in a Bulk file")
        menuReadSurface = filemenu.Append(wx.ID_ANY,"&Read Surfacefile"," Read in a Surface file")
        menuReadParameter = filemenu.Append(wx.ID_ANY,"&Read Parameterfile"," Read in a Parameter file")
        menuReadRigidbody = filemenu.Append(wx.ID_ANY,"&Read Rigidbodyfile"," Read in a Rigid Body file")
        filemenu.AppendSeparator()
        menuExit = filemenu.Append(wx.ID_EXIT,"E&xit"," Terminate the program")

        #Setting up the WriteFiles
        writemenu = wx.Menu()
        menuWriteCif = writemenu.Append(wx.ID_ANY, "&Write .cif file", " Write surface structure to a .cif file")
        menuWritePar = writemenu.Append(wx.ID_ANY, "&Write Parameter file", " Write parameters to a .par file")
        menuWriteSurf = writemenu.Append(wx.ID_ANY, "&Write Surface file", " Write surface in fractional coordinates to a .sur file")
        menuWriteBulk = writemenu.Append(wx.ID_ANY, "&Write Bulk file", " Write bulk structure to a .bul file")
      
        # Creating the menubar.
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&Read Files") # Adding the "filemenu" to the MenuBar
        menuBar.Append(writemenu,"&Write Files")
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        
        # define events
        self.Bind(wx.EVT_MENU, self.OnReadData, menuReadData)
        self.Bind(wx.EVT_MENU, self.OnReadBulk, menuReadBulk)
        self.Bind(wx.EVT_MENU, self.OnReadSurface, menuReadSurface)
        self.Bind(wx.EVT_MENU, self.OnReadParameter, menuReadParameter)
        self.Bind(wx.EVT_MENU, self.OnReadRigidbody, menuReadRigidbody)
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)

        self.Bind(wx.EVT_MENU, self.OnWriteCif, menuWriteCif)
        self.Bind(wx.EVT_MENU, self.OnWritePar, menuWritePar)
        self.Bind(wx.EVT_MENU, self.OnWriteSurf, menuWriteSurf)
        self.Bind(wx.EVT_MENU, self.OnWriteBulk, menuWriteBulk)

        self.nb = CtrNotebook(self)

        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.Show(True)
        

    def OnReadData(self,e):
        """ Read in data"""
        dirname1 = ''
        dlg = wx.FileDialog(self, "Choose a data file", dirname1, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename1 = dlg.GetFilename()
            dirname1 = dlg.GetDirectory()
            self.nb.data = read_data(dirname1+'/'+self.filename1)
            self.nb.MainControlPage.datafile.SetValue(self.filename1)
            for i in range(len(self.nb.data)):
                wx.StaticText(self.nb.MainControlPage, label = (str(int(self.nb.data[i].H))+' '+str(int(self.nb.data[i].K))+' L'), pos=(350,25*i+67), size=(40,20))
                self.nb.MainControlPage.rodweight.append(wx.TextCtrl(self.nb.MainControlPage,1000+i, pos=(400,25*i+65), size=(30,20)))
                self.nb.MainControlPage.rodweight[i].SetValue('1')
                self.Bind(wx.EVT_TEXT, self.nb.MainControlPage.setrodweight, self.nb.MainControlPage.rodweight[i])
                self.nb.MainControlPage.Rod_weight.append(1)
        dlg.Destroy()

    def OnReadBulk(self,e):
        """ Read in bulk file"""
        dirname2 = ''
        dlg = wx.FileDialog(self, "Choose a bulk file", dirname2, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename2 = dlg.GetFilename()
            dirname2 = dlg.GetDirectory()
            self.nb.bulk, self.nb.cell, self.nb.NLayers = read_bulk(dirname2+'/'+self.filename2)
            self.nb.MainControlPage.bulkfile.SetValue(self.filename2)
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
        dlg.Destroy()

    def OnReadSurface(self,e):
        """ Read in surface file"""
        dirname3 = ''
        dlg = wx.FileDialog(self, "Choose a surface file", dirname3, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename3 = dlg.GetFilename()
            dirname3 = dlg.GetDirectory()
            self.nb.surface, self.nb.parameter_usage = read_surface(dirname3+'/'+self.filename3)
            self.nb.MainControlPage.surfacefile.SetValue(self.filename3)
        dlg.Destroy()

    def OnReadParameter(self,e):
        """ Read in parameter file"""
        dirname4 = ''
        dlg = wx.FileDialog(self, "Choose a parameter file", dirname4, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename4 = dlg.GetFilename()
            dirname4 = dlg.GetDirectory()
            self.nb.parameter, self.nb.param_labels = read_parameters(dirname4+'/'+self.filename4)
            self.nb.MainControlPage.parameterfile.SetValue(self.filename4)

            for i in range(len(self.nb.param_labels)):
                wx.StaticText(self.nb.ParameterPage, label = self.nb.param_labels[i], pos=(20, 23*(i+1)+45), size=(120, 20))
                
                control1_tmp = wx.TextCtrl(self.nb.ParameterPage,2000+i, pos=(150, 23*(i+1)+40), size=(120,20))
                self.nb.ParameterPage.control1.append(control1_tmp)
                self.nb.ParameterPage.control1[i].SetValue(str(round(self.nb.parameter[self.nb.param_labels[i]][0], 12)))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparvalue, self.nb.ParameterPage.control1[i])
                
                control2_tmp = (wx.TextCtrl(self.nb.ParameterPage,2000+i+len(self.nb.param_labels), pos=(280, 23*(i+1)+40), size=(80,20)))
                self.nb.ParameterPage.control2.append(control2_tmp)
                self.nb.ParameterPage.control2[i].SetValue(str(self.nb.parameter[self.nb.param_labels[i]][1]))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparmin, self.nb.ParameterPage.control2[i])
                
                control3_tmp = (wx.TextCtrl(self.nb.ParameterPage,2000+i+2*len(self.nb.param_labels), pos=(370, 23*(i+1)+40), size=(80,20)))
                self.nb.ParameterPage.control3.append(control3_tmp)
                self.nb.ParameterPage.control3[i].SetValue(str(self.nb.parameter[self.nb.param_labels[i]][2]))
                self.Bind(wx.EVT_TEXT, self.nb.ParameterPage.editparmax, self.nb.ParameterPage.control3[i])
                                
                control4_tmp = (wx.CheckBox(self.nb.ParameterPage,2000+i+3*len(self.nb.param_labels), label = '', pos = (460, 23*(i+1)+45)))
                self.nb.ParameterPage.control4.append(control4_tmp)
                self.nb.ParameterPage.control4[i].SetValue(self.nb.parameter[self.nb.param_labels[i]][3])
                wx.EVT_CHECKBOX(self.nb.ParameterPage, self.nb.ParameterPage.control4[i].GetId(), self.nb.ParameterPage.editparstate)

            self.nb.ParameterPage.SetScrollbars(0, 10, 0, int((len(self.nb.param_labels)+4)*2.3)+1)
        dlg.Destroy()

    def OnReadRigidbody(self,e):
        """ Read in rigid body file"""
        dirname5 = ''
        dlg = wx.FileDialog(self, "Choose a rigid body file", dirname5, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename5 = dlg.GetFilename()
            dirname5 = dlg.GetDirectory()
            self.nb.rigid_bodies = read_rigid_bodies(dirname5+'/'+self.filename5)
            self.nb.MainControlPage.rigidbodyfile.SetValue(self.filename5)
        dlg.Destroy()

    def OnExit(self,e):
        self.Close(True)  # Close the frame.

    def OnWriteCif(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Write surface structure to .cif file", dirname, ".cif", "*.cif", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            write_cif(self.nb.cell, self.nb.surface, self.nb.parameter,self.nb.parameter_usage, self.nb.rigid_bodies, self.nb.MainControlPage.use_bulk_water, filename)
        dlg.Destroy()

    def OnWritePar(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Write parameters to .par file", dirname, ".par", "*.par", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            write_par(self.nb.parameter, self.nb.param_labels, filename)
        dlg.Destroy()
        
    def OnWriteSurf(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Write surface structure to .sur file", dirname, ".sur", "*.sur", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            write_surface(self.nb.cell, self.nb.surface, self.nb.parameter, self.nb.parameter_usage, self.nb.rigid_bodies, self.nb.MainControlPage.use_bulk_water, filename)
        dlg.Destroy()
        
    def OnWriteBulk(self,e):
        dirname = ''
        dlg = wx.FileDialog(self, "Write bulk structure to .bul file", dirname, ".bul", "*.bul", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            os.chdir(dirname)
            write_bulk(self.nb.bulk, self.nb.parameter,self.nb.MainControlPage.use_bulk_water, filename)
        dlg.Destroy()
############################################################################################################
class CtrNotebook(wx.Notebook):
    def __init__(self, parent):
        wx.Notebook.__init__(self,parent)
        self.MainControlPage = MainControlPanel(self)
        self.ParameterPage = ParameterPanel(self)

        self.AddPage(self.MainControlPage, " Main Controls " )
        self.AddPage(self.ParameterPage, " Parameters ")

        self.parameter_usage = []
        self.surface = []
        self.bulk = []
        self.rigid_bodies = []
        self.data= []
        self.NLayers = 1
        self.cell = []
        self.g_inv = []

        self.param_best = {}
        self.parameter = {}
        self.param_labels = []
        

        
############################################################################################################
class MainControlPanel(wx.Panel):
    def __init__(self,parent):
        wx.Panel.__init__(self,parent)
        self.doplotbulk = False
        self.doplotsurf = False
        self.doplotrough = False
        self.doplotwater = False
        self.doplotRMStrack = False
        self.plotdims = [2,3]

        self.rodweight = []# List of TextControls to read rod weights
        self.sim_an_params = [50, 20, 0.7, 10, 0.01, 500000, False]
        self.Rod_weight = [] #List of Rod weights
        self.RMS = -1
        self.RMS_best = -1
        self.use_bulk_water = True

        self.nb = self.GetParent()

        #Panel Headings
        wx.StaticText(self, label = 'File Information: ', pos=(80, 12), size=(100, 20))
        wx.StaticText(self, label = 'Data Information: ', pos=(350, 12), size=(100, 20))
        wx.StaticText(self, label = 'Fit options: ', pos=(610, 12), size=(100, 20))

        wx.StaticText(self, label = 'Rods:', pos=(350,42), size=(40,20))
        wx.StaticText(self, label = 'weight', pos=(400,42), size=(30,20))

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

        #CTR plotting
        wx.StaticText(self, label = 'Plotting Options: ', pos=(80, 272), size=(100, 20))
        wx.StaticText(self, label = 'CTR plot options (Fig. 1): ', pos=(20, 297), size=(200, 20))

        wx.StaticText(self, label = 'CTR plot dimensions', pos=(40, 322), size=(110, 20))
        self.getplotdims = wx.TextCtrl(self, pos=(160, 320), size=(40,-1))
        self.getplotdims.SetValue(str(self.plotdims[0])+' '+str(self.plotdims[1]))
        self.Bind(wx.EVT_TEXT, self.setplotdims, self.getplotdims)

        self.plot_bulk = wx.CheckBox(self, label =  '  plot bulk (green)', pos = (40, 347))
        self.plot_bulk.SetValue(False)
        wx.EVT_CHECKBOX(self, self.plot_bulk.GetId(), self.setplotbulk)
        
        self.plot_surf = wx.CheckBox(self, label =  '  plot surface (red)', pos = (40, 372))
        self.plot_surf.SetValue(False)
        wx.EVT_CHECKBOX(self, self.plot_surf.GetId(), self.setplotsurf)
        
        self.plot_rough = wx.CheckBox(self, label = '  plot roughness (cyan)', pos = (40, 397))
        self.plot_rough.SetValue(False)
        wx.EVT_CHECKBOX(self, self.plot_rough.GetId(), self.setplotrough)
        
        self.plot_water = wx.CheckBox(self, label = '  plot water (magenta)', pos = (40, 422))
        self.plot_water.SetValue(False)
        wx.EVT_CHECKBOX(self, self.plot_water.GetId(), self.setplotwater)

        self.button = wx.Button(self, label = 'Plot CTRs', pos =(170,415))
        self.Bind(wx.EVT_BUTTON, self.OnClick, self.button)

        #edensity plotting
        wx.StaticText(self, label = 'plot e-density', pos=(20, 460), size=(140, 20))
        wx.StaticText(self, label = '(Fig. 2, current surface)', pos=(20, 480), size=(140, 20))
        self.plotedens = wx.Button(self, label = 'Plot edens', pos =(170,470))
        self.Bind(wx.EVT_BUTTON, self.OnPlotedens, self.plotedens)

        #Fitting optio use bulk water
        self.douse_bulk_water = wx.CheckBox(self, label = '  use bulk water', pos = (550, 37))
        self.douse_bulk_water.SetValue(True)
        wx.EVT_CHECKBOX(self, self.douse_bulk_water.GetId(), self.set_bulk_water)

        #Simulated annealing parameters and options
        wx.StaticText(self, label = 'Simulated Annealing:  ', pos=(590, 67), size=(140, 20))

        wx.StaticText(self, label = 'start Temperature:    ', pos=(550, 92), size=(140, 20))
        self.Tstart = wx.TextCtrl(self, pos=(700,90), size=(60,20))
        self.Tstart.SetValue(str(self.sim_an_params[0]))
        self.Bind(wx.EVT_TEXT, self.setTstart, self.Tstart)

        wx.StaticText(self, label = 'end Temperature:      ', pos=(550, 117), size=(140, 20))
        self.Tend = wx.TextCtrl(self, pos=(700,115), size=(60,20))
        self.Tend.SetValue(str(self.sim_an_params[1]))
        self.Bind(wx.EVT_TEXT, self.setTend, self.Tend)

        wx.StaticText(self, label = 'cooling factor:       ', pos=(550, 142), size=(140, 20))
        self.cool = wx.TextCtrl(self, pos=(700,140), size=(60,20))
        self.cool.SetValue(str(self.sim_an_params[2]))
        self.Bind(wx.EVT_TEXT, self.setcool, self.cool)

        wx.StaticText(self, label = 'iterations per cycle: ', pos=(550, 167), size=(140, 20))
        self.iter = wx.TextCtrl(self, pos=(700,165), size=(60,20))
        self.iter.SetValue(str(self.sim_an_params[3]))
        self.Bind(wx.EVT_TEXT, self.setiter, self.iter)

        wx.StaticText(self, label = 'Boltzmann weight:     ', pos=(550, 192), size=(140, 20))
        self.Boltz = wx.TextCtrl(self, pos=(700,190), size=(60,20))
        self.Boltz.SetValue(str(self.sim_an_params[5]))
        self.Bind(wx.EVT_TEXT, self.setBoltz, self.Boltz)
        
        self.random_pars = wx.CheckBox(self, label = '  start Fit with random parameters', pos = (550, 217))
        self.random_pars.SetValue(False)
        wx.EVT_CHECKBOX(self, self.random_pars.GetId(), self.setrandom_pars)
        
        self.plot_RMS_track = wx.CheckBox(self, label = '  plot RMS track (Fig. 3, after Fit)', pos = (550, 242))
        self.plot_RMS_track.SetValue(False)
        wx.EVT_CHECKBOX(self, self.plot_RMS_track.GetId(), self.setplotRMStrack)

        self.SimAn1button = wx.Button(self, label = 'Simulated Annealing 1', pos =(550,270), size=(210,40))
        self.Bind(wx.EVT_BUTTON, self.OnClickSimAn1, self.SimAn1button)

        self.SimAn2button = wx.Button(self, label = 'Simulated Annealing 2', pos =(550,315), size=(210,40))
        self.Bind(wx.EVT_BUTTON, self.OnClickSimAn2, self.SimAn2button)

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
        self.nb.data, self.RMS = calc_CTRs(self.nb.parameter,self.nb.parameter_usage, self.nb.data, self.nb.cell, self.nb.bulk, \
                               self.nb.surface, self.nb.NLayers, database, self.nb.g_inv, self.Rod_weight, self.nb.rigid_bodies, self.use_bulk_water)
        plot_rods(self.nb.data, self.plotdims, self.doplotbulk, self.doplotsurf, self.doplotrough,\
                                                 self.doplotwater, self.RMS)

    def OnPlotedens(self, e):
        plot_edensity(self.nb.surface, self.nb.parameter, self.nb.parameter_usage, self.nb.cell, database, self.nb.rigid_bodies, self.use_bulk_water)

    def set_bulk_water(self,e): self.use_bulk_water = self.douse_bulk_water.GetValue()

    def setTstart(self,event):
        if (event.GetString() == '') or (event.GetString() == '-'):
            None
        else:
            a = float(event.GetString())
            if a <= self.sim_an_params[1]:
                print 'Tstart must be > Tend'
            else:
                self.sim_an_params[0] = a
            
    def setTend(self,event):
        if (event.GetString() == '') or (event.GetString() == '-'):
            None
        else:
            a = float(event.GetString())
            if (a >= self.sim_an_params[0]) or (a <=0):
                print 'Tend must be < Tstart and > 0'
            else:
                self.sim_an_params[1] = a
            
    def setcool(self,event):
        if (event.GetString() == '') or (event.GetString() == '-'):
            None
        else:
            a = float(event.GetString())
            if (a <=0) or (a >=1):
                print 'cooling factor must be > 0 and < 1'
            else:
                self.sim_an_params[2] = a
            
    def setiter(self,event):
        if (event.GetString() == '') or (event.GetString() == '-'):
            None
        else:
            a = int(event.GetString())
            if a <=0:
                print 'number of iterations must be positive'
            else:
                self.sim_an_params[3] = a
    
    def setBoltz(self,event):
        if (event.GetString() == '') or (event.GetString() == '-'):
            None
        else:
            a = float(event.GetString())
            if a < 1000:
                print 'Boltzmann weight must be > 1000'
            else:
                self.sim_an_params[5] = a
            
    def setrandom_pars(self,e): self.sim_an_params[6] = self.random_pars.GetValue()
    def setplotRMStrack(self,e): self.doplotRMStrack = self.plot_RMS_track.GetValue()

    def OnClickSimAn1(self,e):
        flag = check_model_consistency(self.nb.param_labels, self.nb.parameter, self.nb.parameter_usage, self.nb.rigid_bodies, self.use_bulk_water)
        if flag:
            self.nb.param_best = {}
            self.RMS_best = -1
            self.nb.data, self.nb.param_best, self.RMS_best = simulated_annealing01(self.nb.data, self.nb.cell, self.nb.NLayers, self.nb.bulk, self.nb.surface,\
                    database, self.Rod_weight, self.sim_an_params, self.nb.parameter, self.nb.parameter_usage, self.doplotRMStrack, self.nb.rigid_bodies, self.use_bulk_water)
            plot_rods(self.nb.data, self.plotdims, self.doplotbulk, self.doplotsurf, self.doplotrough, self.doplotwater, self.RMS_best)
            dlg = wx.MessageDialog(self, "Keep refined parameters ?","", wx.YES_NO | wx.STAY_ON_TOP)
            if dlg.ShowModal() == wx.ID_OK:
                for i in range(len(self.nb.param_labels)):
                    self.nb.parameter[self.nb.param_labels[i]][0] = self.nb.param_best[self.nb.param_labels[i]][0]
                    self.nb.ParameterPage.OnClick(True)
            dlg.Destroy()
            
    def OnClickSimAn2(self,e):
        flag = check_model_consistency(self.nb.param_labels, self.nb.parameter, self.nb.parameter_usage, self.nb.rigid_bodies, self.use_bulk_water)
        if flag:
            self.nb.param_best = {}
            self.RMS_best = -1
            self.nb.data, self.nb.param_best, self.RMS_best = simulated_annealing02(self.nb.data, self.nb.cell, self.nb.NLayers, self.nb.bulk, self.nb.surface,\
                    database, self.Rod_weight, self.sim_an_params, self.nb.parameter, self.nb.parameter_usage, self.doplotRMStrack, self.nb.rigid_bodies, self.use_bulk_water)
            plot_rods(self.nb.data, self.plotdims, self.doplotbulk, self.doplotsurf, self.doplotrough, self.doplotwater, self.RMS_best)
            dlg = wx.MessageDialog(self, "Keep refined parameters ?","", wx.YES_NO | wx.STAY_ON_TOP)
            if dlg.ShowModal() == wx.ID_OK:
                for i in range(len(self.nb.param_labels)):
                    self.nb.parameter[self.nb.param_labels[i]][0] = self.nb.param_best[self.nb.param_labels[i]][0]
                    self.nb.ParameterPage.OnClick(True)
            dlg.Destroy()    

##########################################################################################################################
class ParameterPanel(wx.ScrolledWindow):
    def __init__(self,parent):
        wx.ScrolledWindow.__init__(self,parent)

        self.control1 = []
        self.control2 = []
        self.control3 = []
        self.control4 = []

        wx.StaticText(self, label = 'value', pos=(170, 40), size=(100, 20))
        wx.StaticText(self, label = 'min', pos=(300, 40), size=(70, 20))
        wx.StaticText(self, label = 'max', pos=(390, 40), size=(70, 20))
        wx.StaticText(self, label = 'refine', pos=(455, 40), size=(60, 20))

        self.button = wx.Button(self, label = 'Update', pos =(20,20))
        self.Bind(wx.EVT_BUTTON, self.OnClick, self.button)

        self.nb = self.GetParent()
        
    def OnClick(self,e):
        for i in range(len(self.nb.param_labels)):
            self.control1[i].SetValue(str(round(self.nb.parameter[self.nb.param_labels[i]][0], 12)))
        

    def editparvalue(self,event):
        item = event.GetId()-2000
        if (event.GetString() == '') or (event.GetString() == '-'):
            None
        else:
            self.nb.parameter[self.nb.param_labels[item]][0] = float(event.GetString())
        
    def editparmin(self,event):
        item = event.GetId()-len(self.nb.param_labels)-2000
        if (event.GetString() == '') or (event.GetString() == '-'):
            None
        else:
            self.nb.parameter[self.nb.param_labels[item]][1] = float(event.GetString())
        
    def editparmax(self,event):
        item = event.GetId()-2*len(self.nb.param_labels)-2000
        if (event.GetString() == '') or (event.GetString() == '-'):
            None
        else:
            self.nb.parameter[self.nb.param_labels[item]][2] = float(event.GetString())
        
    def editparstate(self,event):
        item = event.GetId()- 3*len(self.nb.param_labels)-2000
        self.nb.parameter[self.nb.param_labels[item]][3] = self.control4[item].GetValue()

############################################################################################################
############################################################################################################
def start_ctr_fitting():
    frame = wxCtrFitFrame(parent = None, title = " CtrFit", size = (800,800))


