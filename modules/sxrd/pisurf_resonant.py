"""
Functions and classes used to build the resonant data handling part of the
pi-surf GUI

Authors/modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""
###############################################################################
import numpy as Num
from pylab import *
import wx

from tdl.modules.sxrd.ctrfitcalcs import param_unfold,RB_update
from tdl.modules.sxrd.resonant_simplex import res_simplex, res_param_statistics
from tdl.modules.sxrd.fourierframe import createfourierframe
from tdl.modules.sxrd.resonant_calcs import *

from tdl.modules.xtab.atomic import f0data as database
#################################################################################
class ResonantDataPanel(wx.ScrolledWindow):
    def __init__(self,parent):
        wx.ScrolledWindow.__init__(self,parent)
        
        self.nb = self.GetParent()
        self.resel = ''
        self.plot_norm = True
        self.Qmax = 6.

        self.allrasd = RasdList()
        self.allrasd.g_inv = self.nb.g_inv
        self.allrasd.cell = self.nb.cell

        self.fframe = None
        self.fig7 = None
        self.fig8 = None

        self.simplex_params = [1.0,0.5,2.0,0.5,1e-6,10000]

        self.datacontrol0 = []
        self.datacontrol1 = []
        self.datacontrol2 = []
        self.datacontrol3 = []
        self.datacontrol4 = []
        self.datacontrol5 = []
        self.datacontrol6 = []
        self.datacontrol7 = []
        self.datacontrol8 = []

        wx.StaticText(self, label = ' Resonant Data ', pos=(290, 17), size=(120, 20))
        wx.StaticText(self, label = ' AR    +/- dAr', pos=(420, 17), size=(80, 20))
        wx.StaticText(self, label = ' PR    +/- dPr', pos=(515, 17), size=(80, 20))
        wx.StaticText(self, label = 'UIF', pos=(610, 17), size=(20, 20))
        wx.StaticText(self, label = 'UIR', pos=(630, 17), size=(20, 20))
        wx.StaticText(self, label = ' AR_m ', pos=(655, 17), size=(60, 20))
        wx.StaticText(self, label = ' PR_m ', pos=(700, 17), size=(60, 20))

        wx.StaticLine(self, pos = (270,0), size = (5,650), style = wx.LI_VERTICAL)

        wx.StaticText(self, label = 'Resonant Element = ', pos=(20, 12), size=(100, 20))
        self.getresel = wx.TextCtrl(self, pos=(130,10), size=(30,20))
        self.Bind(wx.EVT_TEXT, self.setresel, self.getresel)

        wx.StaticText(self, label = 'Z = ', pos=(170, 12), size=(20, 20))
        self.showZ = wx.TextCtrl(self, pos=(200,10), size=(50,20), style= wx.TE_READONLY)

        wx.StaticText(self, label = 'Edge Energy (eV) = ', pos=(20, 37), size=(100, 20))
        self.getE0 = wx.TextCtrl(self, pos=(130,35), size=(120,20))
        self.Bind(wx.EVT_TEXT, self.setE0, self.getE0)

        wx.StaticText(self, label = 'f1f2 file = ', pos=(20, 62), size=(100, 20))
        self.f1f2file = wx.TextCtrl(self, pos=(130,60), size=(120,20), style= wx.TE_READONLY)

        self.getplotnorm = wx.CheckBox(self, label = '    plot data backgroud corrected', pos = (20, 87))
        self.getplotnorm.SetValue(self.plot_norm)
        wx.EVT_CHECKBOX(self, self.getplotnorm.GetId(), self.setplotnorm)

        self.updatebutton = wx.Button(self, label = 'Update non resonant structure factors', pos =(20,110), size=(230,25))
        self.Bind(wx.EVT_BUTTON, self.OnClickUpdate, self.updatebutton)
        
        wx.StaticLine(self, pos = (0,140), size = (270,5), style = wx.LI_HORIZONTAL)
        
        self.getdoabscorr = wx.CheckBox(self, label = '    correct data for absorption effects', pos = (20, 152))
        self.getdoabscorr.SetValue(self.allrasd.do_abs_corr)
        wx.EVT_CHECKBOX(self, self.getdoabscorr.GetId(), self.setdoabscorr)

        wx.StaticText(self, label = 'water film thickness (m) ', pos=(20, 177), size=(150, 20))
        self.getdwater = wx.TextCtrl(self, pos=(180,175), size=(70,20))
        self.Bind(wx.EVT_TEXT, self.setdwater, self.getdwater)

        wx.StaticText(self, label = 'Res. Element conc. (mol/m**3) ', pos=(20, 202), size=(150, 20))
        self.getconc = wx.TextCtrl(self, pos=(180,200), size=(70,20))
        self.Bind(wx.EVT_TEXT, self.setconc, self.getconc)

        wx.StaticText(self, label = 'Kapton film thickness (m) ', pos=(20, 227), size=(150, 20))
        self.getdkapton = wx.TextCtrl(self, pos=(180,225), size=(70,20))
        self.Bind(wx.EVT_TEXT, self.setdkapton, self.getdkapton)

        self.abscorrbutton = wx.Button(self, label = 'Calculate absorption correction', pos =(20,250), size=(230,25))
        self.Bind(wx.EVT_BUTTON, self.OnClickabscorr, self.abscorrbutton)

        wx.StaticLine(self, pos = (0,280), size = (270,5), style = wx.LI_HORIZONTAL)

        wx.StaticText(self, label = 'cells in x ', pos=(20, 292), size=(60, 20))
        self.getxf = wx.TextCtrl(self, pos=(90,290), size=(40,20))
        self.getxf.SetValue(str(self.allrasd.xf))
        self.Bind(wx.EVT_TEXT, self.setxf, self.getxf)

        wx.StaticText(self, label = '          y ', pos=(20, 317), size=(60, 20))
        self.getyf = wx.TextCtrl(self, pos=(90,315), size=(40,20))
        self.getyf.SetValue(str(self.allrasd.yf))
        self.Bind(wx.EVT_TEXT, self.setyf, self.getyf)

        wx.StaticText(self, label = '          z ', pos=(20, 342), size=(60, 20))
        self.getzf = wx.TextCtrl(self, pos=(90,340), size=(40,20))
        self.getzf.SetValue(str(self.allrasd.zf))
        self.Bind(wx.EVT_TEXT, self.setzf, self.getzf)

        wx.StaticText(self, label = '       zmin ', pos=(20, 367), size=(60, 20))
        self.getzmin = wx.TextCtrl(self, pos=(90,365), size=(40,20))
        self.getzmin.SetValue(str(self.allrasd.zmin))
        self.Bind(wx.EVT_TEXT, self.setzmin, self.getzmin)

        wx.StaticText(self, label = 'pixels in x ', pos=(140, 292), size=(60, 20))
        self.getan = wx.TextCtrl(self, pos=(210,290), size=(40,20))
        self.getan.SetValue(str(self.allrasd.an))
        self.Bind(wx.EVT_TEXT, self.setan, self.getan)

        wx.StaticText(self, label = '           y ', pos=(140, 317), size=(60, 20))
        self.getbn = wx.TextCtrl(self, pos=(210,315), size=(40,20))
        self.getbn.SetValue(str(self.allrasd.bn))
        self.Bind(wx.EVT_TEXT, self.setbn, self.getbn)

        wx.StaticText(self, label = '           z ', pos=(140, 342), size=(60, 20))
        self.getcn = wx.TextCtrl(self, pos=(210,340), size=(40,20))
        self.getcn.SetValue(str(self.allrasd.cn))
        self.Bind(wx.EVT_TEXT, self.setcn, self.getcn)

        self.FourierSynthButton = wx.Button(self, label = 'Calculate Fourier Synthesis', pos =(20,390), size=(230,25))
        self.Bind(wx.EVT_BUTTON, self.OnClickFSB, self.FourierSynthButton)
        
        wx.StaticLine(self, pos = (0,420), size = (270,5), style = wx.LI_HORIZONTAL)

        self.getRD = wx.CheckBox(self, label = ' check to refine data, uncheck for AR and PR', pos = (20, 435))
        self.getRD.SetValue(self.allrasd.Refine_Data)
        wx.EVT_CHECKBOX(self, self.getRD.GetId(), self.setRD)
		
        self.ulem = wx.CheckBox(self, label = ' check to use layered resonant element model', pos = (20, 460))
        self.ulem.SetValue(self.nb.MainControlPage.use_lay_el)
        wx.EVT_CHECKBOX(self, self.ulem.GetId(), self.setulem)
		
        wx.StaticText(self, label = 'Qmax ', pos=(20, 487), size=(60, 20))
        self.getqmax = wx.TextCtrl(self, pos=(90,485), size=(40,20))
        self.getqmax.SetValue(str(self.Qmax))
        self.Bind(wx.EVT_TEXT, self.setqmax, self.getqmax)
        self.PlotArPrQButton = wx.Button(self, label = 'Ar(Q) Pr(Q)', pos =(140,485), size=(110,25))
        self.Bind(wx.EVT_BUTTON, self.OnClickPlotArPrQ, self.PlotArPrQButton)
		
        self.RefineButton = wx.Button(self, label = 'Run Refinement', pos =(20,520), size=(230,80))
        self.Bind(wx.EVT_BUTTON, self.OnClickRefine, self.RefineButton)

        self.ResStatisticsButton = wx.Button(self, label = 'Calculate resonant parameter statistics', pos =(20,610), size=(230,30))
        self.Bind(wx.EVT_BUTTON, self.OnClickResStatistics, self.ResStatisticsButton)
        
        wx.StaticLine(self, pos = (0,645), size = (270,5), style = wx.LI_HORIZONTAL)
        
    def setE0(self,event):
        try:
            self.allrasd.E0 = float(event.GetString())
        except ValueError:
            pass

    def setresel(self,event):
        if event.GetString() == '':
            pass
        try:
            self.resel = str(event.GetString())
            if str.lower(self.resel) not in database.keys():
                print 'unknown element: '+self.resel
            else:
                f_par = database[str.lower(self.resel)]
                f = (f_par[0] + f_par[2] + f_par[4] + f_par[6] + f_par[8])
                self.showZ.SetValue(str(round(f,0)))
                self.allrasd.ZR = round(f,0)
                self.nb.MainControlPage.el = str.lower(self.resel)
                if str.lower(self.resel) not in self.nb.runningDB.keys():
                    self.nb.runningDB[str.lower(self.resel)] = database[str.lower(self.resel)]
        except ValueError:
            pass

    def setplotnorm(self,e): self.plot_norm = self.getplotnorm.GetValue()

    def OnClickUpdate(self,e):
        global_parms, surface_new = param_unfold(self.nb.parameter,self.nb.parameter_usage,self.nb.surface,\
                                                 self.nb.MainControlPage.UBW_flag, self.nb.MainControlPage.use_lay_el)
        bulk = self.nb.bulk 
        surface = RB_update(self.nb.rigid_bodies, surface_new, self.nb.parameter, self.allrasd.cell)
        for Rasd in self.allrasd.list:
            Rasd.re_FNR ,Rasd.im_FNR = calcFNR(Rasd.Q,global_parms, self.allrasd.cell, bulk, surface, database,\
                                               self.allrasd.g_inv, self.nb.MainControlPage.UBW_flag,\
                                               self.nb.MainControlPage.use_lay_el, str.lower(self.resel))

        for item in range(len(self.allrasd.list)):
            self.allrasd.list[item] = RASD_Fourier(self.allrasd, item)
            self.datacontrol1[item].SetValue(str(round(self.allrasd.list[item].AR,3)))
            self.datacontrol2[item].SetValue(str(round(self.allrasd.list[item].AR_err,3)))
            self.datacontrol3[item].SetValue(str(round(self.allrasd.list[item].PR,3)))
            self.datacontrol8[item].SetValue(str(round(self.allrasd.list[item].PR_err,3)))
            
    
    def ClickDataButton(self,e):
        item = e.GetId()-10000
        self.allrasd.list[item] = RASD_Fourier(self.allrasd, item)
        self.datacontrol1[item].SetValue(str(round(self.allrasd.list[item].AR,3)))
        self.datacontrol2[item].SetValue(str(round(self.allrasd.list[item].AR_err,3)))
        self.datacontrol3[item].SetValue(str(round(self.allrasd.list[item].PR,3)))
        self.datacontrol8[item].SetValue(str(round(self.allrasd.list[item].PR_err,3)))
        self.allrasd.list[item].plot(norm = self.plot_norm, fig = 4)

    def editPR(self,e):
        item = e.GetId()-10000-3*len(self.allrasd.list)
        try:
            self.allrasd.list[item].PR = float(e.GetString())
        except ValueError:
            pass

    def editUIF(self,e):
        item = e.GetId()-10000-4*len(self.allrasd.list)
        self.allrasd.list[item].use_in_Fourier = self.datacontrol4[item].GetValue()

    def editUIR(self,e):
        item = e.GetId()-10000-5*len(self.allrasd.list)
        self.allrasd.list[item].use_in_Refine = self.datacontrol5[item].GetValue()
        
    def setdoabscorr(self,e): self.allrasd.do_abs_corr = self.getdoabscorr.GetValue()

    def setdwater(self,e):
        try:
            self.allrasd.d_water = float(e.GetString())
        except ValueError:
            pass

    def setconc(self,e):
        try:
            self.allrasd.conc = float(e.GetString())
        except ValueError:
            pass

    def setdkapton(self,e):
        try:
            self.allrasd.d_kapton = float(e.GetString())
        except ValueError:
            pass

    def OnClickabscorr(self,e):
        abs_correct(self.allrasd)
        for item in range(len(self.allrasd.list)):
            self.allrasd.list[item] = RASD_Fourier(self.allrasd, item)
            self.datacontrol1[item].SetValue(str(round(self.allrasd.list[item].AR,3)))
            self.datacontrol2[item].SetValue(str(round(self.allrasd.list[item].AR_err,3)))
            self.datacontrol3[item].SetValue(str(round(self.allrasd.list[item].PR,3)))
            self.datacontrol8[item].SetValue(str(round(self.allrasd.list[item].PR_err,3)))

    def setxf(self,e):
        try:
            self.allrasd.xf = float(e.GetString())
        except ValueError:
            pass
    def setyf(self,e):
        try:
            self.allrasd.yf = float(e.GetString())
        except ValueError:
            pass
    def setzf(self,e):
        try:
            self.allrasd.zf = float(e.GetString())
        except ValueError:
            pass
    def setzmin(self,e):
        try:
            self.allrasd.zmin = float(e.GetString())
        except ValueError:
            pass
    def setan(self,e):
        try:
            self.allrasd.an = int(e.GetString())
        except ValueError:
            pass
    def setbn(self,e):
        try:
            self.allrasd.bn = int(e.GetString())
        except ValueError:
            pass
    def setcn(self,e):
        try:
            self.allrasd.cn = int(e.GetString())
        except ValueError:
            pass

    def OnClickFSB(self,e):
        Fourier = Num.ndarray((0,5),float)
        for Rasd in self.allrasd.list:
            if Rasd.use_in_Fourier:
                F_comp = Num.array([[Rasd.Q[0],Rasd.Q[1],Rasd.Q[2],Rasd.AR,Rasd.PR]])
                Fourier = Num.append(Fourier,F_comp,axis=0)
        if len(Fourier) == 0:
            dlg = wx.MessageDialog(self, " No Fourier Components specified \n for use in Fourier Synthesis ","", wx.OK | wx.STAY_ON_TOP)
            if dlg.ShowModal() == wx.OK:
                dlg.Destroy()        
        else:
            if self.fframe != None:
                self.fframe.Close()
                self.fframe = None
            rho, Fcell = Fourier_synthesis(Fourier, self.allrasd.cell, self.allrasd.ZR, self.allrasd.xf, self.allrasd.yf, self.allrasd.zf,\
                          self.allrasd.an, self.allrasd.bn, self.allrasd.cn, self.allrasd.zmin)
            self.fframe = createfourierframe(self, rho, Fcell, self.allrasd.zmin*self.allrasd.cell[2])
            self.fframe.Show(True)
        
    def setRD(self,e): self.allrasd.Refine_Data = self.getRD.GetValue()
    def setulem(self, e): self.nb.MainControlPage.use_lay_el = self.ulem.GetValue()
    def setqmax(self,e):
        try:
            self.Qmax = float(e.GetString())
        except ValueError:
            pass
    def OnClickPlotArPrQ(self,e):
        res_param_usage = []
        res_surface = []
        for i in range(len(self.nb.surface)):
            if self.nb.surface[i][0] == self.resel:
                res_surface.append(self.nb.surface[i])
                res_param_usage.append(self.nb.parameter_usage[i])

        res_params = {}

        global_params = ['z_el','d_el','sig_el','sig_el_bar','K','occ_el_0','zwater','sig_water','sig_water_bar','d_water','beta','Scale','specScale']
        for i in range(len(res_param_usage)):
            for j in range(1,20,2):
                if res_param_usage[i][j] not in res_params and res_param_usage[i][j] != 'None':
                    res_params[res_param_usage[i][j]] = self.nb.parameter[res_param_usage[i][j]]
        for key in global_params:
            res_params[key] = self.nb.parameter[key]
			
        self.fig8 = plot_Ar_Pr_Q(self.fig8,  self.allrasd, self.Qmax, res_surface, res_params, res_param_usage, self.nb.MainControlPage.UBW_flag, self.nb.MainControlPage.use_lay_el)
        self.fig8.canvas.draw()
        
    def OnClickRefine(self,e):
        res_param_usage = []
        res_surface = []
        for i in range(len(self.nb.surface)):
            if self.nb.surface[i][0] == self.resel:
                res_surface.append(self.nb.surface[i])
                res_param_usage.append(self.nb.parameter_usage[i])

        res_params = {}

        global_params = ['z_el','d_el','sig_el','sig_el_bar','K','occ_el_0','zwater','sig_water','sig_water_bar','d_water','beta','Scale','specScale']
        for i in range(len(res_param_usage)):
            for j in range(1,20,2):
                if res_param_usage[i][j] not in res_params and res_param_usage[i][j] != 'None':
                    res_params[res_param_usage[i][j]] = self.nb.parameter[res_param_usage[i][j]]
        for key in global_params:
            res_params[key] = self.nb.parameter[key]
   
        self.allrasd, res_params = res_simplex(self.nb.frame.statusbar,res_params,res_param_usage, self.allrasd, res_surface, self.simplex_params,\
                                               self.nb.MainControlPage.UBW_flag, self.nb.MainControlPage.use_lay_el)

        if self.allrasd.Refine_Data:
            if self.fig7 == None:
                self.fig7 = figure(7)
            else:
                self.fig7.clf()
            self.fig7.suptitle('Overall chi**2: '+str(self.allrasd.RMS), fontsize = 20)
            rasdplot = self.fig7.add_subplot(111)
            rasdplot.set_xlabel('energy')
            rasdplot.set_ylabel('|F|**2')
            
            i = 0
            j = 0
            for Rasd in self.allrasd.list:
                if Rasd.use_in_Refine:
                    PR = Num.arctan(Rasd.im_Fq/Rasd.re_Fq)/(2*Num.pi)
                    AR = Rasd.re_Fq /Num.cos(2*Num.pi*PR)
                    if AR < 0 and PR < 0.5 and PR > 0.:
                        PR = PR + 0.5
                        AR = Rasd.re_Fq /Num.cos(2*Num.pi*PR)
                    elif AR < 0 and PR > 0.5:
                        PR = PR - 0.5
                        AR = Rasd.re_Fq /Num.cos(2*Num.pi*PR)
                    elif PR < 0 and AR > 0:
                        PR = PR + 1.
                    elif AR < 0 and PR < 0 and PR > -0.5:
                        PR = PR + 0.5
                        AR = Rasd.re_Fq /Num.cos(2*Num.pi*PR)
                    self.datacontrol6[j].SetValue(str(round(AR,3)))
                    self.datacontrol7[j].SetValue(str(round(PR,3)))
                    
                    Rasd.AR_refine = AR
                    Rasd.PR_refine = PR
                
                    rasdplot.errorbar(Rasd.E,(Rasd.F/Rasd.norm)+i*0.5, (Rasd.Ferr/Rasd.norm),fmt = 'bo')
                    rasdplot.plot(Rasd.E, (Rasd.F_calc/Rasd.norm)+i*0.5, 'g-' )
                    rasdplot.text(Rasd.E0,i*0.5+1,(Rasd.file +', chi**2 = ' +str(round(Rasd.RMS,4))))
                    i = i+1
                j = j+1
            self.fig7.canvas.draw()
        else:
            j = 0
            for Rasd in self.allrasd.list:
                if Rasd.use_in_Refine:
                    self.datacontrol6[j].SetValue(str(round(Rasd.AR_refine,3)))
                    self.datacontrol7[j].SetValue(str(round(Rasd.PR_refine,3)))
                j = j+1
        for key in res_params.keys():
            self.nb.parameter[key] = res_params[key]
        for i in range(len(self.nb.param_labels)):
            self.nb.ParameterPage.control1[i].SetValue(str(round(self.nb.parameter[self.nb.param_labels[i]][0], 12)))

    def OnClickResStatistics(self,e):
        self.nb.frame.SetStatusText(' computing parameter statistics ', 0)
        res_param_usage = []
        res_surface = []
        for i in range(len(self.nb.surface)):
            if self.nb.surface[i][0] == self.resel:
                res_surface.append(self.nb.surface[i])
                res_param_usage.append(self.nb.parameter_usage[i])

        res_params = {}

        global_params = ['z_el','d_el','sig_el','sig_el_bar','K','occ_el_0','zwater','sig_water','sig_water_bar','d_water','beta','Scale','specScale']
        for i in range(len(res_param_usage)):
            for j in range(1,20,2):
                if res_param_usage[i][j] not in res_params and res_param_usage[i][j] != 'None':
                    res_params[res_param_usage[i][j]] = self.nb.parameter[res_param_usage[i][j]]
        for key in global_params:
            res_params[key] = self.nb.parameter[key]
        res_params, correl_matrix, used_params = res_param_statistics(self.nb.MainControlPage.fpc,res_params,res_param_usage,\
                                                                                       self.allrasd, res_surface,self.nb.MainControlPage.UBW_flag, \
                                                                                       self.nb.MainControlPage.use_lay_el)
        
        for i in res_params.keys():
            if res_params[i][3]:
                self.nb.parameter[i][4] = res_params[i][4]
        for i in range(len(self.nb.param_labels)):    
            self.nb.ParameterPage.control8[i].SetValue(str(round(self.nb.parameter[self.nb.param_labels[i]][4], 8)))
            
        if correl_matrix[0][0] != 0:
            b = len(correl_matrix)
            print '\nParameter Correlations: \n'
            print '\nParameter Correlations between 0.95 and 1.00: \n'
            n = 0
            for i in range(b):
                for j in range(b):
                    if j > i:
                        if correl_matrix[i][j] >= 0.95 and correl_matrix[i][j] <= 1.0:
                            n = n+1
                            print used_params[i] + ' & ' + used_params[j] + ': ' + str(round(correl_matrix[i][j],5))
            if n == 0: print 'None \n'
            else: print '\n'
            n = 0
            print '\nParameter Correlations between 0.8 and 0.95: \n'
            for i in range(b):
                for j in range(b):
                    if j > i:
                        if correl_matrix[i][j] >= 0.8 and correl_matrix[i][j] < 0.95:
                            n = n+1
                            print used_params[i] + ' & ' + used_params[j] + ': ' + str(round(correl_matrix[i][j],5))
            if n == 0: print 'None \n'
            else: print '\n'
            n = 0
            print '\nParameter Correlations between 0.5 and 0.8: \n'
            for i in range(b):
                for j in range(b):
                    if j > i:
                        if correl_matrix[i][j] >= 0.5 and correl_matrix[i][j] < 0.5:
                            n = n+1
                            print used_params[i] + ' & ' + used_params[j] + ': ' + str(round(correl_matrix[i][j],5))
            if n == 0: print 'None \n'
            else: print '\n'
        self.nb.frame.SetStatusText(' statistics calculation finished, chi**2 = '+str(round(self.allrasd.RMS,3)), 0)
        pass
########################################################################################################

