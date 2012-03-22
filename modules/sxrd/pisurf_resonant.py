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
from scipy.optimize import leastsq
import wx
import os

from tdl.modules.sxrd.ctrfitcalcs import param_unfold,RB_update, calc_g_inv
from tdl.modules.sxrd.resonant_simplex import res_simplex
from tdl.modules.sxrd.fourierframe import createfourierframe

from tdl.modules.xtab.atomic import f0data as database
#################################################################################
class ResonantDataPanel(wx.ScrolledWindow):
    def __init__(self,parent):
        wx.ScrolledWindow.__init__(self,parent)
        
        self.nb = self.GetParent()
        self.resel = ''
        self.plot_norm = True

        self.allrasd = RasdList()
        self.allrasd.g_inv = self.nb.g_inv
        self.allrasd.cell = self.nb.cell

        self.fframe = None
        self.fig7 = None

        self.simplex_params = [1.0,0.5,2.0,1.0,1e-4,1e-6,10000]

        self.datacontrol0 = []
        self.datacontrol1 = []
        self.datacontrol2 = []
        self.datacontrol4 = []
        self.datacontrol5 = []
        self.datacontrol6 = []
        self.datacontrol7 = []

        wx.StaticText(self, label = ' Resonant Data ', pos=(290, 12), size=(120, 20))
        wx.StaticText(self, label = ' AR ', pos=(420, 12), size=(60, 20))
        wx.StaticText(self, label = ' PR ', pos=(490, 12), size=(60, 20))
        wx.StaticText(self, label = 'UIF', pos=(560, 12), size=(20, 20))
        wx.StaticText(self, label = 'UIR', pos=(590, 12), size=(20, 20))
        wx.StaticText(self, label = ' AR ', pos=(620, 12), size=(60, 20))
        wx.StaticText(self, label = ' PR ', pos=(690, 12), size=(60, 20))

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
        self.RefineButton = wx.Button(self, label = 'Run Refinement', pos =(20,485), size=(230,135))
        self.Bind(wx.EVT_BUTTON, self.OnClickRefine, self.RefineButton)
        
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
            self.datacontrol2[item].SetValue(str(round(self.allrasd.list[item].PR,3)))
            
    
    def ClickDataButton(self,e):
        item = e.GetId()-3000
        self.allrasd.list[item] = RASD_Fourier(self.allrasd, item)
        self.datacontrol1[item].SetValue(str(round(self.allrasd.list[item].AR,3)))
        self.datacontrol2[item].SetValue(str(round(self.allrasd.list[item].PR,3)))
        self.allrasd.list[item].plot(norm = self.plot_norm, fig = 4)

    def editPR(self,e):
        item = e.GetId()-3000-2*len(self.allrasd.list)
        try:
            self.allrasd.list[item].PR = float(e.GetString())
        except ValueError:
            pass

    def editUIF(self,e):
        item = e.GetId()-3000-4*len(self.allrasd.list)
        self.allrasd.list[item].use_in_Fourier = self.datacontrol4[item].GetValue()

    def editUIR(self,e):
        item = e.GetId()-3000-5*len(self.allrasd.list)
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
            self.datacontrol1[item].SetValue(str(round(self.allrasd.list[item].AR,3)))
            self.datacontrol2[item].SetValue(str(round(self.allrasd.list[item].PR,3)))

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
   
        self.allrasd, res_params = res_simplex(res_params,res_param_usage, self.allrasd, res_surface, self.simplex_params,\
                                               self.nb.MainControlPage.UBW_flag, self.allrasd.Refine_Data, self.nb.MainControlPage.use_lay_el)

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
            

        dlg = wx.MessageDialog(self, "Keep refined parameters ?","", wx.YES_NO | wx.STAY_ON_TOP)
        if dlg.ShowModal() == wx.ID_YES:
            for key in res_params.keys():
                self.nb.parameter[key] = res_params[key]
            for i in range(len(self.nb.param_labels)):
                self.nb.ParameterPage.control1[i].SetValue(str(round(self.nb.parameter[self.nb.param_labels[i]][0], 12)))
        dlg.Destroy()
########################################################################################################
###################################### calculations ####################################################
def calc_FucNR(hkl,bulk,g_inv,database):
    a = 0
    b = 0
    q = (Num.dot(Num.dot(hkl,g_inv),hkl))**0.5
    for i in range(shape(bulk)[0]):
        f_par = database[str.lower(bulk[i][0])]
        f = (f_par[0]*Num.exp(-(q/4/Num.pi)**2*f_par[1]) + f_par[2]*Num.exp(-(q/4/Num.pi)**2*f_par[3]) +\
            f_par[4]*Num.exp(-(q/4/Num.pi)**2*f_par[5]) + f_par[6]*Num.exp(-(q/4/Num.pi)**2*f_par[7]) + f_par[8])*\
            Num.exp(-2 * Num.pi**2 * q**2 * bulk[i][4])
        a = a + (f * Num.cos(2*Num.pi*(hkl[0]*bulk[i][1] + hkl[1]*bulk[i][2] + hkl[2]*bulk[i][3])))
        b = b + (f * Num.sin(2*Num.pi*(hkl[0]*bulk[i][1] + hkl[1]*bulk[i][2] + hkl[2]*bulk[i][3])))
    return a, b

def calc_FsurfNR(hkl,surface,g_inv,database):
    a = 0
    b = 0
    q = (Num.dot(Num.dot(hkl,g_inv),hkl))**0.5
    q_Ang = [hkl[0]*g_inv[0][0]**0.5, hkl[1]*g_inv[1][1]**0.5, hkl[2]*g_inv[2][2]**0.5]
    for i in range(shape(surface)[0]):
        f_par = database[str.lower(surface[i][0])]

        U = Num.ndarray((3,3),float)
        U[0][0] = surface[i][4]
        U[0][1] = surface[i][7]*(surface[i][4])**0.5*(surface[i][5])**0.5
        U[0][2] = surface[i][8]*(surface[i][4])**0.5*(surface[i][6])**0.5
        U[1][0] = U[0][1]
        U[1][1] = surface[i][5]
        U[1][2] = surface[i][9]*(surface[i][5])**0.5*(surface[i][6])**0.5
        U[2][0] = U[0][2]
        U[2][1] = U[1][2]
        U[2][2] = surface[i][6]
        
        f = (f_par[0]*Num.exp(-(q/4/Num.pi)**2*f_par[1]) + f_par[2]*Num.exp(-(q/4/Num.pi)**2*f_par[3]) +\
            f_par[4]*Num.exp(-(q/4/Num.pi)**2*f_par[5]) + f_par[6]*Num.exp(-(q/4/Num.pi)**2*f_par[7]) + f_par[8])*\
            Num.exp(-2* Num.pi**2*(Num.dot(q_Ang,Num.dot(U,q_Ang)))) * surface[i][10]
        a = a + (f * Num.cos(2*Num.pi*(hkl[0]*surface[i][1] + hkl[1]*surface[i][2] + hkl[2]*surface[i][3])))
        b = b + (f * Num.sin(2*Num.pi*(hkl[0]*surface[i][1] + hkl[1]*surface[i][2] + hkl[2]*surface[i][3])))
    return a, b

def calc_Fwater_layeredNR(hkl, sig, sig_bar, d,zwater, g_inv, database, cell):
    f_par = database['o2-.']
    q = hkl[2]* g_inv[2][2]**0.5
    Auc = cell[0]* Num.sin(Num.radians(cell[5]))* cell[1]
    f = Auc * d * 0.033456 * (f_par[0]*Num.exp(-(q/4/Num.pi)**2*f_par[1]) + f_par[2]*Num.exp(-(q/4/Num.pi)**2*f_par[3]) +\
            f_par[4]*Num.exp(-(q/4/Num.pi)**2*f_par[5]) + f_par[6]*Num.exp(-(q/4/Num.pi)**2*f_par[7]) + f_par[8])*\
            Num.exp(-2 * Num.pi**2 * q**2 * sig)
    x = Num.pi * q * d
    al = 2 * Num.pi**2 * q**2 * sig_bar
    a = Num.exp(al)*Num.cos(2*x)-1
    b = Num.exp(al)*Num.sin(-2*x)
    c = 4 * Num.cos(x)**2 * Num.sinh(al/2)**2 - 4 * Num.sin(x)**2 * Num.cosh(al/2)**2
    d = -2 * Num.sin(2*x) * Num.sinh(al)
    rez = Num.cos(2*Num.pi*q*zwater/g_inv[2][2]**0.5)
    imz = Num.sin(2*Num.pi*q*zwater/g_inv[2][2]**0.5)
    relayer = (a*c + b*d)/(c**2 + d**2)
    imlayer = (b*c - a*d)/(c**2 + d**2)
    re = f* (relayer * rez - imlayer * imz)
    im = f* (relayer * imz + imlayer * rez)
    return re, im
        
def calc_F_layered_el_NR(hkl, occ, K, sig, sig_bar, d, d0, g_inv, database, el):
    f_par = database[el]
    zinv = g_inv[2][2]**0.5
    q = hkl[2]* zinv
    qd4pi = q/4/Num.pi
    f = (f_par[0]*Num.exp(-(qd4pi)**2*f_par[1]) + f_par[2]*Num.exp(-(qd4pi)**2*f_par[3]) +\
            f_par[4]*Num.exp(-(qd4pi)**2*f_par[5]) + f_par[6]*Num.exp(-(qd4pi)**2*f_par[7]) + f_par[8])*\
            Num.exp(-2 * Num.pi**2 * q**2 * sig)*occ
    x = Num.pi * q * d
    al = 2 * Num.pi**2 * q**2 * sig_bar + K * d
    a = Num.exp(al)*Num.cos(2*x)-1
    b = Num.exp(al)*Num.sin(-2*x)
    c = 4 * Num.cos(x)**2 * Num.sinh(al/2)**2 - 4 * Num.sin(x)**2 * Num.cosh(al/2)**2
    d = -2 * Num.sin(2*x) * Num.sinh(al)
    wert = 2*Num.pi*hkl[2]*d0
    rez = Num.cos(wert)
    imz = Num.sin(wert)
    wert = c**2 + d**2
    relayer = (a*c + b*d)/(wert)
    imlayer = (b*c - a*d)/(wert)
    re = f* (relayer * rez - imlayer * imz)
    im = f* (relayer * imz + imlayer * rez)
    return re, im

def calcFNR(hkl, global_parms, cell, bulk, surface, database, g_inv, use_bulk_water, use_lay_el, el):
    
    occ_el, K,sig_el,sig_el_bar,d_el,d0_el,sig_water, sig_water_bar, d_water, zwater, Scale, specScale, beta= global_parms

    zeta = hkl[2]+ hkl[0]*cell[6]+ hkl[1]*cell[7]

    re_ctr = 0.5
    im_ctr = -1/(2*Num.tan(Num.pi*zeta))

    re_bulk , im_bulk = calc_FucNR(hkl,bulk,g_inv,database)
    re_surf, im_surf = calc_FsurfNR(hkl,surface,g_inv,database)
            
    re_bc = re_ctr*re_bulk - im_ctr*im_bulk
    im_bc = re_bulk*im_ctr + re_ctr*im_bulk
    
    if hkl[0] == 0.0 and hkl[1] == 0.0:
        if not use_bulk_water:
            re_water = 0
            im_water = 0
        else:
            re_water, im_water = calc_Fwater_layeredNR(hkl, sig_water, sig_water_bar, d_water,zwater, g_inv, database, cell)
        if not use_lay_el:
            re_el = 0
            im_el = 0
        else:
            re_el, im_el = calc_F_layered_el_NR(hkl, occ_el, K, sig_el, sig_el_bar, d_el, d0_el, g_inv, database, el)
                    
        re_FNR = (re_bc + re_surf + re_water + re_el)
        im_FNR = (im_bc + im_surf + im_water + im_el)
    else:
        re_FNR = (re_bc + re_surf)
        im_FNR = (im_bc + im_surf)
    return re_FNR, im_FNR
##########################################################################################
##################################################################################        
class RasdList:
    """
    RasdList is a container that holds RasdAna objects and Parameters to analyze them
    """
    def __init__(self):
        #List of RasdAna Objects
        self.list = []
        #General stuff
        self.dims = 0
        self.E0 = 0.
        self.E = Num.array([], float)
        self.f1 = Num.array([], float)
        self.f2 = Num.array([], float)
        self.cell = []
        self.g_inv = Num.array([[]],float)
        self.Refine_Data = True
        #Parameters for absorption correction
        self.do_abs_corr = False
        self.d_water = 0.
        self.conc = 0.
        self.d_kapton = 0.
        #Parameters for Fourier synthesis
        self.Fourier = []
        self.ZR = 1
        self.xf = 1
        self.yf = 1
        self.zf = 1
        self.zmin = 0
        self.an = 11
        self.bn = 11
        self.cn = 61
        
#####################################################################################    
class RasdAna:
    """
    Rasd Data object without image data to be used in Modell Refinement and Fourier component extraction 
    """
    def __init__(self):
        self.E = Num.array([], float) #Energies to be changed by e0shift
        self.Eorig = Num.array([], float) #Energies to remain unchanged
        self.F = Num.array([], float) # #this is (in case) the absorption corrected structure factor
        self.Ferr = Num.array([], float) 
        self.f1 = Num.array([], float) #f1 to be shifted by e0shift
        self.f2 = Num.array([], float) #f2 to be shifted by e0shift
        self.file = '' #*.rsd file

        self.Q = Num.ndarray((3), float) #q in fractional reciprocal coordinates
        self.mod_Q = float #norm Q in Ang**-1
        self.re_FNR = float #real part of non resonant Structure Factor
        self.im_FNR = float #imaginary part of non resonant Structure Factor

        self.re_FR = Num.array([], float) #real part of resonant Structure Factor
        self.im_FR = Num.array([], float) #imaginary part of resonant Structure Factor
        self.F_calc = Num.array([], float) #modulus of calculated structure factor
        self.a = 0.1
        self.b = 0.1
        self.E0 = float
        self.e0shift = 0.
        self.Emin = float
        self.Emax = float
        self.norm = Num.array([], float)
        self.delta_F = Num.array([], float)
        self.AR = 0.2
        self.PR = 0.5
        self.AR_refine = 0.
        self.PR_refine = 0.
        self.re_Fq = 0.
        self.im_Fq = 0.
        self.use_in_Fourier = False
        self.use_in_Refine = False
        self.RMS = 0.
        self.ndata = 0
        #Parameters to calculate absorption correction
        self.abs_corr = Num.array([], float) # absorption correction factor
        self.Alpha = Num.array([], float)
        self.Beta = Num.array([], float)

    ########################################################################
    def plot(self, norm = True, fig = 1):
        if self.F_calc == []:
            figure(fig)
            clf()
            title(str(self.file))
            errorbar(self.E, self.F, self.Ferr, fmt = 'o')
        else:
            if norm:
                figure(fig)
                clf()
                title(str(self.file)+'  chi**2: '+str(self.RMS))
                errorbar(self.E,self.F/self.norm, self.Ferr/self.norm,fmt = 'o')
                plot(self.E, self.F_calc/self.norm)
            else:
                figure(fig)
                clf()
                title(str(self.file)+'  chi**2: '+str(self.RMS))
                errorbar(self.E,self.F, self.Ferr,fmt = 'o')
                plot(self.E, self.F_calc)

#################################################################################
def read_f1f2(allrasd, f1f2_file):
    f = file(f1f2_file, 'r')
    data = f.readlines()
    f.close()
    f1f2 = Num.ndarray((0,3),float)
    for i in range(len(data)):
        if '#' not in data[i]:
            tmp = str.rsplit(data[i])
            f1f2 = Num.append(f1f2, [[int(float(tmp[0])),float(tmp[1]),float(tmp[2])]], axis = 0)
    f1f2 = f1f2.transpose()
    allrasd.E = f1f2[0]
    allrasd.f1 = f1f2[1]
    allrasd.f2 = f1f2[2]
    return allrasd
#############################################################################################################################
def read_RSD(allrasd, bulk_tmp, surface_tmp, parameter, param_usage, rigid_bodies, database, rasddata, use_bulk_water, use_lay_el, el, dirname, parent):
    """    
    read in data in RasdList of RasdAna objects for fitting
   """
    if len(allrasd.E)== 0 and len(allrasd.f1) == 0:
        dlg = wx.MessageDialog(parent, "Please read in f1f2 data first","", wx.OK | wx.STAY_ON_TOP)
        if dlg.ShowModal() == wx.OK:
            dlg.Destroy()
        return allrasd
    elif allrasd.E0 == 0:
        dlg = wx.MessageDialog(parent, "Please specify Edge Energy first","", wx.OK | wx.STAY_ON_TOP)
        if dlg.ShowModal() == wx.OK:
            dlg.Destroy()
        return allrasd
    else:
        global_parms, surface_new = param_unfold(parameter,param_usage, surface_tmp, use_bulk_water, use_lay_el)
        bulk = bulk_tmp 
        surface = RB_update(rigid_bodies, surface_new, parameter, allrasd.cell)

        cell = allrasd.cell
        g_inv = allrasd.g_inv
        allrasd.ndata = 0
        allrasd.RMS = 0
        
        for dat in rasddata:
            n = len(dat)
            if dat[n-1] == '\n':
                dat = str.rstrip(dat,'\n')
            Q = None
            f = file(dirname+'/'+dat, 'r')
            data = f.readlines()
            f.close()
            Rasd = RasdAna()
            for i in range(len(data)):
                if '#' not in data[i]:
                    tmp = str.rsplit(data[i])
                    if Q == None:
                        Q = Num.array([tmp[1],tmp[2],tmp[3]], float)
                        Rasd.Q = Q
                        Rasd.mod_Q = Num.dot(Num.dot(Q,g_inv),Q)**0.5
                        Rasd.re_FNR ,Rasd.im_FNR = calcFNR(Q,global_parms, cell, bulk, surface, database, g_inv, use_bulk_water, use_lay_el, el)
                    
                    Rasd.E = Num.append(Rasd.E, int(round(float(tmp[0]))))
                    Rasd.F = Num.append(Rasd.F, float(tmp[4])**2)
                    Rasd.Ferr = Num.append(Rasd.Ferr, float(tmp[5])**2)
                    Rasd.Alpha = Num.append(Rasd.Alpha, float(tmp[6]))
                    Rasd.Beta = Num.append(Rasd.Beta, float(tmp[7]))
            Rasd.E0 = allrasd.E0
            for x in Rasd.E:
                Rasd.Eorig = Num.append(Rasd.Eorig, x)
            Rasd.abs_corr = Num.ones((len(Rasd.E)),float)
            Rasd.ndata = len(Rasd.E)
            Rasd.file = dat
            allrasd.ndata = allrasd.ndata + Rasd.ndata
        
            i=0
            Rasd.f1 = Num.ndarray((Rasd.ndata), float)
            Rasd.f2 = Num.ndarray((Rasd.ndata), float)
            for i in range(Rasd.ndata):
                j = 0
                for j in range(len(allrasd.E)):
                    if Rasd.E[i] == allrasd.E[j]:
                        Rasd.f1[i] = allrasd.f1[j]
                        Rasd.f2[i] = allrasd.f2[j]
        
            allrasd.list.append(Rasd)
        allrasd.dims = len(allrasd.list)
        allrasd.RMS = allrasd.RMS / allrasd.ndata

        return allrasd
####################################################################################################
def write_rids(allrasd, filename):
    f = file(filename, 'w')
    for Rasd in allrasd.list:
        if rasd.use_in_refine:
            f.write('#\n RIDS scan at: '+str(round(Rasd.Q[0],1))+' '+str(round(Rasd.Q[1],1))+' '+str(round(Rasd.Q[2],3))+'\n')
            f.write('AR = '+str(round(Rasd.AR,3))+' PR = '+str(round(Rasd.PR,3))+'\n')
            f.write('AR_model = '+str(round(Rasd.AR_refine,3))+' PR_model = '+str(round(Rasd.PR_refine,3))+'\n')
            f.write('chi**2 = '+str(round(Rasd.RMS,5))+'\n\n')
            f.write('# E    F**2           F_err**2        F_model**2 \n')
            for i in range(len(Rasd.E)):
                line = "%7.2f %15.5f %15.5f %15.5f\n" % (Rasd.E[i], Rasd.F[i]/Rasd.norm[i], Rasd.Ferr[i]/Rasd.norm[i], Rasd.F_calc[i]/Rasd.norm[i])
                f.write(line)
    f.close()
########################## absorption correction #####################################################
def abs_correct(allrasd):
    if allrasd.do_abs_corr:
        for rasd in allrasd.list:
            rasd.F = rasd.F * rasd.abs_corr
            rasd.abs_corr = Num.ones((len(rasd.F)),float)
            for i in range(len(rasd.F)):
                z_sol = allrasd.d_water/ Num.sin(Num.radians(rasd.Alpha[i])) + allrasd.d_water/ Num.sin(Num.radians(rasd.Beta[i]))
                mu_sol = (allrasd.conc * rasd.f2[i]  + 55555 * 8.2e-3)* 6.022e23 * 2 * 2.82e-15 *12398e-10 /rasd.E[i]
                z_kapt = allrasd.d_kapton/ Num.sin(Num.radians(rasd.Alpha[i])) + allrasd.d_kapton/ Num.sin(Num.radians(rasd.Beta[i]))
                mu_kapt = 1/7220e-6
                rasd.abs_corr[i] = Num.exp(-z_sol* mu_sol -z_kapt * mu_kapt)
            rasd.F = rasd.F / rasd.abs_corr
    else:
        for rasd in allrasd.list:
            rasd.F = rasd.F * rasd.abs_corr
            rasd.abs_corr = Num.ones((len(rasd.F)),float)
            
    for i in range(len(allrasd.list)):
        allrasd.list[i] = RASD_Fourier(allrasd, i)            
########################## Fit Fourier Components ####################################################        
def RASD_Fourier(allrasd, pnt):
    Rasd = allrasd.list[pnt]
    for i in range(Rasd.ndata):
        j = 0
        for j in range(len(allrasd.E)):
            if Rasd.E[i] == allrasd.E[j]:
                Rasd.f1[i] = allrasd.f1[j]
                Rasd.f2[i] = allrasd.f2[j]

    vec = Num.array([Rasd.a ,Rasd.b , Rasd.AR, Rasd.PR], float)

    def Wiggle(vec, Rasd):
        a , b, AR, PR = vec
        re_Fq = AR * Num.cos(2*Num.pi*PR)
        im_Fq = AR * Num.sin(2*Num.pi*PR)
        re_FR = Rasd.f1 * re_Fq - Rasd.f2 * im_Fq
        im_FR = Rasd.f1 * im_Fq + Rasd.f2 * re_Fq
        return Rasd.F - ((a + b*(Rasd.E-Rasd.E0))* ((Rasd.re_FNR + re_FR)**2 + (Rasd.im_FNR + im_FR)**2))

    def Jacobi(vec, Rasd):
        J = Num.ndarray((4, Rasd.ndata), float)
        a , b, AR, PR = vec
        re_Fq = AR * Num.cos(2*Num.pi*PR)
        im_Fq = AR * Num.sin(2*Num.pi*PR)
        re_FR = Rasd.f1 * re_Fq - Rasd.f2 * im_Fq
        im_FR = Rasd.f1 * im_Fq + Rasd.f2 * re_Fq

        J[0,:] = -(Rasd.re_FNR + re_FR)**2 - (Rasd.im_FNR + im_FR)**2

        J[1,:] = -(Rasd.E-Rasd.E0) * ((Rasd.re_FNR + re_FR)**2 + (Rasd.im_FNR + im_FR)**2)
        J[2,:] = -2 *(a+b*(Rasd.E-Rasd.E0))* ((Rasd.f1*Num.cos(2*Num.pi*PR)-Rasd.f2*Num.sin(2*Num.pi*PR))*(Rasd.re_FNR + re_FR)+\
                                          (Rasd.f1*Num.sin(2*Num.pi*PR)+Rasd.f2*Num.cos(2*Num.pi*PR))*(Rasd.im_FNR + im_FR))
                 
        J[3,:] = -4 *Num.pi* (a+b*(Rasd.E-Rasd.E0))* ((-Num.sin(2*Num.pi*PR)*AR*Rasd.f1- Num.cos(2*Num.pi*PR)*Rasd.f2*AR)*(Rasd.re_FNR + re_FR) +\
                                                  ( Num.cos(2*Num.pi*PR)*AR*Rasd.f1- Num.sin(2*Num.pi*PR)*Rasd.f2*AR)*(Rasd.im_FNR + im_FR))
        return J
    
    result = leastsq(Wiggle, vec,args = (Rasd), Dfun = Jacobi, col_deriv = 1)
    Rasd.a , Rasd.b, Rasd.AR, Rasd.PR = result[0]
                          
    re_Fq = Rasd.AR * Num.cos(2*Num.pi*Rasd.PR)
    im_Fq = Rasd.AR * Num.sin(2*Num.pi*Rasd.PR)
    Rasd.re_FR = Rasd.f1 * re_Fq - Rasd.f2 * im_Fq
    Rasd.im_FR = Rasd.f1 * im_Fq + Rasd.f2 * re_Fq

    Rasd.F_calc = (Rasd.a + Rasd.b*(Rasd.E-Rasd.E0))* ((Rasd.re_FNR + Rasd.re_FR)**2 + (Rasd.im_FNR + Rasd.im_FR)**2)   
    Rasd.delta_F = (Rasd.F - Rasd.F_calc)**2/Rasd.Ferr**2
    Rasd.norm = (Rasd.a+Rasd.b*(Rasd.E-Rasd.E0))*(Rasd.re_FNR**2+Rasd.im_FNR**2)
    Rasd.RMS = Num.sum(Rasd.delta_F)/Rasd.ndata
    
    return Rasd
#################  Fourier Synthese  ###########################################
def calc_rho(Fourier, r, ginv):
    rho = 0
    for F_comp in Fourier:
        Q = Num.array([F_comp[0]*ginv[0][0]**0.5,F_comp[1]*ginv[1][1]**0.5,\
                      F_comp[2]*ginv[2][2]**0.5], float)
        rho = rho + (F_comp[3] * Num.cos(2*Num.pi*(F_comp[4] - Num.dot(r,Q))))
    return rho

def Fourier_synthesis(Fourier, cell, ZR, xf, yf, zf, an, bn, cn, zmin):
    
    Rho = Num.ndarray((an,bn,cn),float)
    sampx = cell[0]* xf
    sampy = cell[1]* yf
    sampz = cell[2]* zf
    g_inv = calc_g_inv(cell)
    g_inv2 = calc_g_inv([sampx,sampy,sampz,cell[3],cell[4],cell[5]])
    V = Num.linalg.det(Num.linalg.inv(g_inv2))**0.5
    
    for i in range(int(an)):
        x = sampx /an * float(i)
        for j in range(int(bn)):
            y = sampy /bn * float(j)
            for k in range(int(cn)):
                z = sampz /cn * float(k) + sampz/zf *zmin
                R = Num.array([x,y,z],float)
                Rho[i][j][k]= calc_rho(Fourier, R, g_inv)
    Rho = Rho * ZR/(V*2*Num.pi)
    return Rho, [sampx,sampy,sampz]
