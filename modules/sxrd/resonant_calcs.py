"""
Functions used for resonant data analysis in pi_surf

Authors/Modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""
##############################################################################
import numpy as Num
from pylab import *
from scipy.optimize import leastsq
from scipy.interpolate import interp1d

from tdl.modules.sxrd.ctrfitcalcs import param_unfold,RB_update, calc_g_inv
from tdl.modules.sxrd.resonant_simplex import calc_F_lay_el
from tdl.modules.xtab.atomic import f0data as database
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
    rez = Num.cos(2*Num.pi*hkl[2]*zwater)
    imz = Num.sin(2*Num.pi*hkl[2]*zwater)
    relayer = (a*c + b*d)/(c**2 + d**2)
    imlayer = (b*c - a*d)/(c**2 + d**2)
    re = f* (relayer * rez - imlayer * imz)
    im = f* (relayer * imz + imlayer * rez)
    return re, im
        
def calc_F_layered_el_NR(hkl, occ, K, sig, sig_bar, d, d0, g_inv, database, el):
    f_par = database[el]
    q = hkl[2]* g_inv[2][2]**0.5
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
            f1f2 = Num.append(f1f2, [[float(tmp[0]),float(tmp[1]),float(tmp[2])]], axis = 0)
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
        f1 = interp1d(allrasd.E, allrasd.f1, 'cubic')
        f2 = interp1d(allrasd.E, allrasd.f2, 'cubic')
        
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
                    
                    Rasd.E = Num.append(Rasd.E, float(tmp[0]))
                    Rasd.F = Num.append(Rasd.F, float(tmp[4])**2)
                    Rasd.f1 = Num.append(Rasd.f1, f1(Rasd.E[-1]))
                    Rasd.f2 = Num.append(Rasd.f2, f2(Rasd.E[-1]))
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
            allrasd.list.append(Rasd)
        allrasd.dims = len(allrasd.list)
        allrasd.RMS = allrasd.RMS / allrasd.ndata

        return allrasd
####################################################################################################
def write_rids(allrasd, filename):
    f = file(filename, 'w')
    for Rasd in allrasd.list:
        if Rasd.use_in_Refine:
            f.write('#\n#RIDS scan at: '+str(round(Rasd.Q[0],1))+' '+str(round(Rasd.Q[1],1))+' '+str(round(Rasd.Q[2],3))+'\n')
            f.write('#AR = '+str(round(Rasd.AR,3))+' PR = '+str(round(Rasd.PR,3))+'\n')
            f.write('#AR_model = '+str(round(Rasd.AR_refine,3))+' PR_model = '+str(round(Rasd.PR_refine,3))+'\n')
            f.write('#chi**2 = '+str(round(Rasd.RMS,5))+'\n#\n')
            f.write('#E               F**2            F_err**2        F_model**2 \n')
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
	
############################################################################################################################################	
def plot_Ar_Pr_Q(allrasd, Qmax, surface_tmp, parameter, param_usage, use_bulk_water, use_lay_el):
    global_parms, surface = param_unfold(parameter,param_usage, surface_tmp, use_bulk_water, use_lay_el)
    natoms = len(surface)
    occ_el, K,sig_el,sig_el_bar,d_el,d0_el,sig_water, sig_water_bar, d_water, zwater, Scale, specScale, beta= global_parms
    A = []
    P = []
    for l in Num.arange(0,Qmax,0.01):
	Q = Num.array([0,0,l],float)
	Q_ang = Num.array([0,0,l*allrasd.g_inv[2][2]**0.5],float)
	re_Fq = 0
	im_Fq = 0
        for n in range(natoms):
            R = Num.array([surface[n][1],surface[n][2],surface[n][3]],float)

            U = Num.ndarray((3,3),float)
            U[0][0] = surface[n][4]
            U[0][1] = U[1][0] = surface[n][7]*(surface[n][4]*surface[n][5])**0.5
            U[0][2] = U[2][0] = surface[n][8]*(surface[n][4]*surface[n][6])**0.5
            U[1][1] = surface[n][5]
            U[1][2] = U[2][1] = surface[n][9]*(surface[n][5]*surface[n][6])**0.5
            U[2][2] = surface[n][6]

            theta = surface[n][10]
            re_Fq = re_Fq + (theta * Num.cos(2*Num.pi*(Q[0]*R[0]+Q[1]*R[1]+Q[2]*R[2])) * Num.exp(-2* Num.pi**2 * Num.dot(Num.dot(Q_ang,U),Q_ang)))
            im_Fq = im_Fq + (theta * Num.sin(2*Num.pi*(Q[0]*R[0]+Q[1]*R[1]+Q[2]*R[2])) * Num.exp(-2* Num.pi**2 * Num.dot(Num.dot(Q_ang,U),Q_ang)))
	if use_lay_el:
	    re_lay, im_lay = calc_F_lay_el(Q, occ_el, K, sig_el, sig_el_bar, d_el, d0_el, allrasd.g_inv)
            re_Fq = re_Fq + re_lay
            im_Fq = im_Fq + im_lay                                    

        PR = Num.arctan(im_Fq/re_Fq)/(2*Num.pi)
        AR = re_Fq /Num.cos(2*Num.pi*PR)
		
	if AR < 0 and PR < 0.5 and PR > 0.:
            PR = PR + 0.5
            AR = -AR
        elif AR < 0 and PR > 0.5:
            PR = PR - 0.5
            AR = -AR
        elif PR < 0 and AR > 0:
            PR = PR + 1.
        elif AR < 0 and PR < 0 and PR > -0.5:
            PR = PR + 0.5
            AR = -AR
	if PR > 1: 
            PR = PR -1			
        A.append(AR)
        P.append(PR)
    q_data = []
    A_data = []
    P_data = []
    for rasd in allrasd.list:
        if rasd.Q[0] == 0 and rasd.Q[1] == 0:
            q_data.append(rasd.Q[2])
            A_data.append(rasd.AR)
            P_data.append(rasd.PR)
            
    figure(8)
    clf()
    plot(Num.arange(0,Qmax,0.01),A,'b-')
    plot(Num.arange(0,Qmax,0.01),P,'r-')
    plot(q_data,A_data,'bo')
    plot(q_data,P_data,'ro')
		
    return
############################################ f1f2 transformation - Differential Kramers Kroning ####################################
def f1f2(datafile, expfile, e0, e0shift, output ='exp.f1f2', n=30):

    """
    Function to calculate Differential Kramers Kronig transformation from
    experimental f2 (normalized XANES) to experimental f1

    Literature: Ohta and Ishida (1988) ... integration Methods for Kramers Kronig Transformation, Applied Spectroscopy 42,6
                Cross et al. (1998) ... theoretical x-ray resonant scattering amplitudes ...,Physical Review B, 58, 17

    datafile: Hephaestus *.f1f2 file (1eV steps) --> Cromer Liberman (CL) f1 and f2
    expfile: Athena normalized XANES file (same eV grid as .f1f2 file) --> experimental f2
    output: filename experimental E f1 f2 wll be written to (for use in rasd_menu)
    n = 30: number of datapoints to match f2 CL with f2 exp.
    e0: theoretical edge energy
    e0shift: to match theoretical edge with experimental edge
    """
    

    #f1f2 holds
    #[0] - E exp
    #[1] - f1 theo (->E exp)
    #[2] - f2 theo (->E exp)
    #[3] - diff f1theo f1exp 
    #[4] - f1 exp
    #[5] - f2 exp
    
    e0 = e0 + e0shift
    #read Cromer-Liberman calculated f1f2 from HEPHAESTUS
    f = file(datafile, 'r')
    data = f.readlines()
    f.close()
    theo = Num.ndarray((0,3),float)
    for i in range(len(data)):
        if '#' not in data[i]:
            tmp = str.rsplit(data[i])
            if float(tmp[0])< e0:
                theo = Num.append(theo, [[int(round(float(tmp[0]))),float(tmp[1]),float(tmp[2])]], axis = 0)
            else:
                theo = Num.append(theo, [[int(round(float(tmp[0]))),float(tmp[1]),float(tmp[2])]], axis = 0)
    theo = theo.transpose()
    theo[0] = theo[0] + e0shift
    f1 = interp1d(theo[0], theo[1], 'cubic')
    f2 = interp1d(theo[0], theo[2], 'cubic')

    #read experimental f2 (normalized XANES) and interpolate theoretical f1 and f2 to match the E of f2 exp
    f = file(expfile, 'r')
    data = f.readlines()
    f.close()
    f1f2 = Num.ndarray((0,6),float)
    for i in range(len(data)):
        if '#' not in data[i]:
            tmp = str.rsplit(data[i])
            f1f2 = Num.append(f1f2, [[float(tmp[0]),f1(float(tmp[0])), f2(float(tmp[0])), 0., 0., float(tmp[1])]], axis = 0)
    f1f2 = f1f2.transpose()

    low_cl = 0
    up_cl = 0
    low_exp = 0
    up_exp = 0
    for i in range(n):
        low_cl = low_cl + f1f2[2][i]
        low_exp = low_exp + f1f2[5][i]
        up_cl = up_cl + f1f2[2][len(f1f2[0])-i-1]
        up_exp = up_exp + f1f2[5][len(f1f2[0])-i-1]
    
    low_cl= low_cl /n
    low_exp= low_exp /n
    up_cl = up_cl /n
    up_exp = up_exp /n
    rel = (up_cl-low_cl)/(up_exp-low_exp)
    diff = low_cl - low_exp*rel
    f1f2[5] = f1f2[5]*rel + diff

    #calculate f1exp from f2exp (Kramers-Kronig) (see Ohta 1988/ Cross 1998)
    for i in range(Num.shape(f1f2)[1]):
        sum = 0
        if divmod(float(i),2)[1] == 0:
            for j in range(1, len(f1f2[0]),2):
                sum = sum + (f1f2[5][j]-f1f2[2][j])/(f1f2[0][j] - f1f2[0][i])+(f1f2[5][j]-f1f2[2][j])/(f1f2[0][j] + f1f2[0][i])
        else:
            for j in range(0, len(f1f2[0]),2):
                sum = sum + (f1f2[5][j]-f1f2[2][j])/(f1f2[0][j] - f1f2[0][i])+(f1f2[5][j]-f1f2[2][j])/(f1f2[0][j] + f1f2[0][i])
        f1f2[3][i] = (sum * 2 / Num.pi)


    f1f2[4] = f1f2[1] + f1f2[3]

    #write experimental values to .f1f2 file
    f = file(output,'w')
    f.write('# file: "'+output+'" containing experimental f1 f2 values \n')
    f.write('# calculated using Cromer Liberman f1 f2 from: "'+datafile+'"\n')
    f.write('# and experimental f2 from: "'+expfile+'"\n')
    f.write('# E0 = '+str(e0)+', e0shift = '+str(e0shift)+'\n')
    f.write('# Energy f1exp f2exp \n')
    i=0
    for i in range(len(f1f2[0])):
        f.write(str(f1f2[0][i])+'   '+str(f1f2[4][i])+'   '+str(f1f2[5][i])+' \n')
    f.close()

    #plot results
    figure(1)
    clf()
    plot(f1f2[0],f1f2[1],'b-')
    plot(f1f2[0],f1f2[2],'b-')
    plot(f1f2[0],f1f2[4],'r-')
    plot(f1f2[0],f1f2[5],'g-')
