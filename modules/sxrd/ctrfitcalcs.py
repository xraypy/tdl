"""
Ctr fit calculations used in pi-surf

Authors/modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""
####################################################

import numpy as Num
from pylab import *
import random
import wx

from tdl.modules.xtal.bv_params import bv_params
######################  general calculations  ##################################
def calc_g_inv(cell):
    g = Num.ndarray((3,3),float)
    g[0][0] = cell[0]**2
    g[0][1] = g[1][0] = cell[0]*cell[1]*Num.cos(Num.radians(cell[5]))
    g[0][2] = g[2][0] = cell[0]*cell[2]*Num.cos(Num.radians(cell[4]))
    g[1][1] = cell[1]**2
    g[1][2] = g[2][1] = cell[1]*cell[2]*Num.cos(Num.radians(cell[3]))
    g[2][2] = cell[2]**2
    g_inv = Num.linalg.inv(g)
    return g_inv

def calc_g(cell):
    g = Num.ndarray((3,3),float)
    g[0][0] = cell[0]**2
    g[0][1] = g[1][0] = cell[0]*cell[1]*Num.cos(Num.radians(cell[5]))
    g[0][2] = g[2][0] = cell[0]*cell[2]*Num.cos(Num.radians(cell[4]))
    g[1][1] = cell[1]**2
    g[1][2] = g[2][1] = cell[1]*cell[2]*Num.cos(Num.radians(cell[3]))
    g[2][2] = cell[2]**2
    return g

################################################################################
##################  Rigid Body Calculations  ###################################
class rigid_body:
    def __init__(self):
        self.label = ''
        self.atoms = []
        self.angles = []
        
def rigid_body_rotation(atoms, theta, phi, chi, cell):
    Center  = Num.array([atoms[0][1],atoms[0][2],atoms[0][3]],float)
      
    theta = Num.radians(theta)
    phi = Num.radians(phi)
    chi = Num.radians(chi)
    R_theta = Num.array([[Num.cos(theta),Num.sin(theta),0],\
                         [-Num.sin(theta),Num.cos(theta),0],[0,0,1]],float)
    R_phi   = Num.array([[Num.cos(phi),0,Num.sin(phi)],[0,1,0],\
                         [-Num.sin(phi),0,Num.cos(phi)]],float)
    R_chi   = Num.array([[1,0,0],[0,Num.cos(chi),Num.sin(chi)],\
                         [0,-Num.sin(chi),Num.cos(chi)]],float)
    R  = Num.dot(R_theta,Num.dot(R_phi,R_chi))
    P = Num.array([[cell[0],0,0],[0,cell[1],0],[0,0,cell[2]]],float)
    R = Num.dot(Num.linalg.inv(P),Num.dot(R,P))
    new_atoms = []
    for atom in atoms:
        coords = Num.array([atom[1],atom[2],atom[3]],float) - Center
        coords = Num.dot(R, coords) + Center
        new_atom = [atom[0], coords[0], coords[1], coords[2], atom[4], atom[5],\
                    atom[6], atom[7], atom[8], atom[9], atom[10]]
        new_atoms.append(new_atom)
    return new_atoms

def RB_update(rigid_bodies, surface, parameter, cell):
    surface_new = surface[:]
    for RB in rigid_bodies:
        atoms = []
        for i in RB.atoms:
            atoms.append(surface[i])
        theta = RB.angles[0] * parameter[RB.angles[1]][0]
        phi   = RB.angles[2] * parameter[RB.angles[3]][0]
        chi   = RB.angles[4] * parameter[RB.angles[5]][0]
        atoms = rigid_body_rotation(atoms, theta, phi, chi, cell)
        for i in range(len(RB.atoms)):
            surface_new[RB.atoms[i]] = atoms[i]
    return surface_new
################################################################################
###############   Bond Valence Calculations  ###################################
class BVcluster:
    def __init__(self):
        self.center = int
        self.centerxoffset = int
        self.centeryoffset = int
        self.eqval = int
        self.neighbors = []
        self.neighborsxoffset = []
        self.neighborsyoffset = []
        self.ip = [0.,0.]
        self.r0s = []
        self.bs = []
        self.g = Num.ndarray((3,3),float)

    def calc_BVS(self, surface):
        center = Num.array([surface[self.center][1]+self.centerxoffset,\
                            surface[self.center][2]+self.centeryoffset,\
                            surface[self.center][3]],float)
        BVS = 0
        distances = []
        for i in range(len(self.neighbors)):
            neighbor = Num.array([surface[self.neighbors[i]][1]+\
                                  self.neighborsxoffset[i],\
                                  surface[self.neighbors[i]][2]+\
                                  self.neighborsyoffset[i],\
                                  surface[self.neighbors[i]][3]],float)
            vector = neighbor - center
            dist = (Num.dot(vector,Num.dot(self.g,vector)))**0.5
            distances.append(dist)
            BVS = BVS + Num.exp((self.r0s[i]-dist)/self.bs[i])
        return BVS, distances

def BV_impact(BVclusters, surface):
    impact = 0
    for i in BVclusters:
        BVS, dist = i.calc_BVS(surface)
        eqval = (float(i.eqval) **2)**0.5
        BV_offset = ((BVS - eqval)**2)**0.5 / eqval
        impact = impact + (i.ip[0] * BV_offset)**i.ip[1]
    return impact
################################################################################
######################  Parameter handling  ####################################
def param_unfold(param, param_use, surface, use_bulk_water, use_lay_el):
    if use_bulk_water:
        zwater = param['zwater'][0]
        sig_water = param['sig_water'][0]
        sig_water_bar = param['sig_water_bar'][0]
        d_water = param['d_water'][0]
    else:
        zwater = 0
        sig_water = 0
        sig_water_bar = 0
        d_water = 0

    if use_lay_el:
        occ_el = param['occ_el_0'][0]
        K = param['K'][0]
        sig_el = param['sig_el'][0]
        sig_el_bar = param['sig_el_bar'][0]
        d_el = param['d_el'][0]
        d0_el = param['z_el'][0]
    else:
        occ_el = 0
        K = 0
        sig_el = 0
        sig_el_bar = 0
        d_el = 0
        d0_el = 0
    
    Scale = param['Scale'][0]
    specScale = param['specScale'][0]
    beta = param['beta'][0]
    surface_new = []
    for i in range(len(surface)):
        atom = ['a',0,0,0,0,0,0,0,0,0,0]
        atom[0] = surface[i][0]
        if param_use[i][1] != 'None':
            atom[1] = surface[i][1]+ param_use[i][0]* param[param_use[i][1]][0]
        else: atom[1] = surface[i][1]
        if param_use[i][3] != 'None':
            atom[2] = surface[i][2]+ param_use[i][2]* param[param_use[i][3]][0]
        else: atom[2] = surface[i][2]
        if param_use[i][5] != 'None':
            atom[3] = surface[i][3]+ param_use[i][4]* param[param_use[i][5]][0]
        else: atom[3] = surface[i][3]
        if param_use[i][7] != 'None':
            atom[4] = param_use[i][6]* param[param_use[i][7]][0]
        else: atom[4] = surface[i][4]
        if param_use[i][9] != 'None':
            atom[5] = param_use[i][8]* param[param_use[i][9]][0]
        else: atom[5] = surface[i][5]
        if param_use[i][11] != 'None':
            atom[6] = param_use[i][10]* param[param_use[i][11]][0]
        else: atom[6] = surface[i][6]
        if param_use[i][13] != 'None':
            atom[7] = param_use[i][12]* param[param_use[i][13]][0]
        else: atom[7] = surface[i][7]
        if param_use[i][15] != 'None':
            atom[8] = param_use[i][14]* param[param_use[i][15]][0]
        else: atom[8] = surface[i][8]
        if param_use[i][17] != 'None':
            atom[9] = param_use[i][16]* param[param_use[i][17]][0]
        else: atom[9] = surface[i][9]
        if param_use[i][19] != 'None':
            if param_use[i][18] >= 0: atom[10] = param_use[i][18]* \
                                                 param[param_use[i][19]][0]
            else: atom[10] = 1+ param_use[i][18]* param[param_use[i][19]][0]
        else: atom[10] = surface[i][10]
        surface_new.append(atom)
    global_parms = [occ_el, K,sig_el,sig_el_bar,d_el,d0_el,sig_water, \
                    sig_water_bar, d_water, zwater, Scale, specScale, beta]   
    return  global_parms, surface_new
################################################################################
##################  Structure Factor calculations  #############################
class Fitting_Rod:
    def __init__(self):
        self.H = float
        self.K = float
        self.L = Num.array([],float)
        self.F = Num.array([],float)
        self.Ferr = Num.array([],float)
        self.Lb = Num.array([],float)
        self.Db = float

        self.re_bulk = Num.array([],float)
        self.im_bulk = Num.array([],float)
        
        self.bulk = Num.array([],float)
        self.surf = Num.array([],float)
        self.water = Num.array([],float)
        self.rough = Num.array([],float)
        self.Fcalc = Num.array([],float)
        self.difference = Num.array([],float)
        #####################################################
    def calcFbulk(self, cell, bulk, g_inv, database):
        self.re_bulk = Num.ndarray((len(self.L)),float)
        self.im_bulk = Num.ndarray((len(self.L)),float)
        for i in range(len(self.L)):
            hkl = [self.H,self.K,self.L[i]]
            zeta = self.L[i]+ self.H*cell[6]+ self.K*cell[7]
            re_ctr = 0.5
            im_ctr = -1/(2*Num.tan(Num.pi*zeta))
            re_UC , im_UC = calc_Fuc(hkl,bulk,g_inv,database)
            self.re_bulk[i] = re_ctr*re_UC - im_ctr*im_UC
            self.im_bulk[i] = re_UC*im_ctr + re_ctr*im_UC
####################################################################################################    
def calc_Fuc(hkl,bulk,g_inv,database):
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

def calc_Fsurf(hkl,surface,g_inv,database,exp=Num.exp,low=str.lower,dot=Num.dot,pi=Num.pi,cosinus=Num.cos,sinus = Num.sin,U= Num.ndarray((3,3),float)):
    q_Ang = [hkl[0]*g_inv[0][0]**0.5, hkl[1]*g_inv[1][1]**0.5, hkl[2]*g_inv[2][2]**0.5]
    qd4pi = (dot(dot(hkl,g_inv),hkl))**0.5/4/pi
    a = b = 0
    for atom in surface:
        f_par = database[low(atom[0])]
        U[0][0] = atom[4]
        U[0][1] = U[1][0] = atom[7]*(atom[4]*atom[5])**0.5
        U[0][2] = U[2][0] = atom[8]*(atom[4]*atom[6])**0.5
        U[1][1] = atom[5]
        U[1][2] = U[2][1] = atom[9]*(atom[5]*atom[6])**0.5
        U[2][2] = atom[6]  
        f = (f_par[0]*exp(-(qd4pi)**2*f_par[1]) + f_par[2]*exp(-(qd4pi)**2*f_par[3]) +\
            f_par[4]*exp(-(qd4pi)**2*f_par[5]) + f_par[6]*exp(-(qd4pi)**2*f_par[7]) + f_par[8])*\
            exp(-2* pi**2*(dot(q_Ang,dot(U,q_Ang)))) * atom[10]
        x = 2*pi*(hkl[0]*atom[1] + hkl[1]*atom[2] + hkl[2]*atom[3])
        a = a+(f * cosinus(x))
        b = b+(f * sinus(x))
    return a, b

def calc_Fwater_layered(hkl, sig, sig_bar, d,zwater, g_inv, f_par, cell):
    zinv = g_inv[2][2]**0.5
    q = hkl[2]* zinv
    qd4pi = q/4/Num.pi
    Auc = cell[0]* Num.sin(Num.radians(cell[5]))* cell[1]
    f = Auc * d * 0.033456 * (f_par[0]*Num.exp(-(qd4pi)**2*f_par[1]) + f_par[2]*Num.exp(-(qd4pi)**2*f_par[3]) +\
            f_par[4]*Num.exp(-(qd4pi)**2*f_par[5]) + f_par[6]*Num.exp(-(qd4pi)**2*f_par[7]) + f_par[8])*\
            Num.exp(-2 * Num.pi**2 * q**2 * sig)
    x = Num.pi * q * d
    al = 2 * Num.pi**2 * q**2 * sig_bar
    a = Num.exp(al)*Num.cos(2*x)-1
    b = Num.exp(al)*Num.sin(-2*x)
    c = 4 * Num.cos(x)**2 * Num.sinh(al/2)**2 - 4 * Num.sin(x)**2 * Num.cosh(al/2)**2
    d = -2 * Num.sin(2*x) * Num.sinh(al)
    wert = 2*Num.pi*q*zwater/zinv
    rez = Num.cos(wert)
    imz = Num.sin(wert)
    wert = c**2 + d**2
    relayer = (a*c + b*d)/(wert)
    imlayer = (b*c - a*d)/(wert)
    re = f* (relayer * rez - imlayer * imz)
    im = f* (relayer * imz + imlayer * rez)
    return re, im
        
def calc_F_layered_el(hkl, occ, K, sig, sig_bar, d, d0, g_inv, f_par):
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
        
def calcF(ctr,global_parms,cell,surface,g_inv,NLayers,database, use_bulk_water, RMS_flag, use_lay_el, el):
    occ_el, K,sig_el,sig_el_bar,d_el,d0_el,sig_water,sig_water_bar, d_water,zwater, Scale,specScale,beta = global_parms
    ctr.bulk = Num.ndarray((len(ctr.L)),float)
    ctr.Fcalc = Num.ndarray((len(ctr.L)),float)
    ctr.water = Num.ndarray((len(ctr.L)),float)
    ctr.rough = Num.ndarray((len(ctr.L)),float)
    ctr.surf = Num.ndarray((len(ctr.L)),float)
    exp = Num.exp
    low = str.lower
    dot = Num.dot
    pi = Num.pi
    cosinus = Num.cos
    sinus = Num.sin
    U = Num.ndarray((3,3),float)
    f_par_water = database['o2-.']
    if use_lay_el:
        f_par_el = database[el]
    for i in range(len(ctr.L)):
        hkl = [ctr.H,ctr.K,ctr.L[i]]
        re_surf, im_surf = calc_Fsurf(hkl,surface,g_inv,database,exp,low,dot,pi,cosinus,sinus,U)
            
        if ctr.L[i] > 0:
            n = ctr.Lb[i] + round(ctr.L[i]/ctr.Db) * ctr.Db
        else:
            n = - ctr.Lb[i] + round(ctr.L[i]/ctr.Db) * ctr.Db
        rough = (1-beta)/((1-beta)**2 + 4*beta*sinus(pi*(ctr.L[i] - n)/NLayers)**2)**0.5
           
        if hkl[0] == 0.0 and hkl[1] == 0.0:
            if not use_bulk_water:
                re_water = 0
                im_water = 0
                ctr.water[i] = 0
            else:
                re_water, im_water = calc_Fwater_layered(hkl, sig_water, sig_water_bar, d_water,zwater, g_inv, f_par_water, cell)
                ctr.water[i] = specScale * ((re_water)**2 + (im_water)**2)**0.5
            if not use_lay_el:
                re_el = 0
                im_el = 0
            else:
                re_el, im_el = calc_F_layered_el(hkl,occ_el,K, sig_el, sig_el_bar, d_el, d0_el, g_inv, f_par_el) 
                    
            ctr.bulk[i] = specScale * (ctr.re_bulk[i]**2 + ctr.im_bulk[i]**2)**0.5
            ctr.Fcalc[i] = specScale * rough * ((ctr.re_bulk[i] + re_surf + re_water+ re_el)**2 + (ctr.im_bulk[i] + im_surf + im_water+ im_el)**2)**0.5
              
            ctr.rough[i] = rough * specScale
            ctr.surf[i] = (re_surf**2 + im_surf**2)**0.5 * specScale
        else:
            ctr.bulk[i] = Scale * (ctr.re_bulk[i]**2 + ctr.im_bulk[i]**2)**0.5
            ctr.Fcalc[i] = Scale * rough * ((ctr.re_bulk[i] + re_surf)**2 + (ctr.im_bulk[i] + im_surf)**2)**0.5
            ctr.water[i] = 0
            ctr.rough[i] = rough * Scale
            ctr.surf[i] = (re_surf**2 + im_surf**2)**0.5 * Scale
         
    if RMS_flag == 1:
        ctr.difference = Num.log(ctr.F) - Num.log(ctr.Fcalc)
        for i in range(len(ctr.difference)):
            if ctr.difference[i] <0:
                ctr.difference[i] = -ctr.difference[i]                      
    elif RMS_flag == 2:
        ctr.difference = ctr.F - ctr.Fcalc
        for i in range(len(ctr.difference)):
            if ctr.difference[i] <0:
                ctr.difference[i] = -ctr.difference[i] 
    elif RMS_flag == 3:
        ctr.difference = (ctr.F - ctr.Fcalc) * ctr.Ferr
        for i in range(len(ctr.difference)):
            if ctr.difference[i] <0:
                ctr.difference[i] = -ctr.difference[i] 
    elif RMS_flag == 4:
        ctr.difference = ((ctr.F - ctr.Fcalc)/ctr.Ferr)**2
################################################################################
def calc_CTRs(parameter,param_usage, dat, cell, surface_tmp, NLayers, database,\
              g_inv, Rod_weight, rigid_bodies, use_bulk_water,\
              use_BVC, BVclusters, RMS_flag, use_lay_el, el):

    global_parms, surface_new = param_unfold(parameter,param_usage,surface_tmp,\
                                             use_bulk_water, use_lay_el)
    surface_new = RB_update(rigid_bodies, surface_new, parameter, cell) 
                      
    for ctr in dat:
        calcF(ctr,global_parms,cell,surface_new,g_inv,NLayers,database, \
              use_bulk_water, RMS_flag, use_lay_el, el)

    RMS = 0
    n = 0
    b = 0
    avgF = 0
    for i in range(len(dat)):
        RMS = RMS + Num.sum(dat[i].difference)*Rod_weight[i]
        if RMS_flag == 1:
            n = n + len(dat[i].L)*Rod_weight[i]
            avgF = avgF + Num.sum(Num.log(dat[i].F))*Rod_weight[i]
        elif RMS_flag == 2:
            avgF = avgF + Num.sum(dat[i].F)*Rod_weight[i]
        elif RMS_flag == 3:
            b = b + Num.sum(dat[i].Ferr)*Rod_weight[i]
            n = n + len(dat[i].L)*Rod_weight[i]
            avgF = avgF + Num.sum(dat[i].F)*Rod_weight[i]
        elif RMS_flag == 4:
            avgF = 100
            n = n + len(dat[i].L)*Rod_weight[i]
            
    if RMS_flag == 1:
        RMS = (RMS / n)
        avgF = (avgF / n)
        if avgF < 0: avgF = -avgF
    elif RMS_flag == 3:
        RMS = (RMS / b)
        avgF = avgF / n
    elif RMS_flag == 4:
        RMS = RMS / n

    RMS = (RMS/avgF)*100
        
    if use_BVC:
        impact = BV_impact(BVclusters, surface_new)
        RMS = RMS * (1 + impact)

    return dat, RMS
################################################################################
def calc_Rdata(data):
    n = 0
    sumlogFerr = 0
    sumlogF = 0
    sumF = 0
    sumFerr = 0
    sumFerr2 = 0
    for rod in data:
        n = n + len(rod.F)
        sumF = sumF + Num.sum(rod.F)
        sumFerr = sumFerr + Num.sum(rod.Ferr)
        sumFerr2 = sumFerr2 + Num.sum(rod.Ferr**2)
        sumlogF = sumlogF + ((Num.sum(Num.log(rod.F)))**2)**0.5
        sumlogFerr = sumlogFerr + ((Num.sum(Num.log(1+rod.Ferr/rod.F)))**2)**0.5
    Rdata1 = sumlogFerr/sumlogF * 100
    Rdata2 = sumFerr / sumF * 100
    Rdata3 = sumFerr2*n / (sumFerr * sumF) * 100
    return ['(:P)', str(round(Rdata1, 2)),str(round(Rdata2, 2)),\
            str(round(Rdata3, 2)), '1.0']
################################################################################
##############################  reading datafiles  #############################
def read_bulk(bulkfile):
    bulk=[]
    cell = []
    Nlayers = 1
    f = open(bulkfile, 'r')
    data = f.readlines()
    f.close()
    for i in data:
        tmp = str.rsplit(i)
        if len(tmp) == 10 and tmp[0] == 'cell':
            for j in range(8):
                cell.append(float(tmp[j+1]))
            Nlayers = int(tmp[9])
        else:
            tmp2 = [tmp[0], float(tmp[1]),float(tmp[2]),float(tmp[3]),\
                    float(tmp[4])]
            bulk.append(tmp2)
    return bulk, cell, Nlayers

def read_surface(surfacefile,database):
    surface=[]
    f = open(surfacefile, 'r')
    data = f.readlines()
    f.close()
    for i in data:
        tmp = str.rsplit(i)
        tmp2 = [tmp[0],float(tmp[1]),float(tmp[2]),float(tmp[3]),float(tmp[4]),\
                float(tmp[5]),float(tmp[6]),float(tmp[7]),float(tmp[8]),\
                float(tmp[9]),float(tmp[10])]
        surface.append(tmp2)
        
    parameter_usage=[]
    for i in data:
        tmp = str.rsplit(i)
        tmp2 = [float(tmp[11]),tmp[12],float(tmp[13]),tmp[14],float(tmp[15]),\
                tmp[16],float(tmp[17]),tmp[18],float(tmp[19]),tmp[20],\
                float(tmp[21]),tmp[22],float(tmp[23]),tmp[24],float(tmp[25]),\
                tmp[26],float(tmp[27]),tmp[28],float(tmp[29]),tmp[30]]
        parameter_usage.append(tmp2)
    runningDB = {}
    runningDB['o2-.'] = database['o2-.']
    for i in range(len(surface)):
        if surface[i][0] not in runningDB.keys():
            key = str.lower(surface[i][0])
            runningDB[key]=database[key]
        
    return surface, parameter_usage, runningDB

def read_data(datafile):
    f = open(datafile, 'r')
    data = f.readlines()
    f.close()
    len_dat = 7
    dat = Num.ndarray((len(data),len_dat),float)
    
    z=0
    for i in range(len(data)):
        tmp = str.rsplit(data[i])
        if tmp[0] == '%':
            dat[i] = [0,0,0,0,0,0,0]
            z = z+1
        else:
            for j in range(len(tmp)):
                dat[i][j] = float(tmp[j])

    dat1 = Num.ndarray((len(data)-z,len_dat),float)
    x = 0
    for i in range(len(data)):
        if dat[i][2] != 0:
            dat1[i-x] = dat[i]
        else:
            x = x+1

    dat=[]
    tmp = Fitting_Rod()
    tmp.H = dat1[0][0]
    tmp.K = dat1[0][1]
    tmp.Db = dat1[0][6]
    for i in range(len(data)-z):
        one_rod = True
        if i>1:
            if(dat1[i][0]!=dat1[i-1][0]) or (dat1[i][1]!=dat1[i-1][1]):
                one_rod = False        
        if one_rod:
            tmp.L = Num.append(tmp.L,dat1[i][2])
            tmp.F = Num.append(tmp.F,dat1[i][3])
            tmp.Ferr = Num.append(tmp.Ferr,dat1[i][4])
            tmp.Lb = Num.append(tmp.Lb,dat1[i][5])
        else:
            dat.append(tmp)
            tmp = Fitting_Rod()
            tmp.H = dat1[i][0]
            tmp.K = dat1[i][1]
            tmp.Db = dat1[i][6]
                
            tmp.L = Num.append(tmp.L,dat1[i][2])
            tmp.F = Num.append(tmp.F,dat1[i][3])
            tmp.Ferr = Num.append(tmp.Ferr,dat1[i][4])
            tmp.Lb = Num.append(tmp.Lb,dat1[i][5])
                    
    dat.append(tmp)
    return dat

def read_parameters(parameterfile):
    parameter={}
    param_labels = []
    f = open(parameterfile, 'r')
    data = f.readlines()
    f.close()
    for i in data:
        tmp = str.rsplit(i)
        if tmp[0] != '%':
            parameter[tmp[0]]= [float(tmp[1]),float(tmp[2]),float(tmp[3]),False]
            if tmp[4] == 'True':
                parameter[tmp[0]][3] = True
            param_labels.append(tmp[0])

    return parameter, param_labels

def read_rigid_bodies(rigidbodyfile):
    rigid_bodies=[]
    f = open(rigidbodyfile, 'r')
    data = f.readlines()
    f.close()
    for i in data:
        tmp = str.rsplit(i)
        if tmp[0] != '%':
            tmp2 = rigid_body()
            tmp2.label = tmp[0]
            atoms = []
            ct = int(tmp[1])
            for a in range(ct):
                atoms.append(int(tmp[a+2]))
            tmp2.atoms = atoms
            tmp2.angles = [float(tmp[ct+2]),tmp[ct+3],float(tmp[ct+4]),\
                           tmp[ct+5],float(tmp[ct+6]),tmp[ct+7]]
            rigid_bodies.append(tmp2)
            del(tmp2)
        
    return rigid_bodies

def read_BV(BVfile, cell):
    BVclusters = []
    f = open(BVfile,'r')
    data = f.readlines()
    f.close()
    for i in data:
        BVC = BVcluster()
        tmp = str.rsplit(i)
        n = int(tmp[0])
        center_label = tmp[1]
        BVC.eqval = int(tmp[2])
        BVC.center = int(tmp[3])
        BVC.centerxoffset = int(tmp[4])
        BVC.centeryoffset = int(tmp[5])
        for j in range(n-1):
            neighbor_label = tmp[6+j*5]
            neighbor_valence = int(tmp[7+j*5])
            key = center_label + str(BVC.eqval) + neighbor_label + \
                  str(neighbor_valence)
            if key in bv_params.keys():
                BVC.r0s.append(bv_params[key][0])
                BVC.bs.append(bv_params[key][1])
            else:
                key = neighbor_label + str(neighbor_valence) + center_label +\
                      str(BVC.eqval)
                if key in bv_params.keys():
                    BVC.r0s.append(bv_params[key][0])
                    BVC.bs.append(bv_params[key][1])
                else:
                    print ' No bond valence parameters found for: '\
                      +center_label+' '+str(BVC.eqval)+' and '+\
                      neighbor_label+' '+str(neighbor_valence)
            BVC.neighbors.append(int(tmp[8+j*5]))
            BVC.neighborsxoffset.append(int(tmp[9+j*5]))
            BVC.neighborsyoffset.append(int(tmp[10+j*5]))
        BVC.ip[0] = float(tmp[11+(n-2)*5])
        BVC.ip[1] = float(tmp[12+(n-2)*5])
        BVC.g = calc_g(cell)
        BVclusters.append(BVC)
    return BVclusters
################################################################################
########################  writing files  #######################################
def write_surface(cell, surface,param,param_use, rigid_bodies, use_bulk_water, use_lay_el, filename = 'surface.sur'):
    global_parms, surface = param_unfold(param,param_use, surface, use_bulk_water, use_lay_el)
    surface = RB_update(rigid_bodies, surface, param, cell)
    f = file(filename, 'w')
    for atom in surface:
        line = "%5s %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f\n" % (atom[0],atom[1],atom[2],atom[3],\
                                                                                      atom[4],atom[5],atom[6],atom[7],atom[8],atom[9],atom[10])
        f.write(line)    

def write_cif(cell4, surface4,param4,param_use, rigid_bodies, use_bulk_water, use_lay_el, filename = 'surface.cif'):
    global_parms, surface4 = param_unfold(param4,param_use, surface4, use_bulk_water, use_lay_el)
    surface4 = RB_update(rigid_bodies, surface4, param4, cell4)
    
    f = file(filename, 'w')
    f.write('data_global\n')
    f.write('loop_\n')
    f.write('_cell_length_a  '+str(cell4[0])+'\n')
    f.write('_cell_length_b  '+str(cell4[1])+'\n')
    f.write('_cell_length_c  '+str(cell4[2])+'\n')
    f.write('_cell_angle_alpha  '+str(cell4[3])+'\n')
    f.write('_cell_angle_beta  '+str(cell4[4])+'\n')
    f.write('_cell_angle_gamma  '+str(cell4[5])+'\n')
    f.write('_symmetry_space_group_name_H-M "P 1"\n')
    f.write('loop_\n')
    f.write('_space_group_symop_operation_xyz\n')
    f.write(' "x,y,z"\n')
    f.write('loop_\n')
    f.write('_atom_site_label\n')
    f.write('_atom_site_fract_x\n')
    f.write('_atom_site_fract_y\n')
    f.write('_atom_site_fract_z\n')
    for i in range(len(surface4)):
        f.write(str(surface4[i][0])+str(i+1)+'  '+str(surface4[i][1])+'  '+str(surface4[i][2])+'  '+str(surface4[i][3])+'\n')
        f.write(str(surface4[i][0])+str(len(surface4)+i+1)+'  '+str(surface4[i][1]+1)+'  '+str(surface4[i][2])+'  '+str(surface4[i][3])+'\n')
        f.write(str(surface4[i][0])+str(2*len(surface4)+i+1)+'  '+str(surface4[i][1])+'  '+str(surface4[i][2]+1)+'  '+str(surface4[i][3])+'\n')
        f.write(str(surface4[i][0])+str(3*len(surface4)+i+1)+'  '+str(surface4[i][1]+1)+'  '+str(surface4[i][2]+1)+'  '+str(surface4[i][3])+'\n')
    f.write('loop_\n')
    f.write('_atom_site_aniso_label\n')
    f.write('_atom_site_aniso_U_11\n')
    f.write('_atom_site_aniso_U_22\n')
    f.write('_atom_site_aniso_U_33\n')
    f.write('_atom_site_aniso_U_12\n')
    f.write('_atom_site_aniso_U_13\n')
    f.write('_atom_site_aniso_U_23\n')
    for i in range(len(surface4)):
        uxy = surface4[i][7] * (surface4[i][4])**0.5 * (surface4[i][5])**0.5
        uxz = surface4[i][8] * (surface4[i][4])**0.5 * (surface4[i][6])**0.5
        uyz = surface4[i][9] * (surface4[i][5])**0.5 * (surface4[i][6])**0.5
        f.write(str(surface4[i][0])+str(i+1)+'  '+str(surface4[i][4])+'  '+str(surface4[i][5])+'  '+str(surface4[i][6])+'  '+str(uxy)+'  '+str(uxz)+'  '+str(uyz)+'\n')
        f.write(str(surface4[i][0])+str(len(surface4)+i+1)+'  '+str(surface4[i][4])+'  '+str(surface4[i][5])+'  '+str(surface4[i][6])+'  '+str(uxy)+'  '+str(uxz)+'  '+str(uyz)+'\n')
        f.write(str(surface4[i][0])+str(2*len(surface4)+i+1)+'  '+str(surface4[i][4])+'  '+str(surface4[i][5])+'  '+str(surface4[i][6])+'  '+str(uxy)+'  '+str(uxz)+'  '+str(uyz)+'\n')
        f.write(str(surface4[i][0])+str(3*len(surface4)+i+1)+'  '+str(surface4[i][4])+'  '+str(surface4[i][5])+'  '+str(surface4[i][6])+'  '+str(uxy)+'  '+str(uxz)+'  '+str(uyz)+'\n')
    f.close()

def write_par(parameter, param_labels, filename = 'parameters.new'):
    f = file(filename, 'w')
    f.write('% param_label              start          min          max    '+\
            'refine_flag\n')
    for i in param_labels:
        line = "%13s %18.12f %12.4f %12.4f %10s\n" % (i,parameter[i][0],
                                                      parameter[i][1],
                                                      parameter[i][2],
                                                      parameter[i][3])
        f.write(line)
    f.close()

def write_data(data, filename = 'result.dat'):
    f = file(filename, 'w')
    f.write('    H     K        L           F        Ferr       Fcalc       '\
            +'Fbulk       Fsurf      Frough      Fwater\n')
    for i in data:
        for j in range(len(i.L)):
            line = "%5.2f %5.2f %8.4f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n" % (i.H, i.K, i.L[j], i.F[j], i.Ferr[j], i.Fcalc[j],
                                                                                             i.bulk[j], i.surf[j], i.rough[j], i.water[j])                                                                                        
            f.write(line)
    f.close()

################################################################################
##############################  Simulated Annealing  ###########################
def simulated_annealing(dat, cell, NLayers, surface, database, parameter, \
                        param_usage, rigid_bodies, panel):

    Tstart,Tend,cool,maxrun,MC,factor,random_parameters = panel.sim_an_params
    g_inv = calc_g_inv(cell)
    Rod_weight = panel.Rod_weight
    use_bulk_water = panel.UBW_flag
    use_lay_el = panel.use_lay_el
    el = panel.el
    use_BVC = panel.use_BVC
    BVclusters = panel.BVclusters
    RMS_flag = panel.RMS_flag
    fig1 = panel.Figure1
    fig3 = panel.Figure3
    plot_dims = panel.plotdims
    plot_bulk = panel.doplotbulk
    plot_surf = panel.doplotsurf
    plot_rough = panel.doplotrough
    plot_water = panel.doplotwater
    statusbar = panel.nb.frame.statusbar
    
    def plot_R(plot3,R_track):
        plot3.plot(range(len(R_track)),R_track, 'ro')
        plot3.figure.canvas.draw()
        
    if random_parameters:
        for i in parameter.keys():
            if parameter[i][3]:
                parameter[i][0] = random.uniform(parameter[i][1],\
                                                 parameter[i][2])
      
    dat, RMS = calc_CTRs(parameter, param_usage, dat, cell, surface, NLayers,\
                         database, g_inv,Rod_weight, rigid_bodies, \
                         use_bulk_water, use_BVC, BVclusters, RMS_flag,\
                         use_lay_el, el)
    guess = (int(Num.log(float(Tend)/float(Tstart))/Num.log(cool))+1)*maxrun
    print '\n approximated number of iterations in Simulated Annealing: '\
          +str( int(guess) )
    print 'R start = '+str(RMS)
    RMS_best = RMS
    param_best = parameter
    R_track =Num.array([RMS],float)
    fig1 = plot_rods(fig1, dat, plot_dims, plot_bulk, plot_surf, plot_rough,\
                     plot_water, RMS)
    fig1.canvas.draw()
                    
    #counters
    Random = 0
    better = 0
    rejected = 0
    if fig3 == None:
        fig3 = figure(3, figsize = [9,6])
    else:
        fig3.clear()
    fig3.suptitle('R track', fontsize = 20)
    plot3 = None
    plot3 = fig3.add_subplot(111)
    plot_R(plot3,R_track)
    while Tstart > Tend:
        z = 0
        while (z < maxrun):
            while wx.GetApp().Pending():
                wx.GetApp().Dispatch()
                wx.GetApp().Yield(True)
            if panel.StopFit:
                z= maxrun
                Tstart = Tend
            RMS_tmp = 0.
            param_tmp = {}
            #permutation of parameters 
            for i in parameter.keys():
                if parameter[i][3]:
                    param_tmp[i] =  [parameter[i][0] + random.uniform((parameter[i][1]-parameter[i][0])*(MC * Tstart), (parameter[i][2]-parameter[i][0])*(MC * Tstart))]
                    if param_tmp[i][0] < parameter[i][1]: param_tmp[i][0] = parameter[i][1]
                    if param_tmp[i][0] > parameter[i][2]: param_tmp[i][0] = parameter[i][2]
                else:
                    param_tmp[i] =  [parameter[i][0]]
   
            dat, RMS_tmp = calc_CTRs(param_tmp,param_usage, dat, cell,surface, NLayers, database, g_inv,\
                                     Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag,\
                                     use_lay_el, el)
            dR = (RMS - RMS_tmp) * factor

            if dR > 0:
                for i in parameter.keys():
                    parameter[i][0] = param_tmp[i][0]
                z = z+1
                better = better+1
                R_track = Num.append(R_track,RMS_tmp)
                RMS = RMS_tmp
                plot_R(plot3,R_track)
                fig1 = plot_rods(fig1, dat,plot_dims, plot_bulk, plot_surf, plot_rough,plot_water, RMS)
                fig1.canvas.draw()
                statusbar.SetStatusText('better    '+str(round(RMS_tmp,5)),0)
                statusbar.SetStatusText('T: '+str(round(Tstart,1))+', iteration: '+str(z),1) 
                if RMS < RMS_best:
                    RMS_best = RMS
                    for i in parameter.keys():
                        param_best[i][0] = parameter[i][0] 
            elif dR <= 0:
                Boltz = exp(dR / Tstart)
                Rand = random.uniform(0,1)
                if(Boltz > Rand):
                    for i in parameter.keys():
                        parameter[i][0] = param_tmp[i][0]
                    z = maxrun
                    Random = Random+1
                    R_track = Num.append(R_track,RMS_tmp)
                    plot_R(plot3,R_track)
                    param_track.append(parameter)
                    RMS = RMS_tmp
                    statusbar.SetStatusText('random    '+str(round(RMS_tmp,5)),0)
                    statusbar.SetStatusText('T: '+str(round(Tstart,1))+', iteration: '+str(z),1)
                elif Boltz <= Rand:
                    z = z+1
                    rejected = rejected+1
                    statusbar.SetStatusText('rejected    '+str(round(RMS_tmp,5)),0)
                    statusbar.SetStatusText('T: '+str(round(Tstart,1))+', iteration: '+str(z),1)
        Tstart = Tstart * cool

    print '****************************'
    print 'Number of cycles: '+str(Random+ better+ rejected)
    print 'random: ' + str(Random)
    print 'better: ' + str(better)
    print 'rejected: '+ str(rejected)

    data_best, RMS_best = calc_CTRs(param_best,param_usage, dat, cell,surface, NLayers, database,\
                                    g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC,\
                                    BVclusters, RMS_flag, use_lay_el, el)
    statusbar.SetStatusText('End of Simulated Annealing, best R: '+str(round(RMS_best,5)),0)
    statusbar.SetStatusText('',1)
    print '\n####################################################\n'
    print 'the best fit R = '+str(RMS_best)+'\n'
    print '*************************************\n'

    return data_best, param_best, RMS_best
################################################################################
############################  Plotting  ########################################
def plot_rods(fig1, dat, plot_dims, plot_bulk, plot_surf, plot_rough,\
              plot_water, RMS):
    if fig1 == None:
        fig1 = figure(1, figsize = [15,9])
    fig1.clear()
    fig1.suptitle('R = '+str(round(RMS,7)), fontsize=20)
    subplots = []
    for i in range(len(dat)):
        pl = str(plot_dims[0])+str(plot_dims[1])+str(i+1)
        subplots.append(fig1.add_subplot(pl))
        tmp = dat[i]
        if plot_bulk: subplots[i].plot(tmp.L,tmp.bulk,'g')
        if plot_surf: subplots[i].plot(tmp.L,tmp.surf,'r')
        if plot_water: subplots[i].plot(tmp.L,tmp.water,'m')
        if plot_rough: subplots[i].plot(tmp.L,tmp.rough,'c')
        subplots[i].plot(tmp.L,tmp.Fcalc,'k')
        subplots[i].errorbar(tmp.L,tmp.F,tmp.Ferr, fmt = 'bo')
        subplots[i].set_title(str(int(tmp.H))+str(int(tmp.K))+'L')
        subplots[i].semilogy()
    return fig1

def plot_edensity(Fig, surface, param, param_use, cell, database, rigid_bodies,\
                  use_bulk_water, use_lay_el, resel, filename = ''):
    edens = Num.zeros((1000),float)
    global_parms, surface = \
            param_unfold(param, param_use, surface, use_bulk_water, use_lay_el)
    occ_el, K,sig_el,sig_el_bar,d_el,d0_el,\
            sig,sig_bar, d_water,zwater, Scale,specScale,beta = global_parms
    surface = RB_update(rigid_bodies, surface, param, cell)
    zmin = 0.
    zmax = 0.
    Resonant = False
    resedens = None
    for i in range(len(surface)):
        if surface[i][3] > zmax: zmax = surface[i][3]
        if surface[i][3] < zmin: zmin = surface[i][3]
        if surface[i][0] == resel: Resonant = True
    if Resonant or use_lay_el:
        resedens = Num.zeros((1000),float)
    zmin = (zmin - 0.1)*cell[2]
    zmax = (zmax + 0.75)*cell[2]
    abscissa = Num.arange(1000)*(zmax-zmin)/999 + zmin
    Auc = cell[0]* Num.sin(Num.radians(cell[5]))* cell[1]
    if use_bulk_water:
        d = zwater*cell[2]
        while d < zmax+cell[2]:
            for i in range(len(abscissa)):
                edens[i] = edens[i] + ((2*Num.pi*sig)**(-1.5)* \
                                       Num.exp(-0.5*(abscissa[i]-d)**2/sig))\
                                       *2*Num.pi*sig * 0.33456*d_water
            d = d+d_water
            sig = sig+ sig_bar
        
    for x in surface:
        f_par = database[str.lower(x[0])]
        f = (f_par[0] + f_par[2] +f_par[4]+ f_par[6]+ f_par[8]) * x[10]
        for i in range(len(abscissa)):
            edens[i] = edens[i] + ((2*Num.pi*x[6])**(-1.5)*\
                                   Num.exp(-0.5*(abscissa[i]-\
                                   cell[2]*x[3])**2/x[6]))*\
                                   f*2*Num.pi*x[6]/Auc
        if Resonant:
            if x[0] == resel:
                f_par = database[str.lower(resel)]
                f = (f_par[0] + f_par[2] +f_par[4]+ f_par[6]+ f_par[8]) * x[10]
                for i in range(len(abscissa)):
                    resedens[i] = resedens[i] + ((2*Num.pi*x[6])**(-1.5)*\
                                                 Num.exp(-0.5*(abscissa[i]\
                                                 -cell[2]*x[3])**2/x[6]))\
                                                 * f*2*Num.pi*x[6]/Auc
    if use_lay_el:
        layresedens = Num.zeros((1000),float)
        f_par = database[str.lower(resel)]
        f = (f_par[0] + f_par[2] +f_par[4]+ f_par[6]+ f_par[8])
        d = d0_el*cell[2]
        while d < zmax+cell[2]:
            for i in range(len(abscissa)):
                layresedens[i] = layresedens[i] + ((2*Num.pi*sig_el)**(-1.5)* \
                                Num.exp(-0.5*(abscissa[i]-d)**2/sig_el))\
                                *2*Num.pi*sig_el* f * occ_el/Auc
            d = d+d_el
            sig_el = sig_el+ sig_el_bar
            occ_el = occ_el * Num.exp(-K*(d-d0_el*cell[2]))
        edens = edens + layresedens
        resedens = resedens + layresedens
                
    
    if Fig == None:
        Fig = figure(2)
    else: Fig.clear()    
    edensplot = Fig.add_subplot(111)
    edensplot.plot(abscissa, edens, 'r')
    if Resonant or use_lay_el: edensplot.plot(abscissa, resedens, 'b')
    edensplot.set_xlabel('z [Angstroem]')
    edensplot.set_ylabel('electron density [Angstroem**(-3)]')

    if filename != '':
        f = file(filename, 'w')
        if resedens != None:
            f.write('% z            total_e-density '+\
                    'resonant_element_e-density\n')
            for i in range(len(abscissa)):
                line = "%15.5f %15.5f %15.5f\n" % (abscissa[i], edens[i],
                                                   resedens[i])
                f.write(line)
        else:
            f.write('% z           total_e-density\n')
            for i in range(len(abscissa)):
                line = "%15.5f %15.5f\n" % (abscissa[i], edens[i])
                f.write(line)
        f.close()
    return Fig
################################################################################
######################  CHECKS  ################################################
def check_model_consistency(param_labels, parameter, parameter_usage,\
                            rigid_bodies, use_bulk_water):
    used_params = []
    for i in param_labels:
        if parameter[i][3]: used_params.append(i)
    if use_bulk_water:
        global_parameters = ['z_el','d_el','sig_el','sig_el_bar','K',\
                             'occ_el_0','zwater','sig_water','Scale',\
                             'specScale','beta']
    else:
        global_parameters = ['Scale','specScale','beta']
    global_flag = True
    for i in global_parameters:
        flag = True
        for j in param_labels:
            if i == j:
                flag = False
                break
        if flag:
            print 'Global parameter: '+i+' is not defined in the parameter list'
            global_flag = False
            
    for i in used_params:
        flag = True
        for j in global_parameters:
            if i == j: flag = False
        if flag:
            for j in range(len(parameter_usage)):
                for k in range(len(parameter_usage[j])):
                    if i == parameter_usage[j][k]:
                        flag = False
                        break
        if flag:
            for j in rigid_bodies:
                for k in j.angles:
                    if i == k:
                        flag = False
                        break
        if flag:
            print 'Model Inconsistency: parameter '+i+\
                  ' is used in the fit but is not defined in the model'
            global_flag = False
    return global_flag

def check_vibes(surface, param, param_use):
    global_parms, surface = \
            param_unfold(param, param_use, surface, False, False)
    U = Num.ndarray((3,3),float)
    i = 0
    bad_indices = []
    bad_labels = []
    for atom in surface:
        U[0][0] = atom[4]
        U[0][1] = U[1][0] = atom[7]*(atom[4]*atom[5])**0.5
        U[0][2] = U[2][0] = atom[8]*(atom[4]*atom[6])**0.5
        U[1][1] = atom[5]
        U[1][2] = U[2][1] = atom[9]*(atom[5]*atom[6])**0.5
        U[2][2] = atom[6]
        if Num.linalg.det(U) <= 0:
            bad_indices.append(str(i))
            bad_labels.append(atom[0])
        i = i+1
    if len(bad_indices) > 0:
        message = 'BAD VIBRATION ALERT: vibrational ellipsoids of the'+\
                  ' following atoms have negative volumes: \n'
        for i in range(len(bad_indices)):
            message = message + 'atom '+bad_indices[i] + \
                      ' (' + bad_labels[i] + '), '
        print message
    return
################################################################################
############################### Fourier Mapping ################################
class Fourier_Rod:
    def __init__(self):
        self.H = float
        self.K = float
        self.L = []
        self.A = []
        self.P = []
        self.Lb = []
        self.Db = 2

def create_F_data(dat):
    F_data = []
    for CTR in dat:
        Rod = Fourier_Rod()
        Rod.H = float(CTR.H)
        Rod.K = float(CTR.K)
        for i in range(len(CTR.L)):
            Rod.L.append(CTR.L[i])
            Rod.A.append(CTR.F[i] / CTR.rough[i])
            Rod.P.append(0)
        F_data.append(Rod)
    return F_data

def create_F_bulk(dat):
    F_bulk = []
    for CTR in dat:
        Rod = Fourier_Rod()
        Rod.H = float(CTR.H)
        Rod.K = float(CTR.K)
        for i in range(len(CTR.L)):
            Rod.L.append(CTR.L[i])
            re = CTR.re_bulk[i]
            im = CTR.im_bulk[i]
            P = Num.arctan(im/re)/(2*Num.pi)
            Rod.A.append(re/Num.cos(2*Num.pi*P))
            Rod.P.append(P)
        F_bulk.append(Rod)
    return F_bulk

def create_F_model(dat, surface,g_inv,cell, database,parameter, param_usage,\
                   rigid_bodies, use_bulk_water, use_lay_el):
    F_model = []
    global_parms, surface_new = param_unfold(parameter,param_usage, surface,\
                                             use_bulk_water, use_lay_el)
    surface_new = RB_update(rigid_bodies, surface_new, parameter, cell)
    occ_el, K,sig_el,sig_el_bar,d_el,d0_el,sig_water,sig_water_bar, d_water,\
            zwater, Scale,specScale,beta = global_parms

    for CTR in dat:
        Rod = Fourier_Rod()
        Rod.H = float(CTR.H)
        Rod.K = float(CTR.K)
        f_par = database['o2-.']
        for i in range(len(CTR.L)):
            Rod.L.append(CTR.L[i])
            hkl = Num.array([CTR.H,CTR.K,CTR.L[i]],float)
            re_surf, im_surf = calc_Fsurf(hkl,surface_new,g_inv,database)
            if use_bulk_water and hkl[0]==0 and hkl[1]==0:
                re_water, im_water = calc_Fwater_layered(hkl, sig_water,\
                                                         sig_water_bar, d_water\
                                                         ,zwater, g_inv, \
                                                         f_par, cell)
            else:
                re_water = 0
                im_water = 0
            re = CTR.re_bulk[i] + re_surf + re_water
            im = CTR.im_bulk[i] + im_surf + im_water
            P = Num.arctan(im/re)/(2*Num.pi)
            Rod.A.append(re/Num.cos(2*Num.pi*P))
            Rod.P.append(P)
        F_model.append(Rod)
    return F_model

def create_F_obs(dat,surface,g_inv,cell, database,parameter, param_usage,\
                   rigid_bodies, use_bulk_water, use_lay_el, NLayers):
    F_model = create_F_model(dat, surface,g_inv,cell, database,parameter,\
                             param_usage,rigid_bodies, use_bulk_water,\
                             use_lay_el)
    F_obs = create_F_data(dat)
    for i in range(len(F_obs)):
        F_obs[i].P = F_model[i].P
    return F_obs

def fdiff2rho(Fk, Fu, cell, xf, yf, zf, an, bn, cn, zmin, flag, parent):
    """
    calculate the difference fourier map from two sets of given
    structure factors, one with known phases (Fk)
    the other one with unknown phases (Fu)
    """          
        
    rho = Num.zeros((an,bn,cn), float)
    sampx = cell[0]* xf
    sampy = cell[1]* yf
    sampz = cell[2]* zf
    g_inv = calc_g_inv(cell)
    g_inv2 = calc_g_inv([sampx,sampy,sampz,cell[3],cell[4],cell[5]])
    V = Num.linalg.det(Num.linalg.inv(g_inv2))**0.5
    maxpro = int(an*bn*cn)
    dialog = wx.ProgressDialog("Progress of Fourier Synthesis",\
                               "Please be patient",maxpro,parent,\
                               style=wx.PD_CAN_ABORT|wx.PD_ELAPSED_TIME|\
                               wx.PD_AUTO_HIDE)
    keepgoing = True
    count = 0
    for i in range(int(an)):
        x = sampx/an * float(i)
        for j in range(int(bn)):
            y = sampy/bn * float(j)
            for k in range(int(cn)):
                z = sampz/cn * float(k) + sampz/zf *zmin
                R = Num.array([x,y,z],float)
                for r in range(len(Fk)):
                    while wx.GetApp().Pending():
                        wx.GetApp().Dispatch()
                        wx.GetApp().Yield(True)
                    for l in range(len(Fk[r].L)):
                        Q = Num.array([Fk[r].H * g_inv[0][0]**0.5,\
                                        Fk[r].K * g_inv[1][1]**0.5,\
                                        Fk[r].L[l]* g_inv[2][2]**0.5],float)
                        Ak = Fk[r].A[l]
                        P = Fk[r].P[l]
                        if flag < 3:
                            rho[i][j][k] = rho[i][j][k] + Ak\
                                       *Num.cos(2*Num.pi*(P - Num.dot(R,Q)))
                        else:
                            rho[i][j][k] = rho[i][j][k] + (Fu[r].A[l] - Ak)\
                                *Num.cos(2*Num.pi*(P - Num.dot(R,Q)))
                keepgoing = dialog.Update(count)
                if not keepgoing:
                    rho = Num.zeros((an,bn,cn), float)
                    break
                count = count+1        
    rho = rho/(2*Num.pi*V)
    dialog.Destroy()
    return rho, [sampx,sampy,sampz]
################################################################################
######################  PSYCHO IMPORT  #########################################
try:
    import psyco
    psyco.bind(calc_Fsurf)
    psyco.bind(fdiff2rho)
    psyco.bind(simulated_annealing)
except ImportError:
    print '\n COULD NOT FIND PSYCO IN YOUR PYTHON SIDEPACKAGES!!'
    print 'THE USE OF PSYCO SPEEDS UP PISURF BY ~25%!'
    pass



    
    
    



