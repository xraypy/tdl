"""
Functions used by rasd_menu to analyze rasd data

Authors/Modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""

import numpy as num
from pylab import *
from scipy.optimize import leastsq
import random
from atomic import f0data as database

###################################### calculations ####################################################
def calc_g_inv(cell):
    g = num.ndarray((3,3),float)
    g[0][0] = cell[0]**2
    g[0][1] = cell[0]*cell[1]*num.cos(num.radians(cell[5]))
    g[0][2] = cell[0]*cell[2]*num.cos(num.radians(cell[4]))
    g[1][0] = cell[1]*cell[0]*num.cos(num.radians(cell[5]))
    g[1][1] = cell[1]**2
    g[1][2] = cell[1]*cell[2]*num.cos(num.radians(cell[3]))
    g[2][0] = cell[2]*cell[0]*num.cos(num.radians(cell[4]))
    g[2][1] = cell[2]*cell[1]*num.cos(num.radians(cell[3]))
    g[2][2] = cell[2]**2

    g_inv = num.linalg.inv(g)
    return g_inv

def calc_Fuc(hkl,bulk,g_inv,database=database):
    a = 0
    b = 0
    for i in range(shape(bulk)[0]):
        el = str.lower(bulk[i][0])
        f_par = database[el]
        q = (num.dot(num.dot(hkl,g_inv),hkl))**0.5
        f = (f_par[0]*num.exp(-(q/4/num.pi)**2*f_par[1]) + f_par[2]*num.exp(-(q/4/num.pi)**2*f_par[3]) +\
            f_par[4]*num.exp(-(q/4/num.pi)**2*f_par[5]) + f_par[6]*num.exp(-(q/4/num.pi)**2*f_par[7]) + f_par[8])*\
            num.exp(-(q/2)**2 *bulk[i][4])
        a = a + (f * num.cos(2*num.pi*(hkl[0]*bulk[i][1] + hkl[1]*bulk[i][2] + hkl[2]*bulk[i][3])))
        b = b + (f * num.sin(2*num.pi*(hkl[0]*bulk[i][1] + hkl[1]*bulk[i][2] + hkl[2]*bulk[i][3])))
    return a, b

def calc_Fsurf(hkl,surface,g_inv,database=database):
    a = 0
    b = 0
    for i in range(shape(surface)[0]):
        el = str.lower(surface[i][0])
        f_par = database[el]
        q = (num.dot(num.dot(hkl,g_inv),hkl))**0.5
        f = (f_par[0]*num.exp(-(q/4/num.pi)**2*f_par[1]) + f_par[2]*num.exp(-(q/4/num.pi)**2*f_par[3]) +\
            f_par[4]*num.exp(-(q/4/num.pi)**2*f_par[5]) + f_par[6]*num.exp(-(q/4/num.pi)**2*f_par[7]) + f_par[8])*\
            num.exp(-(q/2)**2 *surface[i][5])* surface[i][4]
        a = a + (f * num.cos(2*num.pi*(hkl[0]*surface[i][1] + hkl[1]*surface[i][2] + hkl[2]*surface[i][3])))
        b = b + (f * num.sin(2*num.pi*(hkl[0]*surface[i][1] + hkl[1]*surface[i][2] + hkl[2]*surface[i][3])))
    return a, b
##########################################################################################
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
        self.E = num.array([], float)
        self.f1 = num.array([], float)
        self.f2 = num.array([], float)
        self.cell = []
        self.g_inv = num.array([[]],float)
        #Parameters for Fourier synthesis
        self.Fourier = []
        self.ZR = 1
        self.xf = 1
        self.yf = 1
        self.zf = 1
        self.an = 11
        self.bn = 11
        self.cn = 61
        self.Plusminus = False
        #Parameters for Structure Refinement
        self.reflist = []
        self.RMS = 0.
        self.ndata = 0
        self.natoms = 0
        self.Rmin = None
        self.Rmax = None
        self.thetamin = None
        self.thetamax = None
        self.DWmin = None
        self.DWmax = None
        self.Tstart = 100
        self.Tend = 1
        self.cool = 0.9
        self.maxrun = 1000
        self.MC = 0.5
        self.RMS_count_max = 15
        self.factor = 5E5
        
#####################################################################################    
class RasdAna:
    """
    Rasd Data object without image data to be used in Modell Refinement and Fourier component extraction 
    """
    def __init__(self):
        self.E = num.array([], float) #Energies to be changed by e0shift
        self.Eorig = num.array([], float) #Energies to remain unchanged
        self.F = num.array([], float) # measured Structure Factors
        self.Ferr = num.array([], float) 
        self.f1 = num.array([], float) #f1 to be shifted by e0shift
        self.f2 = num.array([], float) #f2 to be shifted by e0shift
        self.file = '' #*.rsd file

        self.Q = num.ndarray((3), float) #q in fractional reziprocal coordinates
        self.mod_Q = float #norm Q in Ang**-1
        self.re_FNR = float #real part of non resonant Structure Factor
        self.im_FNR = float #imaginary part of non resonant Structure Factor

        self.re_FR = num.array([], float) #real part of resonant Structure Factor
        self.im_FR = num.array([], float) #imaginary part of resonant Structure Factor
        self.F_calc = num.array([], float) #modulus of calculated structure factor
        self.a = 0.1
        self.b = 0.1
        self.E0 = float
        self.e0shift = 0.
        self.Emin = float
        self.Emax = float
        self.norm = num.array([], float)
        self.delta_F = num.array([], float)
        self.AR = 0.2
        self.PR = 0.5
        self.use_in_Fourier = False
        self.use_in_Refine = False
        self.RMS = 0.
        self.ndata = 0

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
                title(str(self.file)+'  error: '+str(self.RMS))
                errorbar(self.E,self.F/self.norm, self.Ferr/self.norm,fmt = 'o')
                plot(self.E, self.F_calc/self.norm)
            else:
                figure(fig)
                clf()
                title(str(self.file)+'  error: '+str(self.RMS))
                errorbar(self.E,self.F, self.Ferr,fmt = 'o')
                plot(self.E, self.F_calc)

        
###################################################################################################
#########################  Read *.rsd files into RasdList  ################################################
def read_RSD(cell, bulk, surface, database, rasddata, f1f2, E0):
    """    
    read in data in RasdList of RasdAna objects for fitting
   """
    allrasd = RasdList()
    g_inv = calc_g_inv(cell)
    allrasd.ndata = 0
    for dat in rasddata:
        Q = None
        f = file(dat, 'r')
        data = f.readlines()
        f.close()
        Rasd = RasdAna()
        for i in range(len(data)):
            if '#' not in data[i]:
                tmp = str.rsplit(data[i])
                if Q == None:
                    Q = num.array([tmp[1],tmp[2],tmp[3]], float)
                    Rasd.Q = Q
                    Rasd.mod_Q = num.dot(num.dot(Q,g_inv),Q)**0.5
                    zeta = Q[2] + cell[6]*Q[0] + cell[7]*Q[1]
                    factor = 4* num.sin(num.pi*zeta)**2
    
                    re_ctr = (1 - num.cos(2*num.pi*zeta))/factor
                    im_ctr = -num.sin(2*num.pi*zeta)/factor
                    re_bulk , im_bulk = calc_Fuc(Q,bulk,g_inv,database)
                    re_surf, im_surf = calc_Fsurf(Q,surface,g_inv,database)

                    re_bc = re_ctr*re_bulk - im_ctr*im_bulk
                    im_bc = re_bulk*im_ctr + re_ctr*im_bulk
    
                    Rasd.re_FNR = re_bc + re_surf
                    Rasd.im_FNR = im_bc + im_surf
                Rasd.E = num.append(Rasd.E, int(round(float(tmp[0]))))
                Rasd.F = num.append(Rasd.F, float(tmp[4])**2)
                Rasd.Ferr = num.append(Rasd.Ferr, float(tmp[5])**2)
        Rasd.E0 = E0
        Rasd.Eorig = Rasd.E
        Rasd.ndata = len(Rasd.E)
        Rasd.file = dat
        allrasd.ndata = allrasd.ndata + Rasd.ndata
        
        i=0
        Rasd.f1 = num.ndarray((Rasd.ndata), float)
        Rasd.f2 = num.ndarray((Rasd.ndata), float)
        for i in range(Rasd.ndata):
            j = 0
            for j in range(len(f1f2[0])):
                if Rasd.E[i] == f1f2[0][j]:
                    Rasd.f1[i] = f1f2[1][j]
                    Rasd.f2[i] = f1f2[2][j]
        
        allrasd.list.append(Rasd)
    allrasd.dims = len(allrasd.list)
    allrasd.E = f1f2[0]
    allrasd.f1 = f1f2[1]
    allrasd.f2 = f1f2[2]
    allrasd.E0 = E0
    allrasd.cell = cell
    allrasd.g_inv = g_inv

    return allrasd

########################## Fit Fourier Components ####################################################    
    
def RASD_Fourier(allrasd, pnt):
    Rasd = allrasd.list[pnt]
    for i in range(Rasd.ndata):
        j = 0
        for j in range(len(allrasd.E)):
            if Rasd.E[i] == allrasd.E[j]:
                Rasd.f1[i] = allrasd.f1[j]
                Rasd.f2[i] = allrasd.f2[j]

    vec = num.array([Rasd.a ,Rasd.b , Rasd.AR, Rasd.PR], float)

    def Wiggle(vec, Rasd):
        a , b, AR, PR = vec
        re_Fq = AR * num.cos(2*num.pi*PR)
        im_Fq = AR * num.sin(2*num.pi*PR)
        re_FR = Rasd.f1 * re_Fq - Rasd.f2 * im_Fq
        im_FR = Rasd.f1 * im_Fq + Rasd.f2 * re_Fq
        return Rasd.F - ((a + b*(Rasd.E-Rasd.E0))* ((Rasd.re_FNR + re_FR)**2 + (Rasd.im_FNR + im_FR)**2))

    def Jacobi(vec, Rasd):
        J = num.ndarray((4, Rasd.ndata), float)
        a , b, AR, PR = vec
        re_Fq = AR * num.cos(2*num.pi*PR)
        im_Fq = AR * num.sin(2*num.pi*PR)
        re_FR = Rasd.f1 * re_Fq - Rasd.f2 * im_Fq
        im_FR = Rasd.f1 * im_Fq + Rasd.f2 * re_Fq

        J[0,:] = -(Rasd.re_FNR + re_FR)**2 - (Rasd.im_FNR + im_FR)**2

        J[1,:] = -(Rasd.E-Rasd.E0) * ((Rasd.re_FNR + re_FR)**2 + (Rasd.im_FNR + im_FR)**2)
        J[2,:] = -2 *(a+b*(Rasd.E-Rasd.E0))* ((Rasd.f1*num.cos(2*num.pi*PR)-Rasd.f2*num.sin(2*num.pi*PR))*(Rasd.re_FNR + re_FR)+\
                                          (Rasd.f1*num.sin(2*num.pi*PR)+Rasd.f2*num.cos(2*num.pi*PR))*(Rasd.im_FNR + im_FR))
                 
        J[3,:] = -4 *num.pi* (a+b*(Rasd.E-Rasd.E0))* ((-num.sin(2*num.pi*PR)*AR*Rasd.f1- num.cos(2*num.pi*PR)*Rasd.f2*AR)*(Rasd.re_FNR + re_FR) +\
                                                  ( num.cos(2*num.pi*PR)*AR*Rasd.f1- num.sin(2*num.pi*PR)*Rasd.f2*AR)*(Rasd.im_FNR + im_FR))
        return J
    
    result = leastsq(Wiggle, vec,args = (Rasd), Dfun = Jacobi, col_deriv = 1)
    Rasd.a , Rasd.b, Rasd.AR, Rasd.PR = result[0]
                          
    re_Fq = Rasd.AR * num.cos(2*num.pi*Rasd.PR)
    im_Fq = Rasd.AR * num.sin(2*num.pi*Rasd.PR)
    Rasd.re_FR = Rasd.f1 * re_Fq - Rasd.f2 * im_Fq
    Rasd.im_FR = Rasd.f1 * im_Fq + Rasd.f2 * re_Fq

    Rasd.F_calc = (Rasd.a + Rasd.b*(Rasd.E-Rasd.E0))* ((Rasd.re_FNR + Rasd.re_FR)**2 + (Rasd.im_FNR + Rasd.im_FR)**2)   
    Rasd.delta_F = Rasd.F - Rasd.F_calc
    Rasd.norm = (Rasd.a+Rasd.b*(Rasd.E-Rasd.E0))*(Rasd.re_FNR**2+Rasd.im_FNR**2)
    Rasd.RMS = num.sum((Rasd.delta_F/Rasd.norm)**2)/Rasd.ndata
    
    return Rasd

##################  Refinement of Atom coordinates, occupancies and DW- Factors  ###################################################################

def Wiggle2(vec, Rasd):
    a, b = vec
    delta_F = Rasd.F - ((a+b*(Rasd.E-Rasd.E0))*((Rasd.re_FNR + Rasd.re_FR)**2 + (Rasd.im_FNR + Rasd.im_FR)**2))
    return delta_F

def Jacobi2(vec, Rasd):
    J = num.ndarray((2,len(Rasd.E)),float)
    J[0][0:] = -(Rasd.re_FNR + Rasd.re_FR)**2 - (Rasd.im_FNR + Rasd.im_FR)**2
    J[1][0:] = -(Rasd.E-Rasd.E0) * ((Rasd.re_FNR + Rasd.re_FR)**2 + (Rasd.im_FNR + Rasd.im_FR)**2)
    return J

### main function that calculates the difference between measured and calculated Fs for a set of: R, theta, and DW ###
def Rasd_difference(allrasd, R, theta, DW):
    allrasd.RMS = 0
    allrasd.ndata = 0
    U = num.zeros((3,3),float)
    U[0][0] = DW[1]
    U[1][1] = DW[2]
    U[2][2] = DW[3]
    U[0][1] = DW[4]
    U[1][0] = DW[4]
    U[0][2] = DW[5]
    U[2][0] = DW[5]
    U[1][2] = DW[6]
    U[2][1] = DW[6]
                  
    for Rasd in allrasd.reflist:
        Q = num.zeros((3),float)
        Q[0] = allrasd.g_inv[0][0]**0.5 * Rasd.Q[0]
        Q[1] = allrasd.g_inv[1][1]**0.5 * Rasd.Q[1]
        Q[2] = allrasd.g_inv[2][2]**0.5 * Rasd.Q[2]
        
        j = 0
        re_Fq = 0
        im_Fq = 0
        for j in range(len(theta)):
            re_Fq = re_Fq + (theta[j] * num.cos(2*num.pi*(Rasd.Q[0]*allrasd.g_inv[0][0]**0.5*R[j][0]+Rasd.Q[1]*allrasd.g_inv[1][1]**0.5*R[j][1]+\
                                                          Rasd.Q[2]*allrasd.g_inv[2][2]**0.5*R[j][2]))* num.exp(-0.5 * num.dot(num.dot(Q,U),Q)))
            im_Fq = im_Fq + (theta[j] * num.sin(2*num.pi*(Rasd.Q[0]*allrasd.g_inv[0][0]**0.5*R[j][0]+Rasd.Q[1]*allrasd.g_inv[1][1]**0.5*R[j][1]+\
                                                          Rasd.Q[2]*allrasd.g_inv[2][2]**0.5*R[j][2]))* num.exp(-0.5 * num.dot(num.dot(Q,U),Q)))
        Rasd.re_FR = Rasd.f1 * re_Fq - Rasd.f2 * im_Fq
        Rasd.im_FR = Rasd.f1 * im_Fq + Rasd.f2 * re_Fq

        vec = [Rasd.a, Rasd.b]
        result = leastsq(Wiggle2, vec, args=(Rasd), Dfun = Jacobi2, col_deriv = 1)
        vec = result[0]
        Rasd.a, Rasd.b = vec

        Rasd.F_calc = (Rasd.a+Rasd.b*(Rasd.E - Rasd.E0))*((Rasd.re_FNR + Rasd.re_FR)**2 + (Rasd.im_FNR + Rasd.im_FR)**2) 
        Rasd.delta_F = Rasd.F - Rasd.F_calc
        Rasd.norm = (Rasd.a+ Rasd.b* (Rasd.E - Rasd.E0))*(Rasd.re_FNR**2 + Rasd.im_FNR**2)
        Rasd.RMS = (num.sum((Rasd.delta_F/Rasd.norm)**2)/Rasd.ndata)**0.5
        
        allrasd.RMS = allrasd.RMS + Rasd.RMS**2 * Rasd.ndata
        allrasd.ndata = allrasd.ndata + Rasd.ndata
        
    allrasd.RMS = (allrasd.RMS/allrasd.ndata)**0.5

    return allrasd

def simulated_annealing(allrasd):
    # initialize values
    Rmin = allrasd.Rmin
    Rmax = allrasd.Rmax
    thetamin = allrasd.thetamin
    thetamax = allrasd.thetamax
    DWmin = allrasd.DWmin 
    DWmax = allrasd.DWmax
    Tstart = allrasd.Tstart
    Tend = allrasd.Tend
    f = allrasd.cool
    maxrun = allrasd.maxrun
    MC = allrasd.MC /100
    RMS_count_max = allrasd.RMS_count_max
    factor = allrasd.factor
    
    RP = num.ndarray((len(thetamin),3),float)
    thetaP = num.ndarray((len(thetamin)),float)
    DWP = num.ndarray((len(thetamin)),float)
    for i in range(len(thetamin)):
        for j in range(3):
           RP[i][j] = random.uniform(Rmin[i][j],Rmax[i][j]) 
        thetaP[i] = random.uniform(thetamin[i], thetamax[i])
        for j in range(6):
            DWP[i][j] = random.uniform(DWmin[i][j], DWmax[i][j])
    
    RMS = Rasd_difference(allrasd, RP, thetaP, DWP).RMS
    RMS_test = RMS
    R_track =num.array([RMS],float)
    param_track = [[RP, thetaP, DWP]]

    #counters
    Random = 0
    better = 0
    rejected = 0
    RMS_count = 0
    while Tstart > Tend and RMS_count < RMS_count_max:
        z = 0
        while (z < maxrun):
            ## random perturbation of parameters
            i=0
            R_tmp = num.ndarray((len(DWmin),3),float)
            theta_tmp = num.ndarray((len(DWmin)),float)
            DW_tmp = num.ndarray((len(DWmin)),float)
            for i in range(allrasd.natoms):
                for j in range(3):
                    R_tmp[i][j] = RP[i][j] * random.uniform( 1-(MC * Tstart), 1+(MC * Tstart))
                    if R_tmp[i][j] < Rmin[i][j]: R_tmp[i][j] = Rmin[i][j]
                    elif R_tmp[i][j] > Rmax[i][j]: R_tmp[i][j] = Rmax[i][j]
                theta_tmp[i] = thetaP[i] * random.uniform( 1-(MC * Tstart), 1+(MC * Tstart))
                if theta_tmp[i] < thetamin[i]: theta_tmp[i] = thetamin[i]
                elif theta_tmp[i] > thetamax[i]: theta_tmp[i] = thetamax[i]
                for j in range(6):
                    DW_tmp[i][j] = DWP[i][j] * random.uniform( 1-(MC * Tstart), 1+(MC * Tstart))
                    if DW_tmp[i][j] < DWmin[i][j]: DW_tmp[i][j] = DWmin[i][j]
                    elif DW_tmp[i][j] > DWmax[i][j]: DW_tmp[i][j] = DWmax[i][j]

            RMS_tmp = Rasd_difference(allrasd, R_tmp, theta_tmp, DW_tmp).RMS

            dR = (RMS - RMS_tmp) * factor

            if dR > 0:
                RP = R_tmp
                thetaP = theta_tmp
                DWP = DW_tmp
                RMS = RMS_tmp
                z = z+1
                better = better+1
                R_track = num.append(R_track,RMS)
                param_track.append([R_tmp,theta_tmp,DW_tmp])

            elif dR <= 0:
                Boltz = exp(dR / Tstart)
                Rand = random.uniform(0,1)
                if(Boltz > Rand):
                    RP = R_tmp
                    thetaP = theta_tmp
                    DWP = DW_tmp
                    RMS = RMS_tmp
                    z = maxrun
                    Random = Random + 1
                elif Boltz <= Rand:
                    z = z+1
                    rejected = rejected+1
            del(R_tmp, theta_tmp, DW_tmp, RMS_tmp)
        print 'Temperature: '+str(Tstart)
        print 'RMS error: '+str(RMS)
       
        if RMS == RMS_test:
            RMS_count = RMS_count+1
        else:
            RMS_test = RMS
            RMS_count = 0

        Tstart = Tstart * f

    print '****************************'
    print 'number of cycles: '+str(Random+ better+ rejected)
    print 'random: ' + str(Random)
    print 'better: ' + str(better)
    print 'rejected: '+ str(rejected)

    mini = num.where(R_track == R_track.min())
    R, theta, DW = param_track[int(mini[0][0])]
    allrasd = Rasd_difference(allrasd, R, theta, DW)
    R_track = num.append(R_track, allrasd.RMS)

    print '####################################################\n'
    print 'the best fit RMS error is: '+str(allrasd.RMS)+'\n'
    i=0
    for i in range(len(DW)):
        print '*************************************'
        print 'atom '+str(i+1)+': \n'
        print 'x, y, z: '+ str(R[i])+' \n'
        print 'occupancy = '+str(theta[i])+' \n'
        print 'DW(iso) = '+str(DW[i])+' \n \n'

    figure(4)
    clf()
    title('Overall RMS error: '+str(allrasd.RMS))
    for i in range(len(allrasd.reflist)):
        errorbar(allrasd.reflist[i].E,(allrasd.reflist[i].F/allrasd.reflist[i].norm)+i*0.5, (allrasd.reflist[i].Ferr/allrasd.reflist[i].norm),fmt = 'bo')
        plot(allrasd.reflist[i].E, (allrasd.reflist[i].F_calc/allrasd.reflist[i].norm)+i*0.5, 'g-' )
        text(allrasd.reflist[i].E0,i*0.5+1,(allrasd.reflist[i].file +', RMS err.= ' +str(round(allrasd.reflist[i].RMS,4))))
    figure(5)
    clf()
    title('Development of RMS during fit')
    plot(range(len(R_track)),R_track,'ro')

#################  Fourier Synthese  ############################################
def calc_rho(Fourier, r, ZR, V, cell):
    rho = 0
    ginv = calc_g_inv(cell)
    for F_comp in Fourier:
        Q = num.array([F_comp[0]*2*num.pi*ginv[0][0]**0.5,F_comp[1]*2*num.pi*ginv[1][1]**0.5,F_comp[2]*2*num.pi*ginv[2][2]**0.5], float)
        rho = rho + (F_comp[3] * num.cos(2*num.pi*F_comp[4] - num.dot(r,Q)))
    rho = rho * ZR / V

    return rho

def Fourier_synthesis(Fourier, cell, ZR = 1, xf = 1, yf = 1, zf = 1, an = 11, bn = 11, cn = 61, Plusminus = False):
    
    Rho = num.ndarray((an,bn,cn),float)
    sampx = cell[0]* xf
    sampy = cell[1]* yf
    sampz = cell[2]* zf
    g_inv = calc_g_inv([sampx,sampy,sampz,cell[3],cell[4],cell[5]])
    g = num.linalg.inv(g_inv)
    V = num.linalg.det(g)
    
    i = 0
    for i in range(Rho.shape[0]):
        x = sampx /(Rho.shape[0]-1) * float(i)
        j = 0
        for j in range(Rho.shape[1]):
            y = sampy /(Rho.shape[1]-1) * float(j)
            k = 0
            for k in range(Rho.shape[2]):
                if Plusminus:
                    z = sampz /(Rho.shape[2]-1) * (float(k)-((Rho.shape[2]-1)/2))
                else:
                    z = sampz /(Rho.shape[2]-1) * float(k)
                if k == 0: zmin = z
                if k == Rho.shape[2]-1: zmax = z
                R = num.array([x,y,z],float)
                Rho[i][j][k]= calc_rho(Fourier, R, ZR, V, cell)

    Rho_XY = num.sum(Rho, axis = 2)
    Rho_XZ = num.sum(Rho, axis = 1)
    Rho_YZ = num.sum(Rho, axis = 0)
    Rho_YZ = Rho_YZ.transpose()

    Max = num.max(Rho)
    xaxis = num.arange(Rho.shape[2]) /(float(Rho.shape[2])-1) *(zmax-zmin) +zmin
    figure(2, figsize = (12,9))
    subplot(221)
    plot(xaxis, num.sum(Rho_XZ,axis=0),'-')
    ylabel('rho_e-(Z)')
    xlabel('Z')
    subplot(222)
    title('sum over X')
    im = imshow(Rho_YZ, interpolation='bilinear', cmap=cm.hot, origin='lower', extent = [0,sampy,zmin,zmax])
    xlabel('Y')
    ylabel('Z')
    subplot(223)
    title('sum over Y')
    im = imshow(Rho_XZ, interpolation='bilinear', cmap=cm.hot, origin='lower', extent = [zmin,zmax,0,sampx])
    xlabel('Z')
    ylabel('X')
    subplot(224)
    title('sum over Z')
    imshow(Rho_XY, interpolation='bilinear', cmap=cm.hot, origin='lower', extent = [0,sampy,0,sampx])
    xlabel('Y')
    ylabel('X')

    i = 0
    figure(3, figsize = (19,6))
    for i in range(9):
        pl = '19'+str(i+1)
        subplot(pl)
        a = int(round(float(Rho.shape[0]-1)/8 *i))
        if a > Rho.shape[0]-1: a = Rho.shape[0]-1
        imshow(Rho[a].transpose(), interpolation='bilinear', cmap=cm.hot, origin='lower', vmin = 0, vmax = Max, extent = [0,sampy,zmin,zmax] )
        title('X = '+str(round(sampx /8 * i,4)))
        xlabel('Y')
        if i == 0: ylabel('Z')

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
    #[0] - E
    #[1] - f1 theo
    #[2] - f2 theo
    #[3] - diff f1theo f1exp 
    #[4] - f1 exp
    #[5] - f2 exp
    
    e0 = e0 + e0shift
    #read Cromer-Liberman calculated f1f2 (1eV grid) from HEPHAESTUS
    f = file(datafile, 'r')
    data = f.readlines()
    f.close()
    f1f2 = num.ndarray((0,6),float)
    for i in range(len(data)):
        if '#' not in data[i]:
            tmp = str.rsplit(data[i])
            if float(tmp[0])< e0:
                f1f2 = num.append(f1f2, [[int(round(float(tmp[0]))),float(tmp[1]),float(tmp[2]),0.,0.,0.]], axis = 0)
            else:
                f1f2 = num.append(f1f2, [[int(round(float(tmp[0]))),float(tmp[1]),float(tmp[2]),0.,0.,1.]], axis = 0)
    f1f2 = f1f2.transpose()
    f1f2[0] = f1f2[0] + e0shift

    #read experimental f2 (normalized XANES)
    f = file(expfile, 'r')
    data = f.readlines()
    f.close()
    f2exp = num.ndarray((0,2),float)
    for i in range(len(data)):
        if '#' not in data[i]:
            tmp = str.rsplit(data[i])
            f2exp = num.append(f2exp, [[int(round(float(tmp[0]))),float(tmp[1])]], axis = 0)
    f2exp = f2exp.transpose()

    #associate experimental values to calculated values
    i=0
    for i in range(len(f1f2[0])):
        j=0
        for j in range(len(f2exp[0])):
            if f1f2[0][i]== f2exp[0][j]:
                f1f2[5][i] = f2exp[1][j]


    lower = 0
    upper = 0
    i = 0
    for i in range(n):
        lower = lower + (f1f2[2][i] - f1f2[5][i])
        upper = upper + (f1f2[2][len(f1f2[0])-i-1]- f1f2[5][len(f1f2[0])-i-1]+ 1)
    
    lower = lower /n
    upper = upper /n
    f1f2[5] = f1f2[5]*(upper-lower) + lower

    #calculate f1exp from f2exp (Kramers-Kronig) (see Ohta 1988/ Cross 1998)
    for i in range(num.shape(f1f2)[1]):
        sum = 0
        if divmod(float(i),2)[1] == 0:
            j = 1
            for j in range(1, len(f1f2[0]),2):
                sum = sum + (f1f2[5][j]-f1f2[2][j])/(f1f2[0][j] - f1f2[0][i])+(f1f2[5][j]-f1f2[2][j])/(f1f2[0][j] + f1f2[0][i])
        else:
            j = 0
            for j in range(0, len(f1f2[0]),2):
                sum = sum + (f1f2[5][j]-f1f2[2][j])/(f1f2[0][j] - f1f2[0][i])+(f1f2[5][j]-f1f2[2][j])/(f1f2[0][j] + f1f2[0][i])
        f1f2[3][i] = (sum * 2 / num.pi)


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


