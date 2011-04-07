import numpy as Num
from pylab import *
import random
from atomic import f0data as database
from bv_params import bv_params

###################################### calculations ####################################################
def calc_g_inv(cell):
    g = Num.ndarray((3,3),float)
    g[0][0] = cell[0]**2
    g[0][1] = cell[0]*cell[1]*Num.cos(Num.radians(cell[5]))
    g[0][2] = cell[0]*cell[2]*Num.cos(Num.radians(cell[4]))
    g[1][0] = cell[1]*cell[0]*Num.cos(Num.radians(cell[5]))
    g[1][1] = cell[1]**2
    g[1][2] = cell[1]*cell[2]*Num.cos(Num.radians(cell[3]))
    g[2][0] = cell[2]*cell[0]*Num.cos(Num.radians(cell[4]))
    g[2][1] = cell[2]*cell[1]*Num.cos(Num.radians(cell[3]))
    g[2][2] = cell[2]**2

    g_inv = Num.linalg.inv(g)
    return g_inv

def calc_g(cell):
    g = Num.ndarray((3,3),float)
    g[0][0] = cell[0]**2
    g[0][1] = cell[0]*cell[1]*Num.cos(Num.radians(cell[5]))
    g[0][2] = cell[0]*cell[2]*Num.cos(Num.radians(cell[4]))
    g[1][0] = cell[1]*cell[0]*Num.cos(Num.radians(cell[5]))
    g[1][1] = cell[1]**2
    g[1][2] = cell[1]*cell[2]*Num.cos(Num.radians(cell[3]))
    g[2][0] = cell[2]*cell[0]*Num.cos(Num.radians(cell[4]))
    g[2][1] = cell[2]*cell[1]*Num.cos(Num.radians(cell[3]))
    g[2][2] = cell[2]**2
    return g

#####################################################################################################
class rigid_body:
    def __init__(self):
        self.label = ''
        self.atoms = []
        self.angles = []
###########################################        
def rigid_body_rotation(atoms, theta, phi, chi, cell):
    Center  = Num.array([atoms[0][1],atoms[0][2],atoms[0][3]],float)
      
    theta = Num.radians(theta)
    phi = Num.radians(phi)
    chi = Num.radians(chi)
    R_theta = Num.array([[Num.cos(theta),Num.sin(theta),0],[-Num.sin(theta),Num.cos(theta),0],[0,0,1]],float)
    R_phi   = Num.array([[Num.cos(phi),0,Num.sin(phi)],[0,1,0],[-Num.sin(phi),0,Num.cos(phi)]],float)
    R_chi   = Num.array([[1,0,0],[0,Num.cos(chi),Num.sin(chi)],[0,-Num.sin(chi),Num.cos(chi)]],float)
    R  = Num.dot(R_theta,Num.dot(R_phi,R_chi))
    P = Num.array([[cell[0],0,0],[0,cell[1],0],[0,0,cell[2]]],float)
    R = Num.dot(Num.linalg.inv(P),Num.dot(R,P))
    new_atoms = []
    for atom in atoms:
        coords = Num.array([atom[1],atom[2],atom[3]],float) - Center
        coords = Num.dot(R, coords) + Center
        new_atom = [atom[0], coords[0], coords[1], coords[2], atom[4], atom[5], atom[6], atom[7], atom[8], atom[9], atom[10]]
        new_atoms.append(new_atom)
    return new_atoms
############################################
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
##########################################################################################################################
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
    ######################################
    def calc_BVS(self, surface):
        center = Num.array([surface[self.center][1]+self.centerxoffset,surface[self.center][2]+self.centeryoffset,surface[self.center][3]],float)
        BVS = 0
        distances = []
        for i in range(len(self.neighbors)):
            neighbor = Num.array([surface[self.neighbors[i]][1]+self.neighborsxoffset[i],surface[self.neighbors[i]][2]+self.neighborsyoffset[i],surface[self.neighbors[i]][3]],float)
            vector = neighbor - center
            dist = (Num.dot(vector,Num.dot(self.g,vector)))**0.5
            distances.append(dist)
            BVS = BVS + Num.exp((self.r0s[i]-dist)/self.bs[i])
        return BVS, distances
#########################################################################################
def BV_impact(BVclusters, surface):
    impact = 0
    for i in BVclusters:
        BVS, dist = i.calc_BVS(surface)
        eqval = (float(i.eqval) **2)**0.5
        BV_offset = ((BVS - eqval)**2)**0.5 / eqval
        impact = impact + (i.ip[0] * BV_offset)**i.ip[1]
    return impact
##########################################################################################################################    
def param_unfold(param, param_use, surface, use_bulk_water):
    if use_bulk_water:
        zwater = param['zwater'][0]
        sig_water = param['sig_water'][0]
    else:
        zwater = 0
        sig_water = 1
    
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
            atom[3] = surface[i][3]+ param_use[i][4]* param[param_use[i][5]][0] -zwater
        else: atom[3] = surface[i][3] -zwater
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
            atom[10] = param_use[i][18]* param[param_use[i][19]][0]
        else: atom[10] = surface[i][10]
        surface_new.append(atom)

    return zwater, sig_water, Scale, specScale, beta, surface_new

##########################################################################################################################
def write_bulk(bulk,param, use_bulk_water, filename = 'shifted_bulk.bul'):
    if use_bulk_water:
        zwater = param['zwater'][0]
    else:
        zwater = 0
    bulk = shift_bulk(zwater, bulk)
    f = file(filename, 'w')
    for atom in bulk:
        line = "%5s %6.5f %6.5f %6.5f %6.5f \n" % (atom[0],atom[1],atom[2],atom[3],atom[4])
        f.write(line)
##########################################################################################################################
def write_surface(cell, surface,param,param_use, rigid_bodies, use_bulk_water, filename = 'surface.sur'):
    zwater, sig_water, Scale,specScale, beta, surface = param_unfold(param,param_use, surface, use_bulk_water)
    surface = RB_update(rigid_bodies, surface, param, cell)
    f = file(filename, 'w')
    for atom in surface:
        line = "%5s %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f\n" % (atom[0],atom[1],atom[2],atom[3],\
                                                                                      atom[4],atom[5],atom[6],atom[7],atom[8],atom[9],atom[10])
        f.write(line)    
##########################################################################################################################
def write_cif(cell4, surface4,param4,param_use, rigid_bodies, use_bulk_water, filename = 'surface.cif'):
    zwater, sig_water, Scale,specScale, beta, surface4 = param_unfold(param4,param_use, surface4, use_bulk_water)
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
        f.write(str(surface4[i][0])+str(i+1)+'  '+str(surface4[i][4])+'  '+str(surface4[i][5])+'  '+str(surface4[i][6])+'  '+str(surface4[i][7])+'  '+str(surface4[i][8])+'  '+str(surface4[i][9])+'\n')
        f.write(str(surface4[i][0])+str(len(surface4)+i+1)+'  '+str(surface4[i][4])+'  '+str(surface4[i][5])+'  '+str(surface4[i][6])+'  '+str(surface4[i][7])+'  '+str(surface4[i][8])+'  '+str(surface4[i][9])+'\n')
        f.write(str(surface4[i][0])+str(2*len(surface4)+i+1)+'  '+str(surface4[i][4])+'  '+str(surface4[i][5])+'  '+str(surface4[i][6])+'  '+str(surface4[i][7])+'  '+str(surface4[i][8])+'  '+str(surface4[i][9])+'\n')
        f.write(str(surface4[i][0])+str(3*len(surface4)+i+1)+'  '+str(surface4[i][4])+'  '+str(surface4[i][5])+'  '+str(surface4[i][6])+'  '+str(surface4[i][7])+'  '+str(surface4[i][8])+'  '+str(surface4[i][9])+'\n')
    f.close()

############################### Fitting Rod #########################################################################
class Fitting_Rod:
    def __init__(self):
        self.H = float
        self.K = float
        self.L = Num.array([],float)
        self.F = Num.array([],float)
        self.Ferr = Num.array([],float)
        self.Lb = Num.array([],float)
        self.Db = float
        
        self.bulk = Num.array([],float)
        self.surf = Num.array([],float)
        self.water = Num.array([],float)
        self.rough = Num.array([],float)
        self.Fcalc = Num.array([],float)
        self.difference = Num.array([],float)
    #####################################################
    def calc_Fuc(self,hkl,bulk,g_inv,database):
        a = 0
        b = 0
        for i in range(shape(bulk)[0]):
            f_par = database[str.lower(bulk[i][0])]
            q = (Num.dot(Num.dot(hkl,g_inv),hkl))**0.5 
            f = (f_par[0]*Num.exp(-(q/4/Num.pi)**2*f_par[1]) + f_par[2]*Num.exp(-(q/4/Num.pi)**2*f_par[3]) +\
                f_par[4]*Num.exp(-(q/4/Num.pi)**2*f_par[5]) + f_par[6]*Num.exp(-(q/4/Num.pi)**2*f_par[7]) + f_par[8])*\
                Num.exp(-(q/2)**2 *bulk[i][4])
            a = a + (f * Num.cos(2*Num.pi*(hkl[0]*bulk[i][1] + hkl[1]*bulk[i][2] + hkl[2]*bulk[i][3])))
            b = b + (f * Num.sin(2*Num.pi*(hkl[0]*bulk[i][1] + hkl[1]*bulk[i][2] + hkl[2]*bulk[i][3])))
        return a, b

    def calc_Fsurf(self,hkl,surface,g_inv,database):
        a = 0
        b = 0
        for i in range(shape(surface)[0]):
            f_par = database[str.lower(surface[i][0])]
            U = Num.array([[surface[i][4],surface[i][7],surface[i][8]],[surface[i][7],surface[i][5],surface[i][9]],[surface[i][8],surface[i][9],surface[i][6]]])
            q = (Num.dot(Num.dot(hkl,g_inv),hkl))**0.5
            q_Ang = [hkl[0]*g_inv[0][0]**0.5, hkl[1]*g_inv[1][1]**0.5, hkl[2]*g_inv[2][2]**0.5] 
            f = (f_par[0]*Num.exp(-(q/4/Num.pi)**2*f_par[1]) + f_par[2]*Num.exp(-(q/4/Num.pi)**2*f_par[3]) +\
                f_par[4]*Num.exp(-(q/4/Num.pi)**2*f_par[5]) + f_par[6]*Num.exp(-(q/4/Num.pi)**2*f_par[7]) + f_par[8])*\
                Num.exp(-2* Num.pi**2*(Num.dot(q_Ang,Num.dot(U,q_Ang)))) * surface[i][10]
            a = a + (f * Num.cos(2*Num.pi*(hkl[0]*surface[i][1] + hkl[1]*surface[i][2] + hkl[2]*surface[i][3])))
            b = b + (f * Num.sin(2*Num.pi*(hkl[0]*surface[i][1] + hkl[1]*surface[i][2] + hkl[2]*surface[i][3])))
        return a, b


    def calc_Fwater(self,hkl, sig, g_inv, database, cell):
        f_par = database['o2-.']
        q = hkl[2]* g_inv[2][2]**0.5
        Auc = cell[0]* Num.sin(Num.radians(cell[5]))* cell[1]
        f = (f_par[0]*Num.exp(-(q/4/Num.pi)**2*f_par[1]) + f_par[2]*Num.exp(-(q/4/Num.pi)**2*f_par[3]) +\
                f_par[4]*Num.exp(-(q/4/Num.pi)**2*f_par[5]) + f_par[6]*Num.exp(-(q/4/Num.pi)**2*f_par[7]) + f_par[8])*\
                Num.exp(-(q/2)**2*sig)
        b = f * 0.033456 * Auc / (q*2*Num.pi)
        return b

    def calcF(self,sig_water,Scale,specScale,beta,cell,bulk,surface,g_inv,NLayers,database, use_bulk_water, RMS_flag):
        self.bulk = Num.ndarray((len(self.L)),float)
        self.surf = Num.ndarray((len(self.L)),float)
        self.water = Num.ndarray((len(self.L)),float)
        self.rough = Num.ndarray((len(self.L)),float)
        self.Fcalc = Num.ndarray((len(self.L)),float)
        self.difference = Num.array((len(self.L)),float)
        for i in range(len(self.L)):
            hkl = [self.H,self.K,self.L[i]]
            zeta = self.L[i]+ self.H*cell[6]+ self.K*cell[7]

            factor = 4* Num.sin(Num.pi*zeta)**2
            re_ctr = (1 - Num.cos(2*Num.pi*zeta))/factor
            im_ctr = -Num.sin(2*Num.pi*zeta)/factor

            re_bulk , im_bulk = self.calc_Fuc(hkl,bulk,g_inv,database)
            re_surf, im_surf = self.calc_Fsurf(hkl,surface,g_inv,database)
            
            re_bc = re_ctr*re_bulk - im_ctr*im_bulk
            im_bc = re_bulk*im_ctr + re_ctr*im_bulk

            if self.L[i] > 0:
                n = self.Lb[i] + round(self.L[i]/self.Db) * self.Db
            else:
                n = - self.Lb[i] + round(self.L[i]/self.Db) * self.Db
            self.rough[i] = (1-beta)/((1-beta)**2 + 4*beta*Num.sin(Num.pi*(self.L[i] - n)/NLayers)**2)**0.5
            
            if hkl[0] == 0.0 and hkl[1] == 0.0:
                if use_bulk_water:
                    im_water = self.calc_Fwater(hkl, sig_water, g_inv, database, cell)
                    self.water[i] = specScale * (im_water**2)**0.5
                else:
                    im_water = 0
                    self.water[i] = 0
                self.bulk[i] = specScale * (re_bc**2 + im_bc**2)**0.5
                self.Fcalc[i] = specScale * self.rough[i] * ((re_bc + re_surf)**2 + (im_bc + im_surf + im_water)**2)**0.5
                
                self.rough[i] = self.rough[i] * specScale
                self.surf[i] = (re_surf**2 + im_surf**2)**0.5 * specScale
            else:
                self.bulk[i] = Scale * (re_bc**2 + im_bc**2)**0.5
                self.Fcalc[i] = Scale * self.rough[i] * ((re_bc + re_surf)**2 + (im_bc + im_surf)**2)**0.5
                self.water[i] = 0
                self.rough[i] = self.rough[i] * Scale
                self.surf[i] = (re_surf**2 + im_surf**2)**0.5 * Scale

        if RMS_flag == 1:
            self.difference = ((Num.log(self.F) - Num.log(self.Fcalc))**2)**0.5
        elif RMS_flag == 2:
            self.difference = ((self.F - self.Fcalc)**2)**0.5
        elif RMS_flag == 3:
            self.difference = ((self.F - self.Fcalc)**2)**0.5 * self.Ferr
        
############################### reading datafiles ###################################################################
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
            tmp2 = [tmp[0], float(tmp[1]),float(tmp[2]),float(tmp[3]),float(tmp[4])]
            bulk.append(tmp2)
    return bulk, cell, Nlayers

def read_surface(surfacefile):
    surface=[]
    f = open(surfacefile, 'r')
    data = f.readlines()
    f.close()
    for i in data:
        tmp = str.rsplit(i)
        tmp2 = [tmp[0], float(tmp[1]),float(tmp[2]),float(tmp[3]),float(tmp[4]),float(tmp[5]),float(tmp[6]),float(tmp[7]),float(tmp[8]),float(tmp[9]),float(tmp[10])]
        surface.append(tmp2)
        
    parameter_usage=[]
    for i in data:
        tmp = str.rsplit(i)
        tmp2 = [float(tmp[11]),tmp[12],float(tmp[13]),tmp[14],float(tmp[15]),tmp[16],float(tmp[17]),tmp[18],float(tmp[19]),tmp[20],\
                float(tmp[21]),tmp[22],float(tmp[23]),tmp[24],float(tmp[25]),tmp[26],float(tmp[27]),tmp[28],float(tmp[29]),tmp[30]]
        parameter_usage.append(tmp2)

    return surface, parameter_usage

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
            tmp2.angles = [float(tmp[ct+2]),tmp[ct+3],float(tmp[ct+4]),tmp[ct+5],float(tmp[ct+6]),tmp[ct+7]]
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
            key = center_label + str(BVC.eqval) + neighbor_label + str(neighbor_valence)
            if key in bv_params.keys():
                BVC.r0s.append(bv_params[key][0])
                BVC.bs.append(bv_params[key][1])
            else:
                key = neighbor_label + str(neighbor_valence) + center_label + str(BVC.eqval)
                if key in bv_params.keys():
                    BVC.r0s.append(bv_params[key][0])
                    BVC.bs.append(bv_params[key][1])
                else: print ' No bond valence parameters found for: '+center_label+' '+str(BVC.eqval)+' and '+neighbor_label+' '+str(neighbor_valence)
            BVC.neighbors.append(int(tmp[8+j*5]))
            BVC.neighborsxoffset.append(int(tmp[9+j*5]))
            BVC.neighborsyoffset.append(int(tmp[10+j*5]))
        BVC.ip[0] = float(tmp[11+(n-2)*5])
        BVC.ip[1] = float(tmp[12+(n-2)*5])
        BVC.g = calc_g(cell)
        BVclusters.append(BVC)
    return BVclusters
                
#########################################################################################
def shift_bulk(zwater, bulk):
    bulk_new =[]
    for i in bulk:
        atom = [i[0],i[1],i[2],i[3]-zwater,i[4]]
        bulk_new.append(atom)
    return bulk_new
    
#########################################calculate rods##################################
def calc_CTRs(parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water,\
              use_BVC, BVclusters, RMS_flag):

    zwater, sig_water, Scale,specScale, beta, surface_new = param_unfold(parameter,param_usage, surface_tmp, use_bulk_water)
    bulk_new = shift_bulk(zwater, bulk_tmp)
    surface_new = RB_update(rigid_bodies, surface_new, parameter, cell) 
                                  
    for x in dat:
        x.calcF(sig_water,Scale,specScale,beta,cell,bulk_new,surface_new,g_inv,NLayers,database, use_bulk_water, RMS_flag)

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
            n = n + len(dat[i].L)*Rod_weight[i]
            avgF = avgF + Num.sum(dat[i].F)*Rod_weight[i]
        elif RMS_flag == 3:
            b = b + Num.sum(dat[i].Ferr)*Rod_weight[i]
            n = n + len(dat[i].L)*Rod_weight[i]
            avgF = avgF + Num.sum(dat[i].F)*Rod_weight[i]
            
    if RMS_flag == 1:
        RMS = (RMS / n)
        avgF = ((avgF / n)**2)**0.5
    elif RMS_flag == 2:
        RMS = (RMS / n)
        avgF = avgF / n
    elif RMS_flag == 3:
        RMS = (RMS / b)
        avgF = avgF / n

    RMS = (RMS/avgF)*100
        
    if use_BVC:
        impact = BV_impact(BVclusters, surface_new)
        RMS = RMS * (1 + impact)

    return dat, RMS
############################### Simulated Annealing #################################################################
def simulated_annealing01(dat, cell, NLayers, bulk, surface, database, Rod_weight, sim_an_params, parameter, param_usage, plot_RMS_track,\
                          rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag):

    Tstart,Tend,cool,maxrun,MC,factor,random_parameters = sim_an_params
    g_inv = calc_g_inv(cell)
    
    if random_parameters:
        for i in parameter.keys():
            if parameter[i][3]:
                parameter[i][0] = random.uniform(parameter[i][1], parameter[i][2])
      
    dat, RMS = calc_CTRs(parameter, param_usage, dat, cell, bulk, surface, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag)

    guess = (int(Num.log(float(Tend)/float(Tstart))/Num.log(cool))+1)*maxrun
    print 'approximated number of iterations: '+str( int(guess) )

    print 'R start = '+str(RMS)

    R_track =Num.array([RMS],float)
    param_track = [parameter]

    #counters
    Random = 0
    better = 0
    rejected = 0
    while Tstart > Tend:
        z = 0
        while (z < maxrun):
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

                
            dat, RMS_tmp = calc_CTRs(param_tmp,param_usage, dat, cell, bulk, surface, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag)

            dR = (RMS - RMS_tmp) * factor

            if dR > 0:
                for i in parameter.keys():
                    parameter[i][0] = param_tmp[i][0]
                z = z+1
                better = better+1
                R_track = Num.append(R_track,RMS_tmp)
                param_track.append(parameter)
                RMS = RMS_tmp
                print 'better    '+str(RMS_tmp)
            elif dR <= 0:
                Boltz = exp(dR / Tstart)
                Rand = random.uniform(0,1)
                if(Boltz > Rand):
                    for i in parameter.keys():
                        parameter[i][0] = param_tmp[i][0]
                    z = maxrun
                    Random = Random+1
                    R_track = Num.append(R_track,RMS_tmp)
                    param_track.append(parameter)
                    RMS = RMS_tmp
                    print 'random    '+str(RMS_tmp)
                elif Boltz <= Rand:
                    z = z+1
                    rejected = rejected+1
                    print 'rejected    '+str(RMS_tmp)
                    
        Tstart = Tstart * cool
        print '\n##################################################'
        print 'Temperature: '+str(Tstart/cool)
        print 'R = '+str(RMS)
        print '##################################################'+'\n'

    print '****************************'
    print 'Number of cycles: '+str(Random+ better+ rejected)
    print 'random: ' + str(Random)
    print 'better: ' + str(better)
    print 'rejected: '+ str(rejected)

    mini = Num.where(R_track == R_track.min())
    param_best = param_track[int(mini[0][0])]

    data_best, RMS_best = calc_CTRs(param_best,param_usage, dat, cell, bulk, surface, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag)
    R_track = Num.append(R_track, RMS_best)

    print '\n####################################################\n'
    print 'the best fit R = '+str(RMS_best)+'\n'
    print '*************************************\n'

    if plot_RMS_track:
        figure(3)
        clf()
        title('Development of R during fit')
        plot(range(len(R_track)),R_track,'ro')

    return data_best, param_best, RMS_best
############################### Simulated Annealing #################################################################
def simulated_annealing02(dat, cell, NLayers, bulk, surface, database, Rod_weight, sim_an_params, parameter, param_usage, plot_RMS_track,\
                          rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag):

    Tstart,Tend,cool,maxrun,MC,factor,random_parameters = sim_an_params
    g_inv = calc_g_inv(cell)
    
    if random_parameters:
        for i in parameter.keys():
            if parameter[i][3]:
                parameter[i][0] = random.uniform(parameter[i][1], parameter[i][2])
      
    dat, RMS = calc_CTRs(parameter, param_usage, dat, cell, bulk, surface, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag)

    R_track =Num.array([RMS],float)
    param_track = [parameter]

    
    fit_param = []
    for i in parameter.keys():
        if parameter[i][3]: fit_param.append(i)

    guess = (int(Num.log(float(Tend)/float(Tstart))/Num.log(cool))+1)*maxrun*1.25*len(fit_param)
    print 'approximated number of iterations: '+str( int(guess) )
    print 'R start = '+str(RMS)

    #counters
    check_better = False
    Random = 0
    better = 0
    rejected = 0   
    while Tstart > Tend:
        z = 1
        while (z <= maxrun):
            fit_param_tmp = fit_param[:]
            while len(fit_param_tmp) >0:
                RMS_tmp = 0.
                param_tmp = {}

                item = random.randint(0, len(fit_param_tmp)-1 )
                a = fit_param_tmp[item]
                addvector =  random.uniform((parameter[a][1]-parameter[a][0])*(MC * Tstart), (parameter[a][2]-parameter[a][0])*(MC * Tstart))

                if parameter[a][0] + addvector < parameter[a][1]: addvector = parameter[a][1] - parameter[a][0]
                if parameter[a][0] + addvector > parameter[a][2]: addvector = parameter[a][2] - parameter[a][0]

                for i in parameter.keys():
                    if i != a:
                        param_tmp[i] = [parameter[i][0]]
                    else:
                        param_tmp[i] = [parameter[i][0]+addvector]

                dat, RMS_tmp = calc_CTRs(param_tmp, param_usage, dat, cell, bulk, surface, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag)

                dR = (RMS - RMS_tmp) * factor

                if dR > 0.00001:
                    for i in parameter.keys():
                        parameter[i][0] = param_tmp[i][0]
                    better = better+1
                    R_track = Num.append(R_track,RMS_tmp)
                    param_track.append(parameter)
                    RMS = RMS_tmp
                    check_better = True
                    print a+',   better    '+str(RMS_tmp)
                elif dR <= 0.00001:
                    Boltz = exp(dR / Tstart)
                    Rand = random.uniform(0,1)
                    if(Boltz > Rand):
                        for i in parameter.keys():
                            parameter[i][0] = param_tmp[i][0]
                        Random = Random+1
                        R_track = Num.append(R_track,RMS_tmp)
                        param_track.append(parameter)
                        RMS = RMS_tmp
                        check_better = False
                        print a+',   random    '+str(RMS_tmp)
                    elif Boltz <= Rand:
                        check_better = False
                        rejected = rejected+1
                        print a+',   rejected   '+str(RMS_tmp)
                
                if not check_better:
                    del fit_param_tmp[item]
            z = z+1
        Tstart = Tstart * cool
        print '\n##################################################'
        print 'Temperature: '+str(Tstart/cool)
        print 'R = '+str(RMS)
        print '##################################################'+'\n'

    print '****************************'
    print 'Number of iterations: '+str(Random+ better+ rejected)
    print 'random: ' + str(Random)
    print 'better: ' + str(better)
    print 'rejected: '+ str(rejected)

    mini = Num.where(R_track == R_track.min())
    param_best = param_track[int(mini[0][0])]

    data_best, RMS_best = calc_CTRs(param_best, param_usage, dat, cell, bulk, surface, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag)
    R_track = Num.append(R_track, RMS_best)

    print '\n####################################################\n'
    print 'the best fit R = '+str(RMS_best)+'\n'
    print '*************************************\n'

    if plot_RMS_track:
        figure(3)
        clf()
        title('Development of R during fit')
        plot(range(len(R_track)),R_track,'ro')

    return data_best, param_best, RMS_best
################################################################################################################################
def plot_rods(dat, plot_dims, plot_bulk, plot_surf, plot_rough,plot_water, RMS):
    fig1 = figure(1, figsize = [15,9])
    clf()
    fig1.suptitle('R = '+str(round(RMS,7)), fontsize=20)
    for i in range(len(dat)):
        pl = str(plot_dims[0])+str(plot_dims[1])+str(i+1)
        fig1.add_subplot(pl)
        tmp = dat[i]
        if plot_bulk: plot(tmp.L,tmp.bulk,'g')
        if plot_surf: plot(tmp.L,tmp.surf,'r')
        if plot_water: plot(tmp.L,tmp.water,'m')
        if plot_rough: plot(tmp.L,tmp.rough,'c')
        plot(tmp.L,tmp.Fcalc,'k')
        errorbar(tmp.L,tmp.F,tmp.Ferr, fmt = 'bo')
        title(str(int(tmp.H))+str(int(tmp.K))+'L')
        semilogy()
#########################################################################################
def erf(x):
    # save the sign of x
    sign = 1
    if x < 0: 
        sign = -1
    x = abs(x)

    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*Num.exp(-x*x)
    return sign*y # erf(-x) = -erf(x)

##########################################################################################
def plot_edensity(surface, param, param_use, cell, database, rigid_bodies, use_bulk_water):
    edens = Num.zeros((1000),float)
    zwater, sig_water, Scale,specScale, beta, surface = param_unfold(param, param_use, surface, use_bulk_water)
    surface = RB_update(rigid_bodies, surface, param, cell)
    zmin = 0.
    zmax = 0.
    for i in range(len(surface)):
        if surface[i][3] > zmax: zmax = surface[i][3]
        if surface[i][3] < zmin: zmin = surface[i][3]
    zmin = (zmin - 0.1)*cell[2]
    zmax = (zmax + 0.5)*cell[2]
    abscissa = Num.arange(zmin,zmax, ((zmax-zmin)/1000))
    sig = (sig_water/(8*Num.pi**2))
    Auc = cell[0]* Num.sin(Num.radians(cell[5]))* cell[1]
    if use_bulk_water:
        for i in range(len(abscissa)):
            edens[i] = edens[i] + (0.5 *(1+ erf((abscissa[i])/(2*sig)**0.5)))*0.33456
    for x in surface:
        f_par = database[str.lower(x[0])]
        f = (f_par[0] + f_par[2] +f_par[4]+ f_par[6]+ f_par[8]) * x[10]
        for i in range(len(abscissa)):
            edens[i] = edens[i] + ((2*Num.pi*x[6])**(-1.5)* Num.exp(-0.5*(abscissa[i]-cell[2]*x[3])**2/x[6]))* f*2*Num.pi*x[6]/Auc
    
    figure(2)
    clf()
    plot(abscissa, edens, 'r')
    xlabel('z [Angstroem]')
    ylabel('electron density [Angstroem**(-3)]')
##########################################################################################
def write_par(parameter, param_labels, filename = 'parameters.new'):
    f = file(filename, 'w')
    f.write('% param_label              start          min          max    refine_flag\n')
    for i in param_labels:
        line = "%13s %18.12f %12.4f %12.4f %10s\n" % (i,parameter[i][0],
                                                                  parameter[i][1],
                                                                  parameter[i][2],
                                                                  parameter[i][3])
        f.write(line)
    f.close()
###########################################################################################
def write_data(data, filename = 'result.dat'):
    f = file(filename, 'w')
    f.write('    H     K        L           F        Ferr       Fcalc       Fbulk       Fsurf      Frough      Fwater\n')
    for i in data:
        for j in range(len(i.L)):
            line = "%5.2f %5.2f %8.4f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n" % (i.H, i.K, i.L[j], i.F[j], i.Ferr[j], i.Fcalc[j],
                                                                                             i.bulk[j], i.surf[j], i.rough[j], i.water[j])                                                                                        
            f.write(line)
    f.close()
############################################################################################
def check_model_consistency(param_labels, parameter, parameter_usage, rigid_bodies, use_bulk_water):
    used_params = []
    for i in param_labels:
        if parameter[i][3]: used_params.append(i)
    if use_bulk_water:
        global_parameters = ['zwater','sig_water','Scale','specScale','beta']
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
            print 'Model Inconsistency: parameter '+i+' shell be refined in the fit but is not defined in the model'
            global_flag = False
    return global_flag      
############################### Execution of fitting #################################################################
### input file information
#datafile = 'CTR2.dat' # same format as ROD datafile
#bulkfile = 'bul_calcite.bul'
#surfacefile = 'bul_calcite.sur'
#parameterfile = 'parameters.dat'
#rigidbodyfile = 'rigid_bodies.dat'
#do_write_surface = True
#do_write_bulk = True

### plotting options
#rod_plot = True
#plot_dims = [2,3]

#plot_compare = True
#plot_bulk = False
#plot_surf = False
#plot_rough = False
#plot_water = False

#plot_edens = True

#Rod_weight = [1,1,1,1,1]

### simulated annealing parameters 
#sim_an01 = True ##all parameters are chenged at once in each iteration
#sim_an02 = False  ##parameters are changed one by one (all parameters in one iteration)
#Tstart = 20
#Tend = 15
#cool = 0.7
#maxrun = 1
#MC = 1.0/100
#factor = 500000
#random_parameters = False
#plot_RMS_track = False

#bulk, cell, NLayers = read_bulk(bulkfile)
#surface, parameter_usage = read_surface(surfacefile)
#data = read_data(datafile)
#parameter, param_labels = read_parameters(parameterfile)
#rigid_bodies = read_rigid_bodies(rigidbodyfile)

#if sim_an01:
#    sim_an_params = [Tstart, Tend, cool, maxrun, MC, factor,random_parameters]
#    data, param_best, RMS_best = simulated_annealing01(data, cell, NLayers, bulk, surface,  database, Rod_weight, sim_an_params, parameter, parameter_usage, plot_RMS_track, rigid_bodies)

#if sim_an02:
#    sim_an_params = [Tstart, Tend, cool, maxrun, MC, factor,random_parameters]
#    data, param_best, RMS_best = simulated_annealing02(data, cell, NLayers, bulk, surface,  database, Rod_weight, sim_an_params, parameter, parameter_usage, plot_RMS_track, rigid_bodies)


#write_par(param_best, param_labels)
#write_cif(cell, surface, param_best, parameter_usage, rigid_bodies)
#if rod_plot: plot_rods(data, plot_dims, plot_bulk, plot_surf, plot_compare,plot_rough,plot_water, RMS_best)
#if plot_edens: plot_edensity(surface, param_best, parameter_usage, cell, database, rigid_bodies)
#if do_write_surface: write_surface(cell, surface, param_best, parameter_usage, rigid_bodies)
#if do_write_bulk: write_bulk(bulk, param_best)



    
    
    



