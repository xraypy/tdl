##########################################################################
"""
Tom Trainor (fftpt@uaf.edu ) 
Xtal calcs

Modifications:
--------------
- convert matlab to python, Kunal Tanwar

"""
##########################################################################
"""
Module for performing xtal calcs 
List of available functions

* g_matx(cell):
  calculate g matrix from real space lattice parameters

* rlat(cell):
  calculate reciprocal space lattice parameters

* reciprocal_to_real(rcell):
  calculate real space basis from given reciprocal cell

* dot_product(u,v,cell):
  calculate dot product of two vectors u and v in basis defined by cell

* vector_mag(v,cell) :
  calculate magnitude of a vector v in basis defined by cell

* vector_angle(u,v,cell) :
  calculate angle between two vectors u and v in basis defined by cell

* d_spacing(hkl,cell) :
  calculate d-spacing for given (hkl) and cell parameters

* calc_d(hkl,cell):
  calculate vector d in bulk real space basis which has a magnitude of d_spacing

* basis_transform_cart(cell) :
  calculate basis and coordinate transformation matrix i.e. Gc and Fc to
  transform to cartesian space from a basis defined by cell

* transform_to_cart(u,v,w,cell) :
  transform given u,v,w coordinates from basis defined by cell to cartesian basis

* cart_to_cell(x,y,z,cell) :
  transform cartesian x,y,z coordinates to a basis defined by given cell 

* trans_hexa_to_rhombo(Va):
  transfrom any vector Va from hexa to rhombo system

* trans_rhombo_to_hexa(Va):
  transfrom any vector from rhombo to hexa system

* inplane_vectors(hkl,cell):
  calculate a number of lattice vectors lying in plane (hkl) of a lattice defined by
  cell. All the vectors are sorted according to their magnitudes given in column 4
  of the output. The first three column's give coefficients of the vector.

* slab_repeat_vectors(hkl,cell):
  calculate a number of possible slab repeat vectors by given (hkl) and cell sorted
  by the angle they make with surface normal that points into bulk

* basis_transformation_matrix(Va,Vb,Vc):
  transform to a basis defined by Va, Vb, and Vc,
    returns matrix MGFN,
    where MGFN[0] = M, M transforms the basis,
    MGFN[1] = G, does reverse of M
    MGFN[2] = F, transforms any vector in Va,Vb,Vc space
    MGFN[3] = N, does reverse of F
    the vectors Va,Vb,Vc are the real space vectors in the 
    original basis system which describe the new basis
    vectors for the primed system ie;
    a" = Va(1) a + Va(2) b + Va(3) c
    b" = Vb(1) a + Vb(2) b + Vb(3) c
    c" = Vc(1) a + Vc(2) b + Vc(3) c
    therefore M is the matrix describing the basis transform
    this is the same matrix which transforms the recip space indicies from the 
    original to the primed system ie 
    |h"|      |h|
    |k"| = M* |k|
    |l"|      |l|

    the inverse of M i.e. G gives the opposite relation

    |h|      |h"|
    |k| = G* |k"|
    |l|      |l"|
    and
    a = G(1,1) a" + G(1,2) b" + G(1,3) c"
    b = G(2,1) a" + G(2,2) b" + G(2,3) c"
    c = G(3,1) a" + G(3,2) b" + G(3,3) c"

    the matrix which transforms the recip space basis is
    given by the inverse transpose of M i.e. F

    % ie 
    a*" = F(1,1) a* + F(1,2) b* + F(1,3) c*
    b*" = F(2,1) a* + F(2,2) b* + F(2,3) c*
    c*" = F(3,1) a* + F(3,2) b* + F(3,3) c*

    this is the same matrix which transforms the real space
    vector indicies from original to primed system

    |x"|      |x|
    |y"| = F* |y|
    |z"|      |z|

    the inverse of F gives the opposite relation
    N = inv(F);

    ie 
    a* = N(1,1) a*" + N(1,2) b*" + N(1,3) c*"
    b* = N(2,1) a*" + N(2,2) b*" + N(2,3) c*"
    c* = N(3,1) a*" + N(3,2) b*" + N(3,3) c*"

    this is the same matrix which transforms the real space
    vector indicies from primed to original system
    |x|      |x"|
    |y| = N* |y"|
    |z|      |z"|


* calc_surf_cell(F,bulk_name,cell):
  calculate a big sur file from the bulk structure file, the extension of the new
  file is .surf, and it needs to filtered by hand to get the fractional coordinates
                                  
"""
##########################################################################
"""
Todo
 - get xyz file from fit and par files and perform bond valence calcs
 - convert to class
"""
##########################################################################

import numpy as num
import sys
import string

##########################################################################
def g_matx(cell): # calculate g matrix from real space lattice parameters
    a = cell[0]
    b = cell[1]
    c = cell[2]
    alp = cell[3]
    bet = cell[4]
    gam = cell[5]
    Pi  = num.pi
    g = num.array([[a*a, a*b*(num.cos(gam*Pi/180)), a*c*(num.cos(bet*Pi/180))],
                    [b*a*(num.cos(gam*Pi/180)), b*b, b*c*(num.cos(alp*Pi/180))],
                    [c*a*(num.cos(bet*Pi/180)), c*b*(num.cos(alp*Pi/180)), c*c]])
    return g

##########################################################################
def rlat(cell): # calculate reciprocal space lattice parameters
    g = g_matx(cell)
    g_inv = num.linalg.inv(g)
    a_r = num.sqrt(g_inv[0,0])
    b_r = num.sqrt(g_inv[1,1])
    c_r = num.sqrt(g_inv[2,2])
    Pi  = num.pi
    alp_r = (180/Pi)*num.arccos(g_inv[1,2]/(b_r*c_r))
    bet_r = (180/Pi)*num.arccos(g_inv[0,2]/(a_r*c_r))
    gam_r = (180/Pi)*num.arccos(g_inv[0,1]/(a_r*b_r))
    rcell = [a_r, b_r, c_r, alp_r, bet_r, gam_r]
    return rcell

##########################################################################
def dot_product(u,v,cell):# calculate dot product of two vectors in basis defined by cell
    g = g_matx(cell)
    x = num.matrixmultiply(g,num.transpose(v)) 
    dot_product = num.matrixmultiply(u,x) # dot product = u*g*v', from Sands.
    return dot_product

##########################################################################
def vector_mag(v,cell): # calculate magnitude of a vector in basis defined by cell
    mag   = dot_product(v,v,cell)
    v_mag = num.sqrt(mag)
    return v_mag

##########################################################################
def vector_angle(u,v,cell): # calculate angle between two vectors in basis defined by cell
    x = dot_product(u,v,cell)
    u_mag = vector_mag(u,cell)
    v_mag = vector_mag(v,cell)
    Pi = num.pi
    v_angle = (180/Pi)*num.arccos(x/(u_mag*v_mag))
    return v_angle

##########################################################################
def d_spacing(hkl,cell): # calculate d-spacing for given h,k,l and cell parameters 
    rcell = rlat(cell)
    H = vector_mag(hkl,rcell)
    d_spacing = 1/H
    return d_spacing

##########################################################################
def calc_d(hkl,cell): # calculate vector d in bulk real space basis which has a magnitude of d_spacing
    g = g_matx(cell)
    g_inv = num.linalg.inverse(g)
    hkl = num.transpose(hkl)
    d_spc = d_spacing(hkl,cell)
    d = num.array(num.matrixmultiply(g_inv,hkl))
    d = d_spc*d_spc*d
    return d

##########################################################################
def basis_transform_cart(cell): #calculate Gc and Fc for transformations to cartesian basis
    a = cell[0]
    b = cell[1]
    c = cell[2]
    alp = cell[3]
    bet = cell[4]
    gam = cell[5]
    Pi = num.pi
    rcell = rlat(cell)
    a_r = rcell[0]
    b_r = rcell[1]
    c_r = rcell[2]
    alp_r = rcell[3]
    bet_r = rcell[4]
    gam_r = rcell[5]
    # Gc transforms basis to cartesian
    Gc = num.array([[1/a, 0 , 0],
                     [-1/(a*num.tan(gam*Pi/180)), 1/(b*num.sin(gam*Pi/180)), 0],
                     [a_r*num.cos(bet_r*Pi/180), b_r*num.cos(alp_r*Pi/180), c_r]])
    # F transforms coordinates to cartesian
    # F = inverse(G')
    Fc = num.array([[a, b*num.cos(gam*Pi/180), c*num.cos(bet*Pi/180)],
                     [0, b*num.sin(gam*Pi/180), -1*c*num.sin(bet*Pi/180)*num.cos(alp_r*Pi/180)],
                     [0, 0, 1/c_r]])
    return Fc

##########################################################################
def transform_to_cart(u,v,w,cell): # transform u,v,w to cartesian
    Fc = basis_transform_cart(cell) 
    uvw = [u, v, w]
    uvw_xyz = num.matrixmultiply(Fc, num.transpose(uvw))
    xyz = num.transpose(uvw_xyz1)
    return xyz

##########################################################################
def cart_to_cell(x,y,z,cell): # transform x,y,z to basis defined by cell
    Fc = basis_transform_cart(cell)
    xyz = [x, y, z]
    xyz_uvw = num.matrixmultiply((num.linalg.inverse(Fc)),num.transpose(xyz))
    uvw = num.transpose(xyz_uvw)
    return uvw

##########################################################################
def trans_hexa_to_rhombo(Va): # transfrom any vector from hexa to rhombo system
    MM = num.array([[0.6667, 0.3333, 0.3333],
                      [-0.3333, 0.3333, 0.3333],
                      [-0.3333, -0.6667, 0.3333]])
    NN = num.linalg.inv(num.transpose(MM))
    V_rh = num.matrixmultiply(NN,num.transpose(Va))
    return V_rh

##########################################################################
def trans_rhombo_to_hexa(Va): # transfrom any vector from rhombo to hexa system
    MM = num.array([[0.6667, 0.3333, 0.3333],
                      [-0.3333, 0.3333, 0.3333],
                      [-0.3333, -0.6667, 0.3333]])
    NN = num.linalg.inv(num.transpose(MM))
    NN_inv = num.linalg.inv(NN)
    V_hx = num.matrixmultiply(NN_inv,num.transpose(Va))
    return V_hx
    
##########################################################################
def inplane_vectors(hkl, cell): #calculate surface vectors for a given hkl plane through origin
    V_tmp=[]
    for n1 in range(-4,5):
        for n2 in range(-4,5):
            for n3 in range(-4,5):
                temp = n1*hkl[0] + n2*hkl[1] + n3*hkl[2]
                if temp == 0:
                    n1n2n3 = [n1, n2, n3]
                    v_mag = vector_mag(n1n2n3,cell)
                    v = [n1, n2, n3, v_mag] 
                    V_tmp.append(v)
    V_tmp = num.array(V_tmp)
    # V_tmp has all the surface vectors, now sort them according to magnitude
    idx = num.argsort(V_tmp[:,3])
    Vs = []
    for n in range(0,len(idx)):
        x = V_tmp[idx[n],:]
        Vs.append(x)
    Vs = num.array(Vs)
    return Vs
     
##########################################################################
def slab_repeat_vectors(hkl,cell): # calculate number of possible slap repeat vectors
    g = g_matx(cell)
    g_inv = num.linalg.inverse(g)
    hkl = num.transpose(hkl)
    d_spc = d_spacing(hkl,cell)
    d = num.array(num.matrixmultiply(g_inv,hkl))
    d = d_spc*d_spc*d
    V_tmp=[]
    for n1 in range(-4,5):
        for n2 in range(-4,5):
            for n3 in range(-4,5):
                temp = n1*hkl[0] + n2*hkl[1] + n3*hkl[2]
                if temp < 0:
                    n1n2n3 = [n1, n2, n3]
                    v_mag = vector_mag(n1n2n3,cell)
                    v_ang = vector_angle(n1n2n3, -d, cell)
                    v = [n1, n2, n3, v_mag, temp, v_ang] 
                    V_tmp.append(v)
    V_tmp = num.array(V_tmp)
   # V_tmp has lot of possible slab vectors, now sort them according to angle they make with -d
    idx = num.argsort(V_tmp[:,5])
    Vb = []
    for n in range(0,len(idx)):
        x = V_tmp[idx[n],:]
        Vb.append(x)
    Vb = num.array(Vb)
    return Vb

##########################################################################
def basis_transformation_matrix(Va,Vb,Vc): # transform to basis defined by Va, Vb and Vc
    # M transforms the basis to surface basis see Trainor et al. 2002
    #        xas   yas  zas
    # [M] =  xbs   ybs  zbs 
    #        xcs   ycs  zcs
    #
    M = num.array([Va,Vb,Vc])
    # G does the opposite of M
    G = num.linalg.inv(M)
    # F transforms any vector to surface (as,bs,cs) system
    # for example V_rs = F*Vr, Vr in repeat vector in bulk system and Vr_s in surface system
    F = num.linalg.inv(num.transpose(M))
    # N does opposite of F
    N = num.linalg.inv(F)
    MGFN = num.array([M,G,F,N])
    return MGFN

##########################################################################
def calc_surf_cell(F,bulk_name,cell): # calculate the big surface cell
    # read bulk file
    file = open(bulk_name,'r')
    lines = file.readlines()
    file.close()
    # first line is title 
    title = lines.pop(0)
    # second line is cell parameters
    cell_str = lines.pop(0)
    xyz_list = []
    el_list = []
    el = None
    for m in range(0,len(lines)):
        tmp = string.split(lines[m])
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        el = tmp[0]
        el_list.append(el)
        xyz_a = [x,y,z]
        xyz_list.append(xyz_a)
    # use p to check if the file is read and sorted correctly
    # and use it later while writing out the file
    p = len(xyz_list)
    xyz = xyz_list
    el_name = el_list
    xyz_list = num.array(xyz_list)
    xyz_tmp = []
    el = None
    el_tmp=[]
    for xt in range(0,8):
        for yt in range(0,8):
            for zt in range(-1,2):
                for j in range(0,len(el_list)):
                    xyz_tmp = [xyz_list[j,0]+ xt, xyz_list[j,1]+yt, xyz_list[j,2]+zt]
                    xyz.append(xyz_tmp)
                    el = el_list[j]
                    el_tmp.append(el)
    xyz = num.array(xyz)
    for j in range(0,len(el_tmp)):
        el_name.append(el_tmp[j])
    for j in range(0,len(el_name)):
        v_tmp = num.array([xyz[j,0],xyz[j,1],xyz[j,2]])
        v = num.matrixmultiply(F,num.transpose(v_tmp))
        xyz[j,0] = v[0]
        xyz[j,1] = v[1]
        xyz[j,2] = v[2]
    fname_out = bulk_name + ".surf"
    file = open(fname_out,'w')
    s = '%-5i\n' % len(el_name)
    file.write(s)
    s = '%-5s\n' % fname_out
    file.write(s)
    for j in range(0,len(el_name)):
        if 0 <(j+1)%p < p:
            s = '%-5s    %7.5f    %7.5f    %7.5f  %5i\n' % (el_name[j],xyz[j,0],xyz[j,1],xyz[j,2],(j+1)%p)
            file.write(s)
        elif (j+1)%p==0:
            s = '%-5s    %7.5f    %7.5f    %7.5f  %5i\n' % (el_name[j],xyz[j,0],xyz[j,1],xyz[j,2],p)
            file.write(s)
    file.close()

##########################################################################
    

