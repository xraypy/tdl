##########################################################################
"""
Tom Trainor (fftpt@uaf.edu ) 
Surface indexing

Modifications:
--------------
 - clean up and convert to class

"""

##########################################################################

import numpy as num
import xtal_calc

##########################################################################
def inplane_vectors(hkl, cell):
    """
    * inplane_vectors(hkl,cell):
      calculate a number of lattice vectors lying in plane (hkl) of a lattice defined by
      cell. All the vectors are sorted according to their magnitudes given in column 4
      of the output. The first three column's give coefficients of the vector.
    """
    #calculate surface vectors for a given hkl plane through origin
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
def slab_repeat_vectors(hkl,cell):
    """
    * slab_repeat_vectors(hkl,cell):
      calculate a number of possible slab repeat vectors by given (hkl) and cell sorted
      by the angle they make with surface normal that points into bulk
    """
    # calculate number of possible slap repeat vectors
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
def basis_transformation_matrix(Va,Vb,Vc):
    """
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

    """
    """
    # transform to basis defined by Va, Vb and Vc
    # M transforms the basis to surface basis see Trainor et al. 2002
    #        xas   yas  zas
    # [M] =  xbs   ybs  zbs 
    #        xcs   ycs  zcs
    #
    """
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
def calc_surf_cell(F,bulk_name,cell):
    """
    * calc_surf_cell(F,bulk_name,cell):
      calculate a big sur file from the bulk structure file, the extension of the new
      file is .surf, and it needs to filtered by hand to get the fractional coordinates
                                      
    """
    # calculate the big surface cell
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
##########################################################################
if __name__ == "__main__":
    # define cell and hkl:
    cell = [a,b,c,alpha,beta,gam]
    hkl = [h,k,l]
    # calculate d spacing in Angstroms
    d_space = xtal.d_spacing(hkl,cell)

    # calculate vector d in bulk real space basis
    d = xtal.calc_d(hkl,cell)

    # calculate inplane vectors to choose Va and Vb
    Vs = xtal.inplane_vectors(hkl,cell)

    # if indexing hematite check if Va and Vb can be further reduced by converting to rhomohedral from hexagonal
    Va_rh = xtal.trans_hexa_to_rhombo(Va)
    Vb_rh = xtal.trans_hexa_to_rhombo(Vb)

    # calculate slab repeat vectors to choose Vr and the terminating layer (n)
    Vb = xtal.slab_repeat_vectors(hkl,cell) # choose Vr and n from the list

    Vc = n*d

    # calculate all the basis transformation matrices
    MGFN = xtal.basis_transformation_matrix(Va,Vb,Vc)
    F = MGFN[2]
    # find repeat vector in surface basis
    Vr_s = num.matrixmultiply(F,num.transpose(Vr))
    del_1 = -1*Vr_s[0]
    del_2 = -1*Vr_s[1]
    mag_as = xtal.vector_mag(Va,cell)
    mag_bs = xtal.vector_mag(Vb,cell)
    mag_cs = xtal.vector_mag(Vc,cell)
    alp_s = xtal.vector_angle(Vb,Vc,cell)
    bet_s = xtal. vector_angle(Va,Vc,cell)
    gam_s = xtal.vector_angle (Va,Vb,cell)

    scell = [ mag_as, mag_bs, mag_cs, alp_s, bet_S, gam_s]

    # must notice repeat angle while choosing Vr but check again
    repeat_angle = 180 - xtal.vector_angle(Vr_s, [0,0,1], cell)

    # calculate the surface cell using the bulk file bulk_name
    # ouput file will be bulk_name.surf
    xtal.calc_surf_cell(F,bulk_name,cell)

