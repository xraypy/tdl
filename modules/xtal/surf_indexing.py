"""
Surface indexing

Authors/Modifications:
-----------------------
* Tom Trainor (tptrainor@alaska.edu)

Todo:
-----
* clean up and convert to class
* use space group to generate from bulk
  assymetric unit
"""

##########################################################################

import numpy as num
#import xtal_calc
from lattice import Lattice, LatticeTransform

##########################################################################
def surface_vectors(hkl, lat):
    """
    Calculate in-plane lattice vectors, and repeat vectors

    Parameters:
    -----------
    * hkl defines the plane
    * lat defines the lattice (is a Lattice instance)

    Output:
    -------
    The first three column's of the output give coefficients of the
    in plane vectors.
    
    All the vectors are sorted according to their magnitudes given 
    in column 4 of the output. 
    """
    # calculate vector d in bulk real space basis
    d = lat.dvec(hkl)

    #calculate surface vectors for a given hkl plane through origin
    Vs_tmp=[]
    Vr_tmp=[]
    vrange = range(-4,5)
    for n1 in vrange:
        for n2 in vrange:
            for n3 in vrange:
                temp = n1*hkl[0] + n2*hkl[1] + n3*hkl[2]
                # from law of rational indicies, temp is zero 
                # if the lattice vector is in the hkl plane
                if temp == 0:
                    v = [n1, n2, n3]
                    v_mag = lat.mag(v,lat)
                    Vs_tmp.append([n1, n2, n3, v_mag] )
                # from law of rational indicies, temp < zero 
                # if the lattice vector points below the hkl plane
                # and temp is the number of planes below the surface
                # where v terminates
                elif temp < 0:
                    v = [n1, n2, n3]
                    v_mag = lat.mag(v,lat)
                    v_ang = lat.angle(v,-1.*d,lat)
                    Vr_tmp.append([n1, n2, n3, v_mag,temp,v_ang])

    # Vs_tmp has all the surface vectors,
    # now sort them according to magnitude
    Vs_tmp = num.array(Vs_tmp)
    idx = num.argsort(Vs_tmp[:,3])
    Vs = []
    for n in range(0,len(idx)):
        x = Vs_tmp[idx[n],:]
        Vs.append(x)
    Vs = num.array(Vs)

    # Vr_tmp has lot of possible slab vectors,
    # now sort them according to angle they make with -d
    Vr_tmp = num.array(Vr_tmp)
    idx = num.argsort(Vr_tmp[:,5])
    Vr = []
    for n in range(0,len(idx)):
        x = Vr_tmp[idx[n],:]
        Vr.append(x)
    Vr = num.array(Vr)

    return Vs,Vr

##########################################################################
def calc_surf_cell(F,bulk_name,cell):
    """
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
    #cell = [a,b,c,alpha,beta,gam]
    #hkl = [h,k,l]
    
    cell = [10.223, 5.992, 4.761,90,90,90]
    hkl = [1.,0.,1.]
    lat = Lattice(*cell)
    
    # calculate d spacing in Angstroms
    d_space = lat.d(hkl)

    # calculate vector d in bulk real space basis
    d = lat.dvec(hkl)

    # calculate inplane vectors to choose Va and Vb
    Vs,Vr = surface_vectors(hkl,lat)

    #################
    # choose basis vectors for surface system
    Va = [1.,0.,-1.]
    Vb = [0.,1.,0.]
    Vc = 5*d
    Vrpt = [-1.,0.,-4.]
    trans = LatticeTransform(lat,Va=Va,Vb=Vb,Vc=Vc)
    slat = trans.plat()
    
    # surface cell params
    # should be same as those in slat
    """
    mag_as = lat.mag(Va)
    mag_bs = lat.mag(Vb)
    mag_cs = lat.mag(Vc)
    alp_s  = lat.angle(Vb,Vc)
    bet_s  = lat.angle(Va,Vc)
    gam_s  = lat.angle (Va,Vb)
    scell = [ mag_as, mag_bs, mag_cs, alp_s, bet_s, gam_s]
    """
    
    # find repeat vector in surface basis
    Vrpt_s = trans.xp(Vrpt)
    del_1 = -1*Vrpt_s[0]
    del_2 = -1*Vrpt_s[1]

    # check repeat angle 
    repeat_angle = slat.angle(Vrpt_s, [0.,0.,-1.])

    # calculate the surface cell using the bulk file bulk_name
    # ouput file will be bulk_name.surf
    #calc_surf_cell(bulk_name,lat,trans)

    # to compute hkl_s from hkl
    # use the following:
    #hkl_s = trans.hp([h,k,l])
    #
    # visa versa given hkl_s compute
    # hkl in the bulk basis from:
    #hkl = trans.h([hs,ks,ls])
