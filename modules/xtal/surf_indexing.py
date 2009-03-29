##########################################################################
"""
Tom Trainor (fftpt@uaf.edu ) 
Surface indexing

Modifications:
--------------


"""
##########################################################################

import numpy
import xtal

##########################################################################
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
Vr_s = numpy.matrixmultiply(F,numpy.transpose(Vr))
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

