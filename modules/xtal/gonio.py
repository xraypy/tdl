##########################################################################
"""
Tom Trainor (tptrainor@alaska.edu) 
gonio calcs

Modifications:
--------------
- Mostly copied from matlab files, Frank Heberling (Frank.Heberling@ine.fzk.de)

"""
##########################################################################
"""
Todo

"""
##########################################################################

import numpy as num
import sys
import types

##########################################################################

def calc_Z(angles):
    #calculate goniometer rotation Matrix Z
    #angles is an array of diffrractometer angles as in spec file (del,eta,chi,phi,nu,mu)
    phi = num.radians(angles[3])
    chi = num.radians(angles[2])
    eta = num.radians(angles[1])
    mu  = num.radians(angles[5])
    Phi = num.array([[cos(phi),sin(phi),0],[-sin(phi),cos(phi),0],[0,0,1]],float)
    Chi = num.array([[cos(chi),0,sin(chi)],[0,1,0],[-sin(chi),0,cos(chi)]],float)
    Eta = num.array([[cos(eta),sin(eta),0],[-sin(eta),cos(eta),0],[0,0,1]],float)
    Mu  = num.array([[1,0,0],[0,cos(mu),-sin(mu)],[0,sin(mu),cos(mu)]],float)
    Z = num.dot(num.dot(num.dot(Mu,Eta),Chi),Phi)
    return Z

def calc_Q(angles, lam):
    #calculate k_in, k_out, and momentum transfer, Q, in lab frame, and two-theta
    #angles is an array of diffrractometer angles as in spec file (del,eta,chi,phi,nu,mu)
    k_in = 2* num.pi / lam * num.array([0,1,0])
    delta = num.radians(angles[0])
    nu = num.radians(angles[4])
    k_out = 2* num.pi / lam * num.array([sin(delta),cos(nu)*cos(delta),sin(nu)*cos(delta)])
    Q = k_out - k_in
    tth = vec_angle(k_in, k_out, num.eye(3))
    return k_in, k_out, Q, tth

def calc_B(cell):
    #calculate the B matrix
    #
    # the matrix B transforms the indicies of the 
    # reciprocal vector h to cartesian coordinates -->hc
    # where the cartesian system has been chosen so that:
    # x is parrallel to b1, 
    # y is in the plane of b1 and b2, 
    # z is perp to the plane of b1 and b2 
    # (see Busing and Levy)
    # Therefore:   
    #              h_c = B*h
    ##########################################################
    a3 = cell[2]
    A1 = num.radians(cell[3])
    G, G_rez, cell_rez = calc_G(cell)
    b1 = cell_rez[0]
    b2 = cell_rez[1]
    b3 = cell_rez[2]
    B1 = num.radians(cell_rez[3])
    B2 = num.radians(cell_rez[4])
    B3 = num.radians(cell_rez[5])
    B = num.array([[b1, b2*cos(B3), b3*cos(B2)],
                   [0, b2*sin(B3), -b3*sin(B2)*cos(A1)],
                   [0, 0, 1/a3]])
    return B

def calc_U(h1, h2, angles1, angles2, lambda1, lambda2, cell):
    #calculates the orientation matrix, U, from two reflectons h1 ,h2 and corresponding
    #diffractometer angles, lambdas and a given crystal coordinate system
    B = calc_B(cell) #calc B Matrix, carthesian crystal coords.
    #calc goniom. rot. Matrix and Momentum Transfer for both reflection
    Z1 = calc_Z(angles1)
    kin1, kout1, Q1, tth1 = calc_Q(angles1, lambda1)
    Z2 = calc_Z(angles2)
    kin2, kout2, Q2, tth2 = calc_Q(angles2, lambda2)
    #calc phi frame coords for diffraction vectors
    vphi_1 = num.dot(num.linalg.inv(Z1), (Q1/(2*num.pi)))
    vphi_2 = num.dot(num.linalg.inv(Z2), (Q2/(2*num.pi)))
    #calc carthesion coords of reflections
    hc1 = num.dot(B, h1)
    hc2 = num.dot(B, h2)
    #define tc vectors from hc vectors 
    tc1 = hc1 / vec_norm(hc1, num.eye(3))
    tc3 = num.cross(tc1, hc2) / vec_norm(num.cross(tc1, hc2), num.eye(3))
    tc2 = num.cross(tc3, tc1) / vec_norm(num.cross(tc3, tc1), num.eye(3))
    #define tphi vectors from vphi vectors
    tphi_1 = vphi_1 / vec_norm(vphi_1,num.eye(3))
    tphi_3 = num.cross(tphi_1,vphi_2) / vec_norm(num.cross(tphi_1,vphi_2), num.eye(3))
    tphi_2 = num.cross(tphi_3,tphi_1) / vec_norm(num.cross(tphi_3,tphi_1), num.eye(3))
    #define cooresponding matrices
    Tc = num.transpose(num.array([tc1,tc2,tc3]))
    Tphi = num.transpose(num.array([tphi_1,tphi_2,tphi_3]))
    #calc orientation matrix U and UB
    U = num.dot(Tphi, num.linalg.inv(Tc))
    UB = num.dot(U,B)
    return U, UB