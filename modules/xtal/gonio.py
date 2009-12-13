##########################################################################
"""
Tom Trainor (tptrainor@alaska.edu)
Frank Heberling (Frank.Heberling@ine.fzk.de)
gonio calcs

Modifications:
--------------
- Mostly copied from matlab files, FH

"""
##########################################################################
"""
Todo

"""
##########################################################################
"""
Notes

The following procedure can be used to define the orientation
matrix of a xtal mounted on a goniometer and to convert between
motor angles and hkl and visa-versa


Define the matrix B which transforms the indicies of the 
reciprocal lattice vector h=[h,k,l] to cartesian coordinates
hc = [xc,yc,zc] where the cartesian system has been chosen so that:
  x is parrallel to ar, 
  y is in the plane of ar and br, 
  z is perp to the plane of ar and br
  (see Busing and Levy)
Therefore:  
  hc = B*h
Note that hc (and the cartesian basis) is defined with respect to
the crystal lattice therefore is independant of the gonio angles.  





For psic geom from You's paper:
1. phi rotaion (P)
2. chi rotation (X)
3. eta roation (H)
4. mu rotation (M)
To calc the lab frame coords of h --> hm:
    hm = M*H*X*P*h_phi = M*H*X*P*U*B*h
 - or - Let Z = M*H*X*P, 
    hm = Z*h_phi
 hm gives the recip lattice vector
 position in the lab frame
 after the sample is rotated.


Calc Z matrix, Qm and ki and kr (all in the m-frame)
the scattering angle tth



note: calc kappa, keta and kphi
kap_alp = 50.031;
keta = eta - (180/pi)*asin(-tand(chi/2)/tand(kap_alp))
kphi = phi - (180/pi)*asin(-tand(chi/2)/tand(kap_alp))
kappa = (180/pi)*2*asin(sind(chi/2)/sind(kap_alp))


1. H You
2. Busy and Leving

"""
##########################################################################

import numpy as num
import types

from Lattice import lattice

##########################################################################
def angles(angles):
    phi = num.radians(angles[3])
    chi = num.radians(angles[2])
    eta = num.radians(angles[1])
    mu  = num.radians(angles[5])
    delta = num.radians(angles[0])
    nu = num.radians(angles[4])
    return {'phi':phi,'chi':chi,'eta':eta,'mu':mu,
            'delta':delta,'nu':nu}

class Psic:
    def __init__(self,a=1.,b=1.,c=1.,alp=90.,beta=90.,gamma=90.):
        self.lattice = Lattice(a,b,c,alp,bet,gamma)
        # primary reflection
        self.or1={'h':num.array([0.,0.,1.]),'lam':1.0,
                  'phi':0.0,'chi':0.0,'eta':0.0,'mu':0.0,
                  'delta':0.0,'nu':0.0}
        # secondary reflection
        self.or1={'h':num.array([0.,0.,1.]),'lam':1.0,
                  'phi':0.0,'chi':0.0,'eta':0.0,'mu':0.0,
                  'delta':0.0,'nu':0.0}
        
    def set_or1(self,h=None,lam=None,phi=None,chi=None,eta=None,
                mu=None,delta=None,nu=None):
        """
        Set / adjust the primary orientation reflection
        """
        if h!=None:     self.or1['h'] = num.array(h,dtype=float)
        if lam!= None:  self.or1['lam']=float(lam)
        if phi!=None:   self.or1['phi']=float(phi)
        if chi!=None:   self.or1['chi']=float(chi)
        if eta!=None:   self.or1['eta']=float(eta)
        if mu!=None:    self.or1['mu']=float(mu)
        if delta!=None: self.or1['delta']=float(delta)

    def set_or2(self,h=None,lam=None,phi=None,chi=None,eta=None,
                mu=None,delta=None,nu=None):
        """
        Set / adjust the secondary orientation reflection
        """
        if h!=None:     self.or2['h'] = num.array(h,dtype=float)
        if lam!= None:  self.or2['lam']=float(lam)
        if phi!=None:   self.or2['phi']=float(phi)
        if chi!=None:   self.or2['chi']=float(chi)
        if eta!=None:   self.or2['eta']=float(eta)
        if mu!=None:    self.or2['mu']=float(mu)
        if delta!=None: self.or2['delta']=float(delta)

    def swap_or(self,):
        """
        Swap the primary and secondary reflection
        """
        tmp = self.or1
        self.or1 = self.or2
        self.or2 = tmp

    def calc_U(self,h1, h2, angles1, angles2, lambda1, lambda2, cell):
        """
        Calculate the orientation matrix, U, from the primary and secondary
        reflectons and given lattice
        """
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

    def _calc_B(self,):
        """
        Calculate the B matrix
        """
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



