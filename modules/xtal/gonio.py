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
------
The following describes the calculations used to define the orientation
matrix of a xtal mounted on a goniometer and to convert between
motor angles and hkl and visa-versa

Define the matrix B which which transforms the indicies of the 
reciprocal lattice vector, h=[h,k,l], to a cartesian basis,
hc=[xc,yc,zc], where the cartesian system has been chosen so that:
  x is parrallel to ar, 
  y is in the plane of ar and br, 
  z is perp to the plane of ar and br
  (see Busing and Levy)
Therefore:  
  hc = B*h
Note that hc (and the cartesian basis) is defined with respect to
the crystal lattice therefore is independant of the gonio angles.  

U is the matrix that defines how the sample is mounted on the
diffractometer.  For example we define the lab frame coordinate
system such that:
  x is vertical (perpendicular, pointing to the ceiling of the hutch)
  y is directed along the incident beam path
  z make the system right handed and lies in the horizontal scattering plane
Therefore with all gonio angles set to zero the lab xyz axis are coincident
with the goniometer axes in a well defined way (depending on the axis definitions).
Under the instrument settings with all angles zero (phi frame) and with the
sample oriented in an arbitrary maner the matrix U is used
to calculate the lab frame indicies of an arbritrary recip lattice vector,
h=[h,k,l], such that:
   hphi = U*hc = U*B*h  
We can then orient this vector in the lab frame applying the rotation
matricies of the various gonio axes.  The algorithm to determine U
from a set of reflections and angles is given by Busing and Levy

For psic geom from You's paper the order of sample
axis rotations is:
1. phi rotation (matrix = P)
2. chi rotation (matrix = X)
3. eta rotation (matrix = H)
4. mu rotation  (matrix = M)

To calculate the lab frame coords of h --> hm:
    hm = M*H*X*P*hphi = M*H*X*P*U*B*h
Or letting Z = M*H*X*P, 
    hm = Z*hphi

Therefore, hm gives the indicies of the recip lattice vector, h,
in the lab frame cartesian basis after the sample is rotated.

The diffraction condition specifies that:
   Q = (2*pi)*h 
The lab frame coordinates of Q can be calc from:
   Qm = kr - ki
where ki and kr in the lab frame are given (for psic) by :
   ki = (2*pi/lam)*[ 0, 1, 0 ]
   kr = (2*pi/lam)*[ sin(del), cos(nu)*cos(del),sin(nu)*cos(del) ]
   Qm = kr - ki;
Therefore the diffraciton condition in the lab frame is
   Qm = (2*pi)*hm = (2*pi)*Z*hphi = (2*pi)*Z*U*B*h


1. H You
2. Busy and Leving

"""
##########################################################################

import numpy as num
import types

from mathutil import cosd, sind, tand
from mathutil import arccosd, arcsind, arctand
from lattice import Lattice

##########################################################################
def angles(angles):
    """
    Angles defined by spec
    
    Note: calc kappa, keta and kphi
    kap_alp = 50.031;
    keta = eta - (180/pi)*asin(-tand(chi/2)/tand(kap_alp))
    kphi = phi - (180/pi)*asin(-tand(chi/2)/tand(kap_alp))
    kappa = (180/pi)*2*asin(sind(chi/2)/sind(kap_alp))
    """
    phi = num.radians(angles[3])
    chi = num.radians(angles[2])
    eta = num.radians(angles[1])
    mu  = num.radians(angles[5])
    delta = num.radians(angles[0])
    nu = num.radians(angles[4])
    return {'phi':phi,'chi':chi,'eta':eta,'mu':mu,
            'delta':delta,'nu':nu}

class Psic:
    """
    Orientation calculations for Psic geometry.

    The default dummy orientation matrix is set up 
    assuming the sample is mounted such that (001) plane
    is perpendicular to the eta and phi rot axes
    (ie c-axis is parrallel to the eta and phi rot axes)
    and that the a-axis points back towards the incident beam,
    ie. the b-axis is parrallel to the nu and mu rot axes

    """
    def __init__(self,a=10.,b=10.,c=10.,alpha=90.,beta=90.,gamma=90.,lam=1.0):
        """
        Initialize by passing a,b,c in angstroms 
        alpha, beta, gamma in degrees,
        and lambda in angstroms
        """
        # set lattice
        self.lattice = Lattice(a,b,c,alpha,beta,gamma)
        # dummy primary reflection
        tth = self.lattice.tth([0.,0.,1.],lam=lam)
        self.or1={'h':num.array([0.,0.,1.]),
                  'phi':0.0,'chi':0.0,'eta':0.0,'mu':tth/2.,
                  'nu':tth,'delta':0.0,'lam':lam}
        # dummy secondary reflection
        tth = self.lattice.tth([0.,1.,0.],lam=lam)
        self.or2={'h':num.array([0.,1.,0.]),
                  'phi':0.0,'chi':0.0,'eta':tth/2.,'mu':0.0,
                  'nu':0.0,'delta':tth,'lam':lam}
        # compute initial matricies
        self.U = []
        self.B = []
        self._update()

    def __repr__(self,):
        lout = self.lattice.__repr__()
        # add or reflections...
        return lout
    
    def _update(self):
        self._calc_B()
        self._calc_U()

    def set_or1(self,h=None,phi=None,chi=None,eta=None,
                mu=None,nu=None,delta=None,lam=None):
        """
        Set / adjust the primary orientation reflection
        """
        if h!=None:     self.or1['h'] = num.array(h,dtype=float)
        if phi!=None:   self.or1['phi']=float(phi)
        if chi!=None:   self.or1['chi']=float(chi)
        if eta!=None:   self.or1['eta']=float(eta)
        if mu!=None:    self.or1['mu']=float(mu)
        if nu!=None:    self.or1['nu']=float(nu)
        if delta!=None: self.or1['delta']=float(delta)
        if lam!= None:  self.or1['lam']=float(lam)
        self._update()

    def set_or2(self,h=None,phi=None,chi=None,eta=None,
                mu=None,nu=None,delta=None,lam=None):
        """
        Set / adjust the secondary orientation reflection
        """
        if h!=None:     self.or2['h'] = num.array(h,dtype=float)
        if phi!=None:   self.or2['phi']=float(phi)
        if chi!=None:   self.or2['chi']=float(chi)
        if eta!=None:   self.or2['eta']=float(eta)
        if mu!=None:    self.or2['mu']=float(mu)
        if nu!=None:    self.or2['nu']=float(nu)
        if delta!=None: self.or2['delta']=float(delta)
        if lam!= None:  self.or2['lam']=float(lam)
        self._update()

    def swap_or(self,):
        """
        Swap the primary and secondary reflection
        """
        tmp = self.or1
        self.or1 = self.or2
        self.or2 = tmp
        self._update()

    def set_lat(self,a=None,b=None,c=None,alpha=None,
                beta=None,gamma=None,lam=None):
        """
        Update lattice parameters
        """
        self.lattice.update(a=a,b=b,c=c,alpha=alpha,
                            beta=beta,gamma=gamma,lam=lam)
        self._update()
    
    def _calc_B(self,):
        """
        Calculate the B matrix
        """
        (a,b,c,alp,bet,gam) = self.lattice.cell()
        (ar,br,cr,alpr,betr,gamr) = self.lattice.rcell()

        B = num.array([[ar,  br*cosd(gamr),     cr*cosd(betr)        ],
                       [0.,  br*sind(gamr),  -cr*sind(betr)*cosd(alp)],
                       [0.,      0.,                 1./c            ]])
        self.B = B

    def _calc_U(self,):
        """
        Calculate the orientation matrix, U, from the primary and secondary
        reflectons and given lattice
        """
        # use these, note they are used below on vectors
        # defined in the cartesian lab frame basis
        cross = num.cross
        norm = num.linalg.norm
        
        #calc Z and Q for the OR reflections
        Z1 = self.calc_Z(self.or1['phi'],self.or1['chi'],
                         self.or1['eta'],self.or1['mu'])
        (Q1,ki1,kr1) = self.calc_Q(self.or1['nu'],self.or1['delta'],self.or1['lam'],)

        Z2 = self.calc_Z(self.or2['phi'],self.or2['chi'],
                         self.or2['eta'],self.or2['mu'])
        (Q2,ki2,kr2) = self.calc_Q(self.or2['nu'],self.or2['delta'],self.or2['lam'],)
        
        # calc the phi frame coords for diffraction vectors
        # note divide out 2pi since the diffraction is 2pi*h = Q
        vphi_1 = num.dot(num.linalg.inv(Z1), (Q1/(2.*num.pi)))
        vphi_2 = num.dot(num.linalg.inv(Z2), (Q2/(2.*num.pi)))
        
        #calc cartesian coords of h vectors
        hc_1 = num.dot(self.B, self.or1['h'])
        hc_2 = num.dot(self.B, self.or2['h'])

        #So at this point the following should be true:
        #     vphi_1 = U*hc_1
        #     vphi_2 = U*hc_2
        # and we could use these relations to solve for U.
        # But U solved directly from above is likely not to be orthogonal
        # since the angles btwn (vphi_1 and vphi_2) and (hc_1 and hc_2) are 
        # not exactly the same due to exp errors.....
        # Therefore, get an orthogonal solution for U from the below treatment
        
        #define the following normalized vectors from hc vectors 
        tc_1 = hc_1 / norm(hc_1)
        tc_3 = cross(tc_1, hc_2) / norm(cross(tc_1, hc_2))
        tc_2 = cross(tc_3, tc_1) / norm(cross(tc_3, tc_1))

        #define tphi vectors from vphi vectors
        tphi_1 = vphi_1 / norm(vphi_1)
        tphi_3 = cross(tphi_1,vphi_2) / norm(cross(tphi_1,vphi_2))
        tphi_2 = cross(tphi_3,tphi_1) / norm(cross(tphi_3,tphi_1))

        #define the following matrices
        Tc   = num.transpose(num.array([tc_1,tc_2,tc_3]))
        Tphi = num.transpose(num.array([tphi_1,tphi_2,tphi_3]))
        
        # calc orientation matrix U
        # note either of the below work since Tc is orthogonal
        #self.U = num.dot(Tphi, Tc.transpose())
        self.U = num.dot(Tphi, num.linalg.inv(Tc))
        
    def calc_Z(self,phi=0.0,chi=0.0,eta=0.0,mu=0.0):
        """
        Calculate the goniometer rotation matrix Z
        angles are in degrees
        """
        P = num.array([[ cosd(phi), sind(phi), 0.],
                       [-sind(phi), cosd(phi), 0.],
                       [  0.,          0.,     1.]],float)
        X = num.array([[ cosd(chi), 0., sind(chi)],
                       [   0.,      1.,    0.],
                       [-sind(chi), 0., cosd(chi)]],float)
        H = num.array([[ cosd(eta), sind(eta), 0.],
                       [-sind(eta), cosd(eta), 0.],
                       [   0.,         0.,     1.]],float)
        M  = num.array([[  1.,         0.,     0.      ],
                        [  0.,      cosd(mu), -sind(mu)],
                        [  0.,      sind(mu), cosd(mu)]],float)
        Z = num.dot(num.dot(num.dot(M,H),X),P)
        return Z

    def calc_Q(self,nu=0.0,delta=0.0,lam=None):
        """
        Calculate ki, kr, and Q in lab frame.
        Angles are in degrees, lam is in angstroms
        If lam = None, then lambda defined for the lattice
        is used.
        """
        if lam == None: lam = self.lattice.lam
        k  = (2.* num.pi / lam)
        ki = k * num.array([0.,1.,0.],dtype=float)
        kr = k * num.array([sind(delta),
                            cosd(nu)*cosd(delta),
                            sind(nu)*cosd(delta)],dtype=float)
        Q = kr - ki
        return (Q, ki, kr)

##########################################################################
##########################################################################
def test_psic():
    # create a new psic instance
    psic = Psic(5.6,5.6,13,90,90,120,lam=1.2)
    return psic


##########################################################################
##########################################################################
if __name__ == "__main__":
    """
    test 
    """
    psic = test_psic()
    
