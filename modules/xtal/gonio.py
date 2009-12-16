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
    (i.e. z is parallel to the phi axis)
Therefore with all gonio angles set to zero the lab xyz axis are coincident
with the goniometer axes in a well defined way (depending on the axis definitions).
Under the instrument settings with all angles zero (phi frame) and with the
sample oriented in an arbitrary maner the matrix U is used
to calculate the lab frame indicies of an arbritrary recip lattice vector,
h=[h,k,l], such that:
   hphi = U*hc = U*B*h
Therefore, U is a simple rotation matrix which gives
the indicies of h in the phi frame accounting for how 
the sample is mounted. ie this matrix (or its transpose)
would take the cartesian reciprocal basis vectors (hc) and
rotate them to be coincident with the lab frame (phi-frame)
basis vectors

The algorithm to determine U from a set of reflections and angles
is given by Busing and Levy

We can then orient the hphi vector in the lab frame applying the rotation
matricies of the various gonio axes.  

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
or
       |      sin(del)      |
Qm = k*|cos(del)*cos(nu) - 1|
       |  cos(del)*sin(nu)  |

Therefore the diffraction condition in the lab frame is
   Qm = (2*pi)*hm = (2*pi)*Z*hphi = (2*pi)*Z*U*B*h

Now given an orientation matrix and set of gonio
angles we can then solve for hphi
   hphi = inv(Z) * Qm / (2*pi)
The reciprocal lattice indicies (h) are then calc from
   h = inv(UB)*hphi
This gives the hkl values of the vector that is in the
diffraction condition for a given set of angles.  

1. H. You, J. Appl. Cryst. (1999) 32, 614-623
2. Busy and Leving

"""
##########################################################################

import numpy as num
import types

from mathutil import cosd, sind, tand
from mathutil import arccosd, arcsind, arctand
from lattice import Lattice

##########################################################################
def kap_angles(angles):
    """
    Angles defined by spec
    
    Note: calc kappa, keta and kphi
    kap_alp = 50.031;
    keta    = eta - asin(-tan(chi/2)/tan(kap_alp))
    kphi    = phi - asin(-tan(chi/2)/tan(kap_alp))
    kappa   = asin(sin(chi/2)/sin(kap_alp))
    """
    # angles from spec
    phi = num.radians(angles[3])
    chi = num.radians(angles[2])
    eta = num.radians(angles[1])
    mu  = num.radians(angles[5])
    delta = num.radians(angles[0])
    nu = num.radians(angles[4])

    # kappa angles
    kap_alp = 50.031;
    keta    = eta - arcsind(-tand(chi/2.)/tand(kap_alp))
    kphi    = phi - arcsind(-tand(chi/2.)/tand(kap_alp))
    kappa   = asind(sind(chi/2.)/sind(kap_alp))

    return {'phi':phi,'chi':chi,'eta':eta,'mu':mu,
            'delta':delta,'nu':nu,
            'keta':keta,'kphi':kphi,'kappa':kappa}

##########################################################################
class Psic:
    """
    Orientation calculations for Psic geometry.

    The default dummy orientation matrix is set up 
    assuming the sample is mounted such that (001) plane
    is perpendicular to the eta and phi rot axes
    (ie c-axis is parrallel to the eta and phi rot axes)
    and that the b-axis is parrallel to the nu and mu rot axes
    (ie parrallel to the lab frame Z)

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

    def calc_h(self,phi=0.0,chi=0.0,eta=0.0,mu=0.0,
               nu=0.0,delta=0.0,lam=None):
        """
        Given the gonio angles and wavelength, calculate
        the hkl values of the vector that is in the
        diffraction condition for a given set of angles.  

        Solve for hphi:
           hphi = inv(Z) * Qm / (2*pi)
        then calc h from
           h = inv(UB)*hphi
        """
        Z = self.calc_Z(phi=phi,chi=chi,eta=eta,mu=mu)
        Q = self.calc_Q(nu=nu,delta=delta,lam=lam)
        hphi = (2.*pi)*num.dot(num.linalg.inv(Z),Q) 
        UB   = num.dot(self.U,self.B)
        h    = num.dot(num.linalg.inv(UB),hphi)
        return h

    def calc_surf_norm(self,f_chi,f_phi):
        """
        Calculate the surface normal hkl given
        flat phi and flat chi settings
        """
        pass


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
    
    # diffr. angles:
    delta = 46.9587
    eta   = 6.6400
    chi   = 37.1005
    phi   = 65.78
    nu    = 0.000
    mu    = 0.0
    lam   = 1.3756
    # reference vector/surface normal in [h,k,l]
    # n = [0, 0, 1]
    n = [-0.0348357, -0.00243595, 1]

##########################################################################
##########################################################################
##########################################################################
##########################################################################
"""
function n_hkl = calc_surf_norm(f_chi, f_phi, UB);
# calc the surface normal in HKL
# given the flat phi and flat chi values

sig_az = -f_chi;
tau_az = -f_phi;

# this block converts the chi and phi values to correctly 
# defined polar coordinates, ie 0<= sig_az <= 180deg .....
if sig_az < 0;
   sig_az = -1*sig_az;
   if tau_az < 0;
      tau_az = 180 + tau_az;
   elseif tau_az > 0;
      tau_az = tau_az - 180;
   end
end

# n in the unrotated lab frame (ie phi frame):
# this is a unit vector!!!!
n_phi = [  sind(sig_az)*cosd(tau_az) ; 
          -sind(sig_az)*sind(tau_az) ; 
                        cosd(sig_az) ];

# n in HKL
n_hkl = inv(UB)*n_phi;
n_hkl = n_hkl/ max(abs(n_hkl));

# note if l-component is negative, then its
# pointing into the surface (ie asume positive L
# is the positive direction away from the surface
# careful here!!!!
if n_hkl(3) < 0;
   n_hkl = -1*n_hkl;
end;

# note to actually normalize need the cell parameters
# above is just setting the largets component to unity
# ie n_hkl isn't a unit vector!
# also need cell params and (HKL) plane to calc miscut.
# ie whats the angle between n_hkl and H_hkl

# test
#disp('#############test##############');
#n_phi = UB*n_hkl;
#n_phi = n_phi/ abs(v_mag(n_phi))
#
# note result of acosd is between 0 and pi
# get correct sign from the sign of the x-component
# sigma_az = acos(n_phi(3))*180/pi
# sigma_az = sign(n_phi(1))*acos(n_phi(3))*180/pi
# sigma_az = acos(n_phi(3))*180/pi
#
#tau_az = atan2(-n_phi(2), n_phi(1))*180/pi
#tau_az = atan(-n_phi(2)/ n_phi(1))*180/pi
"""
#
"""
###################################################
# Now calc some interesting parameters            #
###################################################

# calc tth, d and magnitude of Q
#[tth, d, q_mag] = calc_tth(h,cell,lam);

# alternate method of calc tth:
tth = acos(cosd(del)*cosd(nu))*180/pi

# calc qaz, the angle btwn Q and the yz plane 
# qaz = atan2(tand(del), sind(nu) ) *180/pi
qaz = atan2(sind(del), cosd(del)*sind(nu) ) *180/pi

# calc n in the lab frame (unrotated) and make a unit vector:
n_phi = UB*n;
n_phi = n_phi/v_mag(n_phi)

# calc n in the rotated lab frame and make a unit vector
n_l = Z*UB*n;
n_l = n_l/v_mag(n_l)

# calc the sigma_az and tau_az angles:
# sigma_az = angle between the z-axis and n
# tau_az = angle between the projection of n
#          in the xy-plane and the x-axis
sigma_az = acos(n_phi(3))*180/pi
tau_az = atan2(-n_phi(2), n_phi(1))*180/pi

# calc naz, this is the angle btwn n and the yz plane
naz = atan2( n_l(1), n_l(3) ) *180/pi

# calc alpha, ie incidence angle or angle btwn 
# k_in (which is along y) and the plane perp to n 
# (note the neg sign):
alpha = asin(dot_product(n_l,[0;-1;0]))*180/pi

# calc tau, this is the angle btwn n and the scatt-plane
# tau = acos( cosd(alpha) * cosd(tth/2) * cosd(naz - qaz) ...
#            + sind(alpha) * sind(tth/2) ) * 180/pi
tau = v_angle(Q_l,n_l)

# calc beta, ie exit angle, or angle btwn k_f and the
# plane perp to n
# beta = asin( 2*sind(tth/2)*cosd(tau) - sind(alpha) ) * 180/pi
kf_l = k*[sind(del); cosd(nu)*cosd(del); sind(nu)*cosd(del)];
beta = asin( ( dot_product(n_l, kf_l)/k) ) *180/pi

# calc psi, this is the azmuthal angle of n wrt Q, 
# ie for tau != 0, psi is the rotation of n about Q
xx = (cosd(tau)*sind(tth/2) - sind(alpha)) / ...
             (sind(tau)*cosd(tth/2));
psi = acos( xx ) * 180/pi

#xx = (-cosd(tau)*sind(tth/2) + sind(beta)) / ...
#             (sind(tau)*cosd(tth/2));
#psi = acos( xx ) * 180/pi


# calc omega, this is the angle between Q and the plane
# which is perpendicular to the axis of the chi circle.
# note for nu=mu=0 this is the same as the four circle def:
# omega = 0.5*TTH - TH, where TTH is the detector motor (=del)
# and TH is the sample circle (=eta).  Therefore, for 
# mu=nu=0 and del=0.5*eta, omega = 0, which means that Q
# is in the plane perpendicular to the chi axis.

# check the mult order here!!!!
# T = (H')*(M');
T = (M')*(H');
#T = inv(XX'); # this just gives back T
Qpp = T*Q_l;
omega = -1*v_angle([Qpp(1), 0, Qpp(3)],Qpp)

"""

