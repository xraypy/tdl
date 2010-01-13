##########################################################################
"""
Tom Trainor (tptrainor@alaska.edu)
Frank Heberling (Frank.Heberling@ine.fzk.de)

Gonio calcs for 6 circle psic geometry

Modifications:
--------------

"""
##########################################################################
"""
Todo
- Test.
- Rename this module to 'geom_psic' ?
- Psic should store ki and kr since area calcs use these
  (no need to recalc)

"""
##########################################################################
"""
Notes
------
The following describes the calculations used to define the orientation
matrix of a xtal mounted on a goniometer and to convert between
motor angles and hkl and visa-versa.  Additional functions are also
provided for computing beam and detector slit aperature vectors
and sample position vectors etc..  

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
system for the psic geometry such that:
  x is vertical (perpendicular, pointing to the ceiling of the hutch)
  y is directed along the incident beam path
  z make the system right handed and lies in the horizontal scattering plane
    (i.e. z is parallel to the phi axis)
Therefore with all gonio angles set to zero the lab xyz axis are coincident
with the goniometer axes in a well defined way (depending on the axis definitions).
Note that the yz plane is the horizontal scattering plane, and yx is the vertical
scattering plane.

With the instrument settings all at zero angles (phi frame), and with the
sample oriented in an arbitrary maner, the matrix U is used to calculate
the lab frame indicies (hphi) of a given reciprocal lattice vector, h=[h,k,l],
according to:
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
import copy

from mathutil import cosd, sind, tand
from mathutil import arccosd, arcsind, arctand
from mathutil import cartesian_mag, cartesian_angle

from lattice import Lattice

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
    ###################################################
    def __init__(self,a=10.,b=10.,c=10.,alpha=90.,beta=90.,gamma=90.,lam=1.0):
        """
        Initialize by passing a,b,c in angstroms 
        alpha, beta, gamma in degrees,
        and lambda in angstroms
        """
        # set lattice and lambda
        self.lattice = Lattice(a,b,c,alpha,beta,gamma,lam)
        
        # hold gonio angles 
        self.angles={'phi':0.0,'chi':0.0,'eta':0.0,'mu':0.0,
                     'nu':0.0,'delta':0.0}

        # hold psuedo angles
        self.pangles = {}
        self.calc_psuedo = True

        # hold n (reference) vector in HKL
        # eg surface normal vector for psuedo angles
        self.n = num.array([0.,0.,1.],dtype=float)
        
        # Z and calc h
        self.Z = []
        self.Q = []
        self.h = [0.,0.,0.]

        # dummy primary reflection
        tth = self.lattice.tth([0.,0.,1.],lam=lam)
        self.or0={'h':num.array([0.,0.,1.]),
                  'phi':0.0,'chi':0.0,'eta':0.0,'mu':tth/2.,
                  'nu':tth,'delta':0.0,'lam':lam}
        
        # dummy secondary reflection
        tth = self.lattice.tth([0.,1.,0.],lam=lam)
        self.or1={'h':num.array([0.,1.,0.]),
                  'phi':0.0,'chi':0.0,'eta':tth/2.,'mu':0.0,
                  'nu':0.0,'delta':tth,'lam':lam}

        # Compute OR matricies
        self.U = []
        self.B = []
        self.UB = []
        self._calc_UB()

    ###################################################
    def __repr__(self,):
        lout = self.lattice.__repr__()
        lout = "%sPrimary:\n   h=%3.2f,k=%3.2f," % (lout,self.or0['h'][0],self.or0['h'][1])
        lout = "%sl=%3.2f, lam=%6.6f\n" % (lout,self.or0['h'][2],self.or0['lam'])
        lout = "%s   phi=%6.3f,chi=%6.3f," % (lout,self.or0['phi'],self.or0['chi'])
        lout = "%seta=%6.3f,mu=%6.3f," % (lout,self.or0['eta'],self.or0['mu'])
        lout = "%snu=%6.3f,delta=%6.3f\n" % (lout,self.or0['nu'],self.or0['delta'])
        #
        lout = "%sSecondary:\n   h=%3.2f,k=%3.2f," % (lout,self.or1['h'][0],self.or1['h'][1])
        lout = "%sl=%3.2f, lam=%6.6f\n" % (lout,self.or1['h'][2],self.or1['lam'])
        lout = "%s   phi=%6.3f,chi=%6.3f," % (lout,self.or1['phi'],self.or1['chi'])
        lout = "%seta=%6.3f,mu=%6.3f," % (lout,self.or1['eta'],self.or1['mu'])
        lout = "%snu=%6.3f,delta=%6.3f\n" % (lout,self.or1['nu'],self.or1['delta'])
        #
        lout = "%sSetting:" % (lout)
        lout = "%s   h=%3.2f,k=%3.2f,l=%3.2f\n" % (lout,self.h[0],self.h[1],self.h[2])
        lout = "%s   phi=%6.3f,chi=%6.3f," % (lout,self.angles['phi'],self.angles['chi'])
        lout = "%seta=%6.3f,mu=%6.3f," % (lout,self.angles['eta'],self.angles['mu'])
        lout = "%snu=%6.3f,delta=%6.3f\n" % (lout,self.angles['nu'],self.angles['delta'])
        #
        if self.calc_psuedo:
            lout = "%s   TTH=%6.3f," % (lout,self.pangles['tth'])
            lout = "%sSIGMA_AZ=%6.3f," % (lout,self.pangles['sigma_az'])
            lout = "%sTAU_AZ=%6.3f," % (lout,self.pangles['tau_az'])
            lout = "%sN_AZ=%6.3f," % (lout,self.pangles['naz'])
            lout = "%sALPHA=%6.3f," % (lout,self.pangles['alpha'])
            lout = "%sBETA=%6.3f\n" % (lout,self.pangles['beta'])
            lout = "%s   TAU=%6.3f," % (lout,self.pangles['tau'])
            lout = "%sPSI=%6.3f," % (lout,self.pangles['psi'])
            lout = "%sQ_AZ=%6.3f," % (lout,self.pangles['qaz'])
            lout = "%sOMEGA=%6.3f," % (lout,self.pangles['omega'])
        #
        return lout
    
    ###################################################
    def set_lat(self,a=None,b=None,c=None,alpha=None,
                beta=None,gamma=None,lam=None):
        """
        Update lattice parameters and lambda

        abc are in angstroms, angles are in degrees,
        lam is in angstroms

        """
        self.lattice.update(a=a,b=b,c=c,alpha=alpha,
                            beta=beta,gamma=gamma,lam=lam)
        self._calc_UB()

    def set_spec_G(self,G):
        """
        Take the spec G array for the psic geometry
        and set all the relevant orientation info...
        """
        (cell,or0,or1,n) = spec_psic_G(G)
        self.n   = n
        self.or0 = or0
        self.or1 = or1
        self.lattice = Lattice(*cell)
        self._calc_UB()

    ################################################### 
    def set_or0(self,h=None,phi=None,chi=None,eta=None,
                mu=None,nu=None,delta=None,lam=None):
        """
        Set / adjust the primary orientation reflection

        Angles are in degrees, lam is in angstroms
        If lam = None, then lambda defined for the lattice
        is used.
        """
        if h!=None:     self.or0['h'] = num.array(h,dtype=float)
        if phi!=None:   self.or0['phi']=float(phi)
        if chi!=None:   self.or0['chi']=float(chi)
        if eta!=None:   self.or0['eta']=float(eta)
        if mu!=None:    self.or0['mu']=float(mu)
        if nu!=None:    self.or0['nu']=float(nu)
        if delta!=None: self.or0['delta']=float(delta)
        if lam!= None:  self.or0['lam']=float(lam)
        self._calc_UB()

    def set_or1(self,h=None,phi=None,chi=None,eta=None,
                mu=None,nu=None,delta=None,lam=None):
        """
        Set / adjust the secondary orientation reflection

        Angles are in degrees, lam is in angstroms
        If lam = None, then lambda defined for the lattice
        is used.
        """
        if h!=None:     self.or1['h'] = num.array(h,dtype=float)
        if phi!=None:   self.or1['phi']=float(phi)
        if chi!=None:   self.or1['chi']=float(chi)
        if eta!=None:   self.or1['eta']=float(eta)
        if mu!=None:    self.or1['mu']=float(mu)
        if nu!=None:    self.or1['nu']=float(nu)
        if delta!=None: self.or1['delta']=float(delta)
        if lam!= None:  self.or1['lam']=float(lam)
        self._calc_UB()

    def swap_or(self,):
        """
        Swap the primary and secondary reflection
        """
        tmp = copy.copy(self.or0)
        self.or0 = copy.copy(self.or1)
        self.or1 = tmp
        self._calc_UB()

    ################################################### 
    def _calc_UB(self,):
        """
        Calculate the orientation matrix, U,
        from the primary and secondary
        reflectons and given lattice

        Note dont really ever use B by itself.  so we
        should combine this and above to calc_UB and
        just store UB??
        
        """
        # use these, note they are used below on vectors
        # defined in the cartesian lab frame basis
        cross = num.cross
        norm  = num.linalg.norm

        #Calculate the B matrix
        (a,b,c,alp,bet,gam)       = self.lattice.cell()
        (ar,br,cr,alpr,betr,gamr) = self.lattice.rcell()
        B = num.array([[ar,  br*cosd(gamr),     cr*cosd(betr)        ],
                       [0.,  br*sind(gamr),  -cr*sind(betr)*cosd(alp)],
                       [0.,      0.,                 1./c            ]])
        self.B = B
        
        # calc Z and Q for the OR reflections
        Z1 = calc_Z(self.or0['phi'],self.or0['chi'],self.or0['eta'],self.or0['mu'])
        Q1 = calc_Q(self.or0['nu'],self.or0['delta'],self.or0['lam'])
        #
        Z2 = calc_Z(self.or1['phi'],self.or1['chi'],self.or1['eta'],self.or1['mu'])
        Q2 = calc_Q(self.or1['nu'],self.or1['delta'],self.or1['lam'])

        # calc the phi frame coords for diffraction vectors
        # note divide out 2pi since the diffraction is 2pi*h = Q
        vphi_1 = num.dot(num.linalg.inv(Z1), (Q1/(2.*num.pi)))
        vphi_2 = num.dot(num.linalg.inv(Z2), (Q2/(2.*num.pi)))
        
        #calc cartesian coords of h vectors
        hc_1 = num.dot(self.B, self.or0['h'])
        hc_2 = num.dot(self.B, self.or1['h'])

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

        # calc UB
        self.UB = num.dot(self.U,self.B)

        #update h and psuedo angles...
        self.set_angles()

    ###################################################
    def set_angles(self,phi=None,chi=None,eta=None,
                   mu=None,nu=None,delta=None):
        """
        Set goniometer angles (all in degrees)
        """
        if phi!=None:   self.angles['phi']=float(phi)
        if chi!=None:   self.angles['chi']=float(chi)
        if eta!=None:   self.angles['eta']=float(eta)
        if mu!=None:    self.angles['mu']=float(mu)
        if nu!=None:    self.angles['nu']=float(nu)
        if delta!=None: self.angles['delta']=float(delta)
        # update h, also calc Z etc..
        self._calc_h()
        # update psuedo
        self._update_psuedo()
        
    def _calc_h(self,):
        """
        Calculate the hkl values of the vector that is in the
        diffraction condition for the given set of angles.  

        Solve for hphi using Z and lab frame Q:
           hphi = inv(Z) * Q / (2*pi)
        then calc h from
           h = inv(UB)*hphi
        """
        self.Z = calc_Z(phi=self.angles['phi'],chi=self.angles['chi'],
                        eta=self.angles['eta'],mu=self.angles['mu'])
        self.Q = calc_Q(self.angles['nu'],self.angles['delta'],self.lattice.lam)
        
        hphi = num.dot(num.linalg.inv(self.Z),self.Q) / (2.*num.pi) 
        h    = num.dot(num.linalg.inv(self.UB),hphi)
        self.h = h
        
    ###################################################
    def set_n(self,n=[0,0,1]):
        """
        Set n, the reference vector used for psuedo angles.
        The n vector is given in hkl values.  see calc_n
        to determine n from chi and phi settings
        """
        self.n = num.array(n,dtype=float)
        self._update_psuedo()

    def calc_n(self,fchi=0.0,fphi=0.0):
        """
        Calculate the hkl values of a reference vector given
        the chi and phi settings that align this
        vector with the eta axis.

        For example this algorith is used
        to compute the surface normal from the (flat) chi and
        (flat) phi angles that leave an optical reflection in
        a fixed position during an eta rotation 

        Note the vector is normalized such that
        the largest component is unity,
        ie n_hkl isn't a unit vector!
        
        """
        # polar angles
        sig_az = -fchi
        tau_az = -fphi

        # this block converts the chi and phi values to correctly 
        # defined polar coordinates, ie 0<= sig_az <= 180deg .....
        if sig_az < 0.:
            sig_az = -1.*sig_az
            if tau_az < 0.:
                tau_az = 180. + tau_az
            elif tau_az > 0.:
                tau_az = tau_az - 180.

        # n in the unrotated lab frame (ie phi frame):
        # this is a unit vector!
        n_phi = num.array([ sind(sig_az)*cosd(tau_az),
                           -sind(sig_az)*sind(tau_az), 
                                  cosd(sig_az)        ], dtype=float)
        # n in HKL
        n_hkl = num.dot(num.linalg.inv(self.UB),n_phi)
        n_hkl = n_hkl/ num.max(num.abs(n_hkl))
        
        # note if l-component is negative, then its
        # pointing into the surface (ie assume positive L
        # is the positive direction away from the surface
        # careful here!!
        if n_hkl[2] < 0.:
            n_hkl = -1.*n_hkl

        # set n which triggers recalc of
        # all the psuedo angles
        self.set_n(n_hkl)

    ################################################### 
    ## Pseudo angles
    ###################################################
    def _update_psuedo(self):
        """
        Compute psuedo angles
        Note use this to compute psuedo angles rather than
        individual calls.  ie some psuedo angles depend on others
        so its important that the calc are executed in the correct
        order.  Also important is that _calc_h is called before this...
        """
        if self.calc_psuedo == True:
            self._calc_tth()
            self._calc_nm()
            self._calc_sigma_az()
            self._calc_tau_az()
            self._calc_naz()
            self._calc_alpha()
            self._calc_beta()
            self._calc_tau()
            self._calc_psi()
            self._calc_qaz()
            self._calc_omega()
        else:
            self.pangles = {}
    
    def _calc_tth(self):
        """
        Calculate 2Theta, the scattering angle
        
        This should be the same as:
          (ki,kr) = self._calc_kvecs()
           tth = cartesian_angle(ki,kr)

        You can also get this given h, the reciprocal lattice
        vector that is in the diffraction condition.  E.g.
          h   = self.calc_h()
          tth = self.lattice.tth(h)
        """
        nu    = self.angles['nu']
        delta = self.angles['delta']
        tth   = arccosd(cosd(delta)*cosd(nu))
        self.pangles['tth'] = tth

    def _calc_nm(self):
        """
        Calculate the rotated cartesian lab indicies
        of the reference vector n

        The reference vector n is given in recip
        lattice indicies (hkl)
        """
        # calc n in the rotated lab frame and make a unit vector
        n  = self.n
        Z  = self.Z
        UB = self.UB
        nm = num.dot(num.dot(Z,UB),n)
        nm = nm/cartesian_mag(nm)
        self.nm = nm
    
    def _calc_sigma_az(self):
        """
        sigma_az = angle between the z-axis and n in the phi frame
        
        The reference vector is given in recip lattice indicies (hkl) 
        """
        # calc n in the lab frame (unrotated) and make a unit vector
        n_phi = num.dot(self.UB,self.n)
        n_phi = n_phi/cartesian_mag(n_phi)
        
        # note result of acosd is between 0 and pi
        # get correct sign from the sign of the x-component
        #sigma_az = num.sign(n_phi[0])*arccosd(n_phi[2])
        sigma_az = arccosd(n_phi[2])
        self.pangles['sigma_az'] = sigma_az

    def _calc_tau_az(self):
        """
        tau_az = angle between the projection of n in the
        xy-plane and the x-axis

        The reference vector is given in recip lattice indicies (hkl) 
        """
        n_phi = num.dot(self.UB,self.n)
        n_phi = n_phi/cartesian_mag(n_phi)
        tau_az = num.arctan2(-n_phi[1], n_phi[0])
        tau_az = tau_az*180./num.pi
        self.pangles['tau_az'] = tau_az

    def _calc_naz(self):
        """
        calc naz, this is the angle btwn the reference vector n 
        and the yz plane at the given angle settings

        The reference vector is given in recip lattice indicies (hkl) 
        """
        # get norm reference vector in cartesian lab frame
        nm  = self.nm
        naz = num.arctan2( nm[0], nm[2] )
        naz = num.degrees(naz)
        self.pangles['naz'] = naz

    def _calc_alpha(self):
        """
        Calc alpha, ie incidence angle or angle btwn 
        k_in (which is along y) and the plane perp to
        the reference vector n.

        The reference vector is given in recip lattice indicies (hkl) 
        """
        # get norm reference vector in cartesian lab frame
        nm = self.nm
        # note the neg sign
        ki = num.array([0.,-1.,0.],dtype=float)
        alpha = arcsind(num.dot(nm,ki))
        self.pangles['alpha'] = alpha

    def _calc_beta(self):
        """
        Calc beta, ie exit angle, or angle btwn k_r and the
        # plane perp to the reference vector n
        """
        # get norm reference vector in cartesian lab frame
        nm = self.nm
        # calc normalized kr
        delta = self.angles['delta']
        nu    = self.angles['nu']
        kr = num.array([sind(delta),
                        cosd(nu)*cosd(delta),
                        sind(nu)*cosd(delta)],dtype=float)
        beta = arcsind(num.dot(nm, kr))
        # beta = arcsind( 2*sind(tth/2)*cosd(tau) - sind(alpha) )
        self.pangles['beta'] = beta

    def _calc_tau(self):
        """
        Calc tau, this is the angle btwn n and the scattering-plane
        defined by ki and kr.  ie the angle between n and Q

        Can also calc from:
         tau = acos( cosd(alpha) * cosd(tth/2) * cosd(naz - qaz) ...
                    + sind(alpha) * sind(tth/2) ) 
        """
        # get norm reference vector in cartesian lab frame
        nm  = self.nm
        tau = cartesian_angle(self.Q,nm)
        self.pangles['tau'] = tau

    def _calc_psi(self):
        """
        calc psi, this is the azmuthal angle of n wrt Q. 
        ie for tau != 0, psi is the rotation of n about Q

        Note this must be calc after tth,tau, and alpha!
        """
        tau   = self.pangles['tau']
        alpha = self.pangles['alpha']
        tth   = self.pangles['tth']
        xx    = (cosd(tau)*sind(tth/2.) - sind(alpha))
        xx    = xx /(sind(tau)*cosd(tth/2.))
        psi = arccosd( xx )
        #beta = self.calc_beta()
        #xx = (-cosd(tau)*sind(tth/2.) + sind(beta))
        # xx = xx /(sind(tau)*cosd(tth/2.))
        #psi = arccosd( xx )
        self.pangles['psi'] = psi

    def _calc_qaz(self):
        """
        Calc qaz, the angle btwn Q and the yz plane 
        """
        nu    = self.angles['nu']
        delta = self.angles['delta']
        qaz = num.arctan2(sind(delta), cosd(delta)*sind(nu) )
        qaz = num.degrees(qaz)
        self.pangles['qaz'] = qaz

    def _calc_omega(self):
        """
        calc omega, this is the angle between Q and the plane
        which is perpendicular to the axis of the chi circle.
        note for nu=mu=0 this is the same as the four circle def:
        omega = 0.5*TTH - TH, where TTH is the detector motor (=del)
        and TH is the sample circle (=eta).  Therefore, for 
        mu=nu=0 and del=0.5*eta, omega = 0, which means that Q
        is in the plane perpendicular to the chi axis.

        Note check sign of results??? 
        """
        phi=self.angles['phi']
        chi=self.angles['chi']
        eta=self.angles['eta']
        mu=self.angles['mu']
        H = num.array([[ cosd(eta), sind(eta), 0.],
                       [-sind(eta), cosd(eta), 0.],
                       [   0.,         0.,     1.]],float)
        M  = num.array([[  1.,         0.,     0.      ],
                        [  0.,      cosd(mu), -sind(mu)],
                        [  0.,      sind(mu), cosd(mu)]],float)
        # check the mult order here!!!!
        # T = num.dot(H.transpose(),M.transpose())
        T     = num.dot(M.transpose(),H.transpose())
        Qpp   = num.dot(T,self.Q)
        #omega = -1.*cartesian_angle([Qpp[0], 0, Qpp[2]],Qpp)
        omega = cartesian_angle([Qpp[0], 0, Qpp[2]],Qpp)
        self.pangles['omega'] = omega

##########################################################################
def psic_from_spec(G,angles={}):
    """
    pass spec G array and dictionary of angles
    returns a psic instance
    """
    gonio = Psic()
    if G != None: gonio.set_spec_G(G)
    gonio.set_angles(**angles)
    return gonio

##########################################################################
def spec_psic_G(G):
    """
    Parse essential lattice and OR data
    from the spec G array for psic geometry
    See specfile.py for details.
    
    """
    #azimuthal reference vector, n (hkl)
    n = num.array(G[3:6],dtype=float)

    #lattice params a,b,c,alp,bet,gam
    cell = G[22:28]
    # add lambda to end of cell
    cell.append(G[66])
    cell = num.array(cell,dtype=float)

    # or0
    or0 = {}
    or0['h']   = num.array(G[34:37],dtype=float)
    or0.update(_spec_or_angles(num.array(G[40:46],dtype=float)))
    or0['lam'] = float(G[52])

    # or1
    or1 = {}
    or1['h']   = num.array(G[37:40],dtype=float)
    or1.update(_spec_or_angles(num.array(G[46:52],dtype=float)))
    or1['lam'] = float(G[53])

    return (cell,or0,or1,n)

##########################################################################
def _spec_or_angles(angles,calc_kappa=False):
    """
    Angles defined by spec for the OR's.
    See specfile.py

    Assume the following.
    
    If parsing angles from the P array:
    (generally shouldnt need this since angles
    are tagged with motor labels on read)
    angles = P
    if psic:
       angles = angles[0:5]
    elif kappa fourc
       angles = [angles[0:3], angles[8], angles[7]]

    If parsing angles from the G array:
      angles = G[x:y]
    where x:y depend on whether you are
    parsing out or0 or or1.
    See spec_G below

    We then assume the following:
      del = angles[0]
      eta = angles[1]
      chi = angles[2]
      phi = angles[3]
      nu  = angles[4]
      mu  = angles[5]

    Note: calc kappa, keta and kphi
    kap_alp = 50.031;
    keta    = eta - asin(-tan(chi/2)/tan(kap_alp))
    kphi    = phi - asin(-tan(chi/2)/tan(kap_alp))
    kappa   = asin(sin(chi/2)/sin(kap_alp))
    """
    # angles from spec
    delta = angles[0]
    eta   = angles[1]
    chi   = angles[2]
    phi   = angles[3]
    nu    = angles[4]
    mu    = angles[5]

    # kappa angles
    if calc_kappa:
        kap_alp = 50.031;
        keta    = eta - arcsind(-tand(chi/2.)/tand(kap_alp))
        kphi    = phi - arcsind(-tand(chi/2.)/tand(kap_alp))
        kappa   = asind(sind(chi/2.)/sind(kap_alp))
        return {'phi':phi,'chi':chi,'eta':eta,'mu':mu,
                'delta':delta,'nu':nu,
                'keta':keta,'kphi':kphi,'kappa':kappa}
    else:
        return {'phi':phi,'chi':chi,'eta':eta,'mu':mu,
                'delta':delta,'nu':nu}

##########################################################################
def calc_Z(phi=0.0,chi=0.0,eta=0.0,mu=0.0):
    """
    Calculate the psic goniometer rotation matrix Z
    for the 4 sample angles. Angles are in degrees

    Z is the matrix that rotates a vector defined in the phi frame
    ie a vector defined with all angles zero => vphi.  After rotation
    the lab frame coordinates of the vector => vm are given by:
         vm = Z*vphi
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

##########################################################################
def calc_Q(nu=0.0,delta=0.0,lam=1.0):
    """
    Calculate psic Q in the cartesian lab frame.
    nu and delta are in degrees, lam is in angstroms
    """
    (ki,kr) = calc_kvecs(nu=nu,delta=delta,lam=lam)
    Q = kr - ki
    return Q

##########################################################################
def calc_kvecs(nu=0.0,delta=0.0,lam=1.0):
    """
    Calculate psic ki, kr in the cartesian lab frame.
    nu and delta are in degrees, lam is in angstroms
    """
    k  = (2.* num.pi / lam)
    ki = k * num.array([0.,1.,0.],dtype=float)
    kr = k * num.array([sind(delta),
                        cosd(nu)*cosd(delta),
                        sind(nu)*cosd(delta)],dtype=float)
    return (ki,kr)

##########################################################################
def calc_D(nu=0.0,delta=0.0):
    """
    Calculate the detector rotation matrix.
    Angles are in degrees

    D is the matrix that rotates a vector defined in the phi frame
    ie a vector defined with all angles zero => vphi.  After rotation
    the lab frame coordinates of the vector => vm are given by:
         vm = D*vphi
    For example 
                            |0|   
         kr_phi = (2pi/lam) |1|
                            |0|
    Since kr is defined by the detector rotation, the lab frame 
    coordinates of the kr vector after detector rotation are
         kr_m = D*kr_phi
    
    """
    D1 = num.array([[cosd(delta),  sind(delta),  0.], 
                    [-sind(delta), cosd(delta),  0.],
                    [     0.     ,     0.     ,  1.]])
          
    D2 = num.array([[    1.,     0.   ,      0.  ],
                    [    0.,  cosd(nu), -sind(nu)], 
                    [    0.,  sind(nu),  cosd(nu)]])
          
    D = num.dot(D2,D1)
    return (D)

##########################################################################
def beam_vectors(w=1.0,h=1.0):
    """
    Compute the beam apperature vectors in lab frame
    
    The slit settings:
       w = beam horz width (total slit width in lab-z)
       h = beam vert hieght (total slit width in lab-x)

    Assume these are centered on the origin
    """
    # beam vectors, [x,y,z], in lab frame
    bh = num.array([   0.,  0., 0.5*w])
    bv = num.array([0.5*h,  0.,    0.])

    # corners of beam apperature
    a =  bv + bh
    b =  bv - bh
    c = -bv - bh
    d = -bv + bh
    beam = [a,b,c,d]

    return beam

##########################################################################
def det_vectors(w=1.0,h=1.0,nu=0.0,delta=0.0):
    """
    Compute detector apperature vectors in lab frame
    
    The slit settings (all in same units, eg mm):
      w = detector horz width (total slit width in lab-z)
      h = detector vert hieght (total slit width in lab-x)

    Assume these are centered on the origin

    """
    # detector vectors, [x,y,z] in lab frame
    # note rotation of the vectors...
    dh = num.array([   0.,  0.,  0.5*w])
    dv = num.array([0.5*h,  0.,  0.   ]) 
    D  = calc_D(nu=nu,delta=delta)
    dh = num.dot(D,dh)
    dv = num.dot(D,dv)

    # corners of detector apperature 
    e =  dv + dh
    f =  dv - dh
    g = -dv - dh
    h = -dv + dh
    det = [e,f,g,h]

    return det

##########################################################################
def sample_vectors(sample,angles={},gonio=None):
    """
    sample = [[x,y,z],[x,y,z],[x,y,z],....]
             is a list of vectors that describe the shape of
             the sample.  They should be given in general lab
             frame coordinates.

    angles = {'phi':0.,'chi':0.,'eta':0.,'mu':0.}
             are the instrument angles at which the sample
             vectors were determined.
    
    The lab frame coordinate systems is defined such that:
        x is vertical (perpendicular, pointing to the ceiling of the hutch)
        y is directed along the incident beam path
        z make the system right handed and lies in the horizontal scattering plane
          (i.e. z is parallel to the phi axis)

    The center (0,0,0) of the lab frame is the rotation center of the instrument.

    If the sample vectors are given at the flat phi and chi values and with
    the correct sample hieght (sample Z set so the sample surface is on the
    rotation center), then the z values of the sample vectors will be zero.
    If 2D vectors are passed we therefore assume these are [x,y,0].  If this
    is the case then make sure:
        angles = {'phi':flatphi,'chi':flatchi,'eta':0.,'mu':0.}

    Note that the sample_poly that is returned is a list of 3D vectors.
    If gonio == None these are defined in the lab phi frame.
    If a gonio instance is passed then they will be rotated to the m-frame

    The easiest way to determine the sample coordinate vectors is to take a picture
    of the sample with a camera mounted such that is looks directly down the omega
    axis and the gonio angles set at the sample flat phi and chi values and
    eta = mu = 0. Then find the sample rotation center and measure the position
    of each corner (in mm) with up being the +x direction, and downstream
    being the +y direction.  

    """
    if len(sample) < 3:
        print "Sample polygon must be 3 or more points"
        return None
    # If angles are provided then we need to compute the phi
    # frame vectors first. i.e. assume p's are defined in lab
    # frame at the given set of angles. therefore we need to 
    # unrotate to get back to phi frame vectors
    # Otherwise we assume the vectors passed are phi frame
    # (in that case they should be 3D vectors)
    if angles == None: angles = {}
    if len(angles) > 0:
        # calc sample rotation matrix
        Z = calc_Z(**angles)
        Zinv = num.linalg.inv(Z)
        polygon_phi = []
        # If p's have anly two components then we assume they are 
        # xy pairs therefore we can add a third value of zero for z
        for p in sample:
            if len(p) == 2: p = [p[0],p[1],0.]
            p_phi = num.dot(Zinv,p)
            polygon_phi.append(p_phi)
    else:
        polygon_phi = sample

    # If gonio is passed then rotate the
    # vectors into the m-frame.  Otherwise
    # just return the phi frame vectors
    polygon = []
    if gonio != None:
        for p in polygon_phi:
            p_m = num.dot(gonio.Z,p)
            polygon.append(p_m)
    else:
        polygon = polygon_phi

    return polygon

##########################################################################
##########################################################################
def test1():
    # create a new psic instance
    psic = Psic(5.6,5.6,13,90,90,120,lam=1.3756)

    # diffr. angles:
    psic.set_angles(phi=65.78,chi=37.1005,eta=6.6400,
                    mu=0.0,nu=0.0,delta=46.9587)
    #calc tth, d and magnitude of Q
    print "\nh=",  psic.h
    print "tth =", psic.pangles['tth']
    print "tth =", psic.lattice.tth(psic.h)

    # reference vector/surface normal in [h,k,l]
    # n = [0, 0, 1]
    n = [-0.0348357, -0.00243595, 1]
    psic.set_n(n)
    #calc miscut
    print "\nmiscut=",psic.lattice.angle([0,0,1],n,recip=True)

    # test surf norm
    n_phi = num.dot(psic.UB,n)
    n_phi = n_phi/ num.fabs(cartesian_mag(n_phi))
    print "\nnphi: ", n_phi

    print "\nsigma_az=",psic.pangles['sigma_az']

    print "\ntau_az=",psic.pangles['tau_az']

    psic.calc_n(-psic.pangles['sigma_az'],-psic.pangles['tau_az'])
    print "calc n from sigma and tau az: ",psic.n 
    
    #calc miscut
    print "\nmiscut=",psic.lattice.angle([0,0,1],psic.n,recip=True)
    
    return psic

##########################################################################
def test2(show=True):
    """
    G parsing from specfile.py
    --------------------------
    67. ALPHA=0     # Q[4], Angle between the reference vector and the xz-plane (incidence angle) 
    68. BETA=0      # Q[5], Exit angle when the reference vector is the surface normal
    69. OMEGA=0     # Q[6], Angle btwn Q and the plane perpendicular to the chi axis
    70. TTH=0       # Q[7], Scattering angle
    71. PSI=0       # Q[8], Azimuthal angle of the reference vector wrt Q and the scatt-plane
    72. TAU=0       # Q[9], Longitudinal angle of the reference vector wrt Q and the scatt-plane
    73. QAZ=0       # Q[10], Angle btwn Q and the yz-plane
    74. NAZ=0       # Q[11], Angle btwn the reference vector and the yz-plane
    75. SIGMA_AZ=0  # Q[12], Angle to specify reference vector, = -flat_chi
    76. TAU_AZ=0    # Q[13], Angle to specify reference vector, = -flat_phi
    """
    psic = Psic()
    G = [0.0, 0.0, 1.0, 0.0039915744589999998, 0.00075650941450000001, 1.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 50.0, 0.0, 0.0, 1.0, 4.0, 4.0, 5.0, 4.0, 0.0, 0.0, 8.0939999999999994,
         4.9880000000000004, 6.0709999999999997, 90.0, 90.0, 90.0, 0.77627690969999996,
         1.2596602459999999, 1.034950635, 90.0, 90.0, 90.0, 0.0, 0.0, 4.0, 2.0, 0.0,
         2.4809999999999999, -0.00089999999999999998, 0.00080000000000000004, -0.1244,
         175.2192, 31.965, 16.27375, 15.090299999999999, 7.5477999999999996, 11.9002,
         -96.048199999999994, 17.594249999999999, 0.35625000000000001, 0.83580100000000002,
         0.83580100000000002, -0.23926186720000001, 1.1983306439999999, -0.0026481987470000001,
         -0.73847260000000003, -0.38826052760000002, -0.0050577973450000001, -0.004221190899,
         0.001168841303, 1.034934888, -0.00011346681140000001, 0.0002801762336, 5.9400141580000003,
         0.83580100000000002, 23.955432519999999, 24.314067479999999, 0.28474999779999999,
         48.269500000000001, 0.075104988760000005, 0.17931763710000001, -0.00040199175790000001,
         -0.00014478299110000001, 0.48309999999999997, 108.00069999999999, 2.0, 0.0, 0.0, 0.0, 0.0,
         12.0, 0.0, 0.0, 2.0802999999999998, 123.1461, 0.0, 0.0, 0.0, 0.0, -180.0, -180.0, -180.0,
         -180.0, -180.0, -180.0, -180.0, -180.0, -180.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0]
    psic.set_spec_G(G)

    # diffr. angles:
    psic.set_angles(phi=178.1354,chi=-0.1344,eta=0.0002,
                    mu=24.4195,nu=48.2695,delta=-0.0003)

    if show:
        print "#########################################"
        print "h calc is:", psic.h
        print "should be: [-0.000113467 0.000280176 5.94001]"
        print "----"
        print "alpha calc is:", psic.pangles['alpha']
        print "should be:", G[67]
        print "----"
        print "beta calc is:", psic.pangles['beta']
        print "should be:", G[68]
        print "----"
        print "omega calc is:", psic.pangles['omega']
        print "should be:", G[69]
        print "----"
        print "tth calc is:", psic.pangles['tth']
        print "should be:", G[70]
        print "----"
        print "psi calc is:", psic.pangles['psi']
        print "should be:", G[71]
        print "----"
        print "tau calc is:", psic.pangles['tau']
        print "should be:", G[72]
        print "----"
        print "qaz calc is:", psic.pangles['qaz']
        print "should be:", G[73]
        print "----"
        print "naz calc is:", psic.pangles['naz']
        print "should be:", G[74]
        print "----"
        print "sigma_az calc is:", psic.pangles['sigma_az']
        print "should be:", G[75]
        print "----"
        print "tau_az calc is:", psic.pangles['tau_az']
        print "should be:", G[76]
        print "#########################################"

    return psic
    
##########################################################################
if __name__ == "__main__":
    """
    test 
    """
    psic = test1()
    psic = test2()
    

