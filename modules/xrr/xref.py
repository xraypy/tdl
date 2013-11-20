"""
Calculate reflectivity / reflection XSW

Authors/Modifications:
----------------------
* T. Trainor (tptrainor@alaska.edu)

Notes:
------
For this module to work correctly you need to ensure that:
* Ifeffit is installed on the system, and the Ifeffit/bin
  directory is in the systems search path
* The tdl/lib directory must also be in the systems path
  so that gsl dll's can be found (on windows)

Notes on the calculation:
-------------------------
This model calculates the electric field intensity at any point
within a multilayer structure as a function of incidence angle, and
the Reflectivity and fluorescent yield profile.

We assume that the multilayer has 'nlayers' and each layer
within the model has a homogeneous composition and density.
Hence composition and density variations are modeled by
discrete layers with different properties.

The functions/data structures assume j = 'nlayers - 1' is the top layer,
which should generally be modelled as air/vaccum, and j = 0 is the substrate.

For e.g. a 4 layer model (nlayers=4) would have:


               Air
          j = nlayers -1 = 3
 -------------------------------- interface 2
          j = nlayers -2 = 2
 -------------------------------- interface 1
          j = nlayers -3 = 1
 -------------------------------- interface 0
          j = nlayers -4 = 0
             Substrate
                 |
                 |
                \ /
               -inf
 ---------------------------------

We use the following conventions:

The layer thickness are in angstroms and the
array is indexed as:
    d[3] = Air layer thickness (Layer 3)
    d[2] = Layer 2 thickness
    d[1] = Layer 1 thickness
    d[0] = Substrate thickness (Layer 0)

Note that the thickness of the Air and Substrate layers are arbitraty
for the calculation of the field intensity and Reflectivity.
The specified layer thicknesses for these layers currently determine how
deep in the layer (ie distance away from the interface) to calculate the
field intensity for fluorescent yield calculations.  Note if the calc
parameter pdepth is passed, d[0] is also ignored in the FY calcs,
rather the depth is calculated as a multiple of the penetration depth...

The densities are in g/cm^3:
    rho[3] = Air layer density (Layer 3)
    rho[2] = Layer 2 density
    rho[1] = Layer 1 density
    rho[0] = Substrate density (Layer 0)

The compositions are given as the elemental stiochiometric
coefficients, or mole fractions.  The are given as:
    comp[k][j] = fraction of element k in layer j

The element Z's are passed in the array (len = nelem)
    elem_z[k] = Z of element k

We also require that fp and fpp are passed in
(which we assume are calculated for the incident energy)
These arrays are of len = nelem:
    fp[k]  = f prime (electrons/atom) of element k
    fpp[k] = f double prime (electrons/atom) of element k

The amu array (len=nelem) are the elements atomic masses:
    amu[k] = mass of element k (g/mole)

The mu_at array (len=nelem) are the atomic mass absorption coefficients
(cm^2/g) at the fluorescence energy (ie used to calc the fluorescence attenuation)
    mu_at[k] = mass absorption coefficient of element k (cm^2/atom)

Note these are related to the atomic cross sections and scattering factors via:
    sig_at = 2*r_e*lambda*(fpp+coh+incoh) --> cm^2/atom
    mu_at = (Na/A)*sig_at

Therefore for a given layer (density rho (g/cm^3)) and formula weight (MW):
    mu_t =  10^-8 *(rho/MW) * Sum( A[k] * mu_at[k] * x[k])
    att  = exp(-mu_t*thickness)

The calculations also allow the inclusion of interface roughness through
a simple Debye-Waller type model.  Note that in the above model there are
3 interfaces.  Therefore the 'sigma' array should have three values, each
being the rms roughness value in angstroms:
    sigma[0] = layer0/layer1 interface roughness
    sigma[1] = layer1/layer2 interface roughness
    sigma[2] = layer2/layer3 interface roughness

Note that the Debye-Waller model fails for high values of
sigma (e.g. > 20 angstroms).  For rougher interfaces, or for
rapidly varying compsition distributions, the best way to model these
is by generating a discretized model with many homogeneous layers
that approximates the continuous distribution (e.g. density or composition
profiles that follow an erf distribution)

References:
-----------
1) Bartels et. al., Acta Cryst (1986) A42 539-545.
2) de Boer, PRB (1991) V44 498-511
3) Krol et al, PRB (1988) V38 8579-8592
4) Vidal and Vincent, Applied Optics (1984) V23 1794-1801
5) Trainor, Templeton and Eng, J Electron Spec Rel Phenom (2006) V150, 66-85

For information about scattering factors see (e.g.)
6) Chantler, et al. (2005), X-Ray Form Factor,
7) Attenuation and Scattering Tables (version 2.1):
   http://physics.nist.gov/ffast
8) Chantler, C.T., J. Phys. Chem. Ref. Data V29, 597-1048 (2000)
9) Chantler, C.T., J. Phys. Chem. Ref. Data V24, 71-643 (1995).

"""
###############################################################################

import numpy as num
import exceptions

from Ifeffit  import Ifeffit
from .wrxrr import _XrayRefl
from tdl.modules.utils import elements

###############################################################################

DTYPE = num.dtype('double')

DEFAULT_PARAMS = {'energy':10000.,'wconv':0.01,
                  'slen':50.,'bvert':0.05,'bhorz':10.0,
                  'aflag':0.,'fyidx':-1.,'fyenergy':7000.,
                  'adet':90.,'tnorm':1.0,'rflag':1.0,
                  'delz':5.0,'pdepth':3.0,'rscale':1.0}

###############################################################################
class RefModel(_XrayRefl):
    """
    Reflectivity model

    See wrxrr._XrayRefl for details
    regarding the initialization - we need to make
    sure that the arrays/data are appropriate
    to pass to the c library.
    """
    def __init__(self,d=[],rho=[],sigma=[],comp=[],elem_z=[],theta=[],params={}):
        """
        Parameters:
        -----------
        All these should be num arrays, dtype = double
        * d array of layer thicknesses (angstroms)
        * rho array of corresponding densities (g/cm^3)
        * sigma array of interface roughnesses (angstroms)
        * comp is the composition array
        * elem_z are the Z values of each element in the model
        * theta array of theta values (degrees)
        * params is a dictionary of model parameters (see set_params)
        """
        self.iff = Ifeffit(screen_echo = 0)
        self.nlayer    = 0
        self.nelem     = 0
        self.nthet     = 0
        self.fp        = []
        self.fpp       = []
        self.amu       = []
        self.mu_at     = []

        self._init_params()
        # some flags
        self._init_en   = True
        self._init_fy   = True
        self._init_carr = True
        self._init_ptr  = True
        # init
        self.init_model(d=d,rho=rho,sigma=sigma,comp=comp,elem_z=elem_z,theta=theta)
        self.set_params(**params)

    ###########################################################################
    def init_model(self, d=None, rho=None, sigma=None, comp=None,
                   elem_z=None, theta=None):
        """
        Initialize/modify model data and arrays.

        Note that arrays should be passed in as numpy arrays, and will
        be assigned by reference!! See wrxrr._XrayRefl._init_model for
        array creation
        """
        if d != None:
            self.nlayer = len(d)
            if (d.dtype == DTYPE): self.d = d
            else: raise exceptions.ValueError
            self._init_carr = True
            self._init_ptr = True
        #
        if rho != None:
            if len(rho) != self.nlayer: raise exceptions.ValueError
            if (rho.dtype == DTYPE): self.rho = rho
            else: raise exceptions.ValueError
            self._init_ptr = True
        #
        if sigma != None:
            if len(sigma) != self.nlayer-1: raise exceptions.ValueError
            if (sigma.dtype == DTYPE): self.sigma = sigma
            else: raise exceptions.ValueError
            self._init_ptr = True
        #
        if comp != None:
            if (comp.dtype != DTYPE): raise exceptions.ValueError
            nel,nz = comp.shape
            if nz != self.nlayer: raise exceptions.ValueError
            self.comp = comp
            if nel != self.nelem:
                self.nelem  = nel
                self._init_en   = True
                self._init_fy   = True
            self._init_ptr  = True
        #
        if elem_z != None:
            if len(elem_z) != self.nelem: raise exceptions.ValueError
            if (elem_z.dtype == DTYPE): self.elem_z = elem_z
            else: raise exceptions.ValueError
            self._init_en   = True
            self._init_fy   = True
            self._init_ptr  = True
        #
        if theta != None:
            self.nthet  = len(theta)
            if (theta.dtype == DTYPE): self.theta = theta
            else: raise exceptions.ValueError
            self._init_carr = True
            self._init_ptr  = True

    ##########################################################
    def set_params(self,**params):
        """
        Parameters:
        -----------
        * energy (eV) is the incident energy        # calc_params[0]
        * wconv (degrees) is the concolution width  # calc_params[1]
        * slen  (mm) is the sample length           # calc_params[2]
        * bvert (mm) is the vertical beam dim       # calc_params[3]
        * bhorz (mm) is the horizontal beam dim     # calc_params[4]
        * aflag is a flag for area calcs            # calc_params[6]
        * fyidx is the index of the element for FY  # calc_params[6]
        * fyenergy (eV) is the energy of FY         # calc_params[7]
        * adet  (deg) is the angle of the detector  # calc_params[8]
        * tnorm (deg) is the theta for normalizatio # calc_params[9]
        * rflag is the roughness flag               # calc_params[10]
        * delz (ang) is the delta z for integration # calc_params[11]
        * pdepth is the substrate FY clac flag      # calc_params[12]
        * rscale is a scale factor for reflectivity # calc_params[13]
        """
        if len(params) == 0: return
        if params.has_key('energy'):
            self.calc_params[0] = params['energy']
            self._init_en = True
        if params.has_key('wconv'):  self.calc_params[1] = params['wconv']
        if params.has_key('slen'):   self.calc_params[2] = params['slen']
        if params.has_key('bvert'):  self.calc_params[3] = params['bvert']
        if params.has_key('bhorz'):  self.calc_params[4] = params['bhorz']
        if params.has_key('aflag'):  self.calc_params[5] = params['aflag']
        if params.has_key('fyidx'):  self.calc_params[6] = params['fyidx']
        if params.has_key('fyenergy'):
            self.calc_params[7] = params['fyenergy']
            self._init_fy = True
        if params.has_key('adet'):   self.calc_params[8] = params['adet']
        if params.has_key('tnorm'):  self.calc_params[9] = params['tnorm']
        if params.has_key('rflag'):  self.calc_params[10] = params['rflag']
        if params.has_key('delz'):   self.calc_params[11] = params['delz']
        if params.has_key('pdepth'): self.calc_params[12] = params['pdepth']
        if params.has_key('rscale'): self.calc_params[13] = params['rscale']

    ##########################################################
    def get_params(self,**params):
        """
        return params
        """
        params = {}
        params['energy']   = self.calc_params[0]
        params['wconv']    = self.calc_params[1]
        params['slen']     = self.calc_params[2]
        params['bvert']    = self.calc_params[3]
        params['bhorz']    = self.calc_params[4]
        params['aflag']    = self.calc_params[5]
        params['fyidx']    = self.calc_params[6]
        params['fyenergy'] = self.calc_params[7]
        params['adet']     = self.calc_params[8]
        params['tnorm']    = self.calc_params[9]
        params['rflag']    = self.calc_params[10]
        params['delz']     = self.calc_params[11]
        params['pdepth']   = self.calc_params[12]
        params['rscale']   = self.calc_params[13]
        return params

    ##########################################################
    def init_energy(self,):
        """
        Change the incident energy of the calc

        get fp,fpp from ifeffit for index of refraction calcs
        """
        if len(self.fp) != self.nelem:
            self.fp    = num.zeros(self.nelem, dtype=num.double)
            self.fpp   = num.zeros(self.nelem, dtype=num.double)
            self.amu   = num.zeros(self.nelem, dtype=num.double)
            self._init_ptr  = True

        energy = self.calc_params[0]
        for j in range(self.nelem):
            en = [energy, energy + 1.]
            self.iff.put_array('calc.en',en)
            cmd = 'f1f2(energy=calc.en,z=%i)' % self.elem_z[j]
            self.iff.ifeffit(cmd)
            fp  = self.iff.get_array('calc.f1')
            fpp = self.iff.get_array('calc.f2')
            #
            self.fp[j]  = fp[0]
            self.fpp[j] = fpp[0]
            self.amu[j] = elements.amu(self.elem_z[j])

        self._init_en = False

    ##########################################################
    def init_fy(self,):
        """
        Get fpp from ifeffit and use this to est mu_atomic for fy.
        This will be a decent estimate for high z and E
        however, these are approx since we are ignoring
        coherent and incoherent scattering cross sections
        to computing total (ie this is photoeffect only)
        """
        if len(self.mu_at) != self.nelem:
            self.mu_at = num.zeros(self.nelem, dtype=num.double)
            self._init_ptr  = True

        #below is 2*Na*r_e*hc (in cm^2*eV*mole/atom)
        con = 4.20792637233e07
        fy_energy = self.calc_params[7]
        if fy_energy == 0.0:
            self.mu_at = num.zeros(self.nelem, dtype=num.double)
        else:
            for j in range(self.nelem):
                en = [fy_energy, fy_energy + 1.]
                self.iff.put_array('calc.en',en)
                cmd = 'f1f2(energy=calc.en,z=%i)' % self.elem_z[j]
                self.iff.ifeffit(cmd)
                fp  = self.iff.get_array('calc.f1')
                fpp = self.iff.get_array('calc.f2')
                fpp = float(fpp[0])
                mu  = con*fpp/(fy_energy*self.amu[j])
                #
                self.mu_at[j] = mu
        self._init_fy = False

    ##########################################################
    def calc_R(self,ret=False):
        """
        Calc reflectivity
        """
        # save fy_idx to reset
        tmp = self.calc_params[6]
        self.calc_params[6] = -1.0
        #
        if self._init_en: self.init_energy()
        if self._init_fy: self.init_fy()
        #
        self._calc(init_ptrs=self._init_ptr,init_arrs=self._init_carr)
        #
        if self._init_ptr == True: self._init_ptr = False
        if self._init_carr == True: self._init_carr = False
        #
        self.calc_params[6] = tmp
        #
        if ret: return (self.R.copy())

    ##########################################################
    def calc_FY(self,ret=False):
        """
        Calc fluorescent yield
        """
        # check some stuff
        fy_idx = self.calc_params[6]
        if (fy_idx < 0.) or (fy_idx > self.nelem-1):
            print "Error, fy_idx out of range"
            return
        delz = self.calc_params[11]
        if delz < 1. :
            print "Error, zint too small, min = 1 ang."
            return
        #
        if self._init_en: self.init_energy()
        if self._init_fy: self.init_fy()
        #
        self._calc(init_ptrs=self._init_ptr, init_arrs=self._init_carr)


        #
        if self._init_ptr == True: self._init_ptr = False
        if self._init_carr == True: self._init_carr = False
        #
        if ret: return (self.Y.copy(), self.R.copy())

    ##########################################################
    def make_mole_fractions(self):
        """
        make composition into mole fractions
        """
        self.comp = self.comp / self.comp.sum(0)
        return


############################################################################
############################################################################
def test():
    import time
    from matplotlib import pyplot
    pyplot.ion()

    # make a simple multilayer and calc R and FY

    #### layer thicknesses
    # note the first layer is the base layer,
    # it is always assumed infiniteley thick
    # for the R calc (ie all phase shift calcs).
    # For the FY calc, it sets the depth to which the yield
    # is calculated, unless calc_params[12], pdepth, is > 0.
    # If pdepth is greater than zero than its assumed to
    # be a multiplier - the depth of the FY integration in the
    # base layer is calc from this multiplier times the
    # penetration depth (which is angle dependant).  This is
    # much faster since it limits the integration range to small
    # values at low angles.
    # In either case you may have to play with these paramter
    # to see when it no longer contributes to the FY
    # (if the element of interest is in the base layer)
    # The top layer should be a low density layer
    # that approximates air
    # The thicknesses are in angstroms
    d = num.array([50000.0, 100.0, 10.0, 1000.0],dtype=num.double)

    # the corresponding densities of the layers
    # these are g/cm^3
    rho = num.array([3.0, 4.0, 1.7, 0.0002],dtype=num.double)

    # Debye-Waller type interface roughness parameters: r = r*exp(-(q*sigma)^2)
    # The first value is the 1/0 interface roughness
    # ie the interface between the first layer and the base layer
    # there are nlayer - 1 total interfaces.
    sigma = num.array([10, 10, 10],dtype=num.double)

    #this is the list of all elements that make up the multilayer composition
    #         N, O, Al, Si, Fe
    elem_z = num.array([7, 8, 13, 14, 26],dtype=num.double)

    # Composition is a matrix giving the mole fraction or stiochiometric
    # coefficients of each of the above elements in each layer.
    # Note if the sum of the comp numbers for a layer does not sum to 1
    # then they will be normalized to mole fraction
    # -> base layer = SiO2
    # -> layer 1    = Al (with trace of Fe)
    # -> layer 2    = SiAlFeO (arb ratios)
    # -> layer 3    = N2
    # therefore, comp[0] = fraction of N in each layer etc..
    comp = [ [0.0,   0.0,   0.0,  1.0],  # N
             [2.0,   0.0,   0.52, 0.0],  # O
             [0.0,   0.999, 0.07, 0.0],  # Al
             [1.0,   0.0,   0.4,  0.0],  # Si
             [0.00001,   0.001, 0.1,  0.0]]  # Fe
    comp = num.array(comp,dtype = num.double)

    # theta is the arrays of angles for the calculation
    theta = num.arange(0.01, 1.0, 0.01, dtype=num.double)

    # calc_params holds a bunch of parameter stuff
    # calc_params[0] = energy (eV)
    # calc_params[1] = wconv (degrees)
    # calc_params[2] = slen  (mm)
    # calc_params[3] = bvert (mm)
    # calc_params[4] = bhorz (mm)
    # calc_params[6] = aflag
    # calc_params[6] = fyidx
    # calc_params[7] = fyenergy (eV)
    # calc_params[8] = adet  (deg)
    # calc_params[9] = tnorm (deg)
    # calc_params[10] = rflag
    # calc_params[11] = delz (ang)
    # calc_params[12] = pdepth
    # calc_params[13] = rscale
    calc_params = {'energy':10000.,'wconv':0.01,'slen':10.,'bvert':0.01,
                   'aflag':1.,'fyidx':4,'fyenergy':7000.,
                   'delz':5.0,'pdepth':3.0}

    # create an instance of the model
    ref = RefModel(d=d,rho=rho,sigma=sigma,comp=comp,elem_z=elem_z,
                   theta=theta,params=calc_params)

    # plot R
    pyplot.subplot(3,1,1)
    pyplot.semilogy()
    ref.calc_R()
    pyplot.plot(ref.theta,ref.R,'g.')

    # run the FY calcs
    t = time.time()
    ref.calc_FY()
    print " Elapsed time = %f seconds\n" % (time.time() - t)

    #FY
    pyplot.subplot(3,1,1)
    pyplot.plot(ref.theta,ref.R,'r')
    pyplot.subplot(3,1,2)
    pyplot.plot(ref.theta,ref.Y)

    # hist of Fe dist
    #print comp
    pyplot.subplot(3,1,3)
    idx = num.arange(len(comp[4]))
    pyplot.bar(idx,comp[4])

    return ref


###############################################################################
###############################################################################
if __name__ == "__main__":
    ref = test()
    raw_input("hit enter to continue")

