#######################################################################
"""
T. Trainor (fftpt@uaf.edu)
Calculate reflectivity / reflection FY as
a function of interface composition

Modifications:
--------------


"""
#######################################################################
"""
Todo
-----

- Check carefully:
  * Pass by ref between slab and ref??  ie arrays getting updated!
  * follow all update flags.... make sure we're efficient!

- Test distro functions
  * setting dist params...
  * test idx ranges (speed up)...

- Add doc!!!

============

- Add plotting specific component/element profiles.  

----------------------------------------------------
- Fitting/modeling (single and mult spectra)
  * new class holds:
        [{'params':{},'R':[],'FY':[]}, ...]
        params = {'energy':10000,'fyidx':'Fe','fyenergy':'Fe Ka',etc}

  * First make a script/function to model/plot the sequence

  * Area correction for data

  * loop through and compute R/FY for each
    - fitting loop compute total (wieghted) residual
    
  * need a parameter map - set of (linked) parameters...
    - list of params with some corresponding list of identifier for where
      they plug into the model (maybe mult. places).  eg. sig of interface
      dist for 2 components may be the same 'parameter'
    - may be able to do using lambda's / eval expressions (inefficient)....
       param[0] = 10.
       param_map[0] = ['self.model.slab.distpar[0].inter[0]['sig']',
                       'self.model.slab.distpar[1].inter[0]['sig']']
       for var in param_map[0]:
           src = "%s=%f" % (var, param[0])
           eval(src)

--------------------------------------------------------------

- Allow constructing stiochiometry of species from wt% numbers
  ==> see species.set in compound.py
  
- Add an I/z/theta matrix generation and plotting routine...

- Create a text file dump (and read) of ref object?

"""
#######################################################################

import numpy as num
import scipy
import types, copy

from reflectivity import RefModel, DEFAULT_PARAMS
from xtab     import xrf_lookup
from mpcutils import element_data as elements
from mpcutils import compound

#######################################################################
class Layer(compound.Material):
    """
    Same as material/species from compound module,
    except we include area, thickness, and roughness
    """
    #######################################################################
    def __init__(self,**kw):
        self.init_layer(**kw)
        
    #######################################################################
    def init_layer(self,comp=[],density=1.0,roughness=0.0,
                   thickness=100.,area=1.0):
        self.roughness  = roughness   # rms, angstroms
        self.thickness  = thickness   # angstroms
        self.area       = area        # cm^2
        self._compute_vol()           # cm^3
        #
        self.init_mat(comp=comp,density=density,volume=self.vol)

    #######################################################################
    def _compute_vol(self):
        z = self.thickness  # angstroms
        a = self.area       # cm^2
        v = z*(1.0e-8)*a
        self.vol = v

###########################################################################
class DistParams:
    """
    Class for component distribution parameters
       dist = {'type':'box','zst':10.0,'zen':100,'CX':1.0}
    etc..
    """
    ########################################################################
    def __init__(self, name=''):
        """
        initialize
        """
        # string name of component
        self.cname = name
        
        # distribution
        self.top   = {'type':'box','CX':0.0}  # top
        self.inter = []                       # interface distribution
        self.subs  = {'type':'box','CX':0.0}  # substrate

        # flags
        self.scale_to_norm = True  # adjust scale factors to norm
        self.norm  = 0             # dist normalization flag
                                   # 0 = None
                                   # 1 = fix scale factors for interface mass balance
                                   # 2 = fix CX[1]    = CX_substrate
                                   # 3 = fix CX[nz-2] = CX_top

        # params
        self.totNX = 0.0    # total moles of component in interface
        self.aveCX = 0.0    # average interface concentration (moles/cm^3)
                            #  = self.totNX[j]/totV
                            # Interface is defined as slabs between top and bottom 

    #######################################################################
    def __repr__(self,):
        lout = '\n** Comp = %s, (norm= %s, scale_to_norm =%s)\n' % (self.cname,
                                                                   repr(self.norm),
                                                                   repr(self.scale_to_norm))
        lout = "%s - top      type=%6s, zst=%6s, zen=%6s, CX=%s\n" %  (lout,
                                                              str(self.top.get('type')),
                                                              str(self.top.get('zst')),
                                                              str(self.top.get('zen')),
                                                              str(self.top.get('CX')))
        idx = range(len(self.inter))
        idx.reverse()
        for j in idx:
            lout = "%s - inter %i  type=%6s, zst=%6s, zen=%6s, CX=%s\n" %  (lout, j,
                                                                   str(self.inter[j].get('type')),
                                                                   str(self.inter[j].get('zst')),
                                                                   str(self.inter[j].get('zen')),
                                                                   str(self.inter[j].get('CX')))
            
        lout = "%s - subs     type=%6s, zst=%6s, zen=%6s, CX=%s\n" %  (lout,
                                                               str(self.subs.get('type')),
                                                               str(self.subs.get('zst')),
                                                               str(self.subs.get('zen')),
                                                               str(self.subs.get('CX')))
        return lout

    ########################################################################
    def init_dist(self,subs=None,top=None,inter=None):
        """
        (re)init distribution parameters
        """
        #
        if top!= None: self.top = top
        if subs!=None: self.subs = subs
        #
        self.inter = []
        if inter != None:
            if type(inter) != types.ListType:
                inter = [inter]
            for d in inter:
                self.inter.append(d)
        self._sort()
        
    ########################################################################
    def add_dist(self,**kw):
        """
        Add a new distribution to interface
        """
        dd = {'type':'box','CX':0.0}
        for key in kw.keys():
            dd[key] = kw[key]
        self.inter.append(dd)
        self._sort()
    
    ########################################################################
    def _sort(self,d0=False):
        """
        loop through and sort by zst
        work_here --> finish ??
        """
        pass
    
###########################################################################
class Slab:
    """
    Create slab model
    * Conventions / indexing chemistry
     - The list self.comp holds component objects.  The list self.distpar
       holds the corresponding distibution parameters. ie distpar[j] are the
       dist parameters for comp[j]

    * Arrray references
     - Note we try to handle the arrays:
           self.d, self.rho, self.sigma, self.fZ, self.elem_z
       such that thier reference (to thier memory location)
       does not change after calling init().  i.e. these arrays
       are used by the Reflectivity module in C, and if we're careful
       we can avoid having to reassign these all the time...
       
    * Conventions / indexing slabs
     - Note that we do not slabify layer[0] and layer[nlayer-1]
       (ie not divied up and not included in mass balance).
     - However, they are included as the first and last entry of z, d etc.
       and thier components are part of total component list
     - Conventions for z:
       * z values are for the bottom of each slab
       * all d's are positive (layer thicknesses)
       * Therefore:  self.z + self.d = top of each slab
     - Note self.delta is the target slab size.  However, 
       actuals d's are adjusted to ensure that layers maintain
       the correct total thickness
     - If delta <= 0, slabs are same thickness as layers
     - If delta > 0 , we turn off roughness at top of
       any interface layer.
     - To index slabs according to original layer use:
       self.z[num.where(self.zidx==layer_idx)]
       
    """
    ########################################################################
    def __init__(self,layer=[],delta=10.):

        ###
        self.layer  = layer   # list of layer objects
        self.delta  = delta   # delta for slabify (if <0 slabs = layers)
        
        ### components
        self.comp    = []     # list of components (component objects), len=numX
        self.elem    = []     # list of elements (symbols), len=numel
        self.elem_z  = []     # list of elements Z, len=numel
        self.totV    = 1.0    # total volume of interface, cm^3
        ### z arrays 
        self.z      = []      # z-values, len=numz
        self.zidx   = []      # which layer z belongs to, len=numz
        self.d      = []      # thickness of each z segment, len=numz
        self.sig    = []      # roughnesses, len=numz-1
        self.layer_range = [] # zrange etc spanned by each layer
        ### distribution
        self.distpar = []     # component dist parameters, len = numX 
        self.CX      = []     # component concetrations, mole/cm^3, shape=(numX,numz)
        self.rho     = []     # densities, g/cm^3, shape=(numz)
        self.rhoflag = 0      # flag for density constraints
        self.rhoscale = True  # flag for rescaling dist ampls given density constraint
        self.CZ      = []     # element concentrations, mole/cm^3, shape=(numEl,numz)
        self.fZ      = []     # element mole fractions, shape=(numEl,numz) 
        ###
        self.init()

    ########################################################################
    def __repr__(self,):
        lout = "**** Slab Model ****\n"
        lout = "%s - components = %s\n" % (lout, repr(self.list_comps()) )
        lout = "%s - elements   = %s\n" % (lout, repr(self.list_elements()) )
        lout = "%s - slab delta = %6.3f, numz = %i\n" % (lout, self.delta, len(self.z))
        for dpar in self.distpar:
            lout = lout + dpar.__repr__()
        return lout

    ########################################################################
    def init(self):
        """
        Create (new) model
        """
        if len(self.layer) < 2:
            print "Two or more layers are required"
            return
        self._init_comp()
        self._init_z()
        self._init_dist()
        
    ########################################################################
    def _init_comp(self):
        """
        Get all the:
        - unique components
        - list of elements
        - total 'interface' volume
        - distparam objects and component totals
          (component totals are just for the interface region!!)
        """
        self.comp    = []   # list of components (component objects), len = numX
        self.distpar = []   # component dist parameters, len = numX 
        self.elem    = []   # list of elements (symbols)
        self.elem_z  = []   # list of elements (Z)
        self.totV    = 0.0  # total volume of interface, cm^3
        ###
        cnames = []
        def _add(mat):
            for (comp,nu) in mat.comp:
                if comp.name not in cnames:
                    self.comp.append(comp)
                    cnames.append(comp.name)
                tmp = comp.elements()
                for elem in tmp:
                    if elem not in self.elem:
                        self.elem.append(elem)
                        self.elem_z.append(elements.number(elem))

        ### add each material
        for mat in self.layer:
            _add(mat)
        self.elem_z = num.array(self.elem_z,dtype=num.double)

        ### Get 'interface' volume
        ### (excluding top and bottom layers)
        nlayer = len(self.layer)
        for j in range(nlayer):
            if (j > 0) and (j < nlayer -1):
                self.totV = self.totV + self.layer[j].vol

        ### Create component dist params and
        ### get total 'interface' component moles
        ### (excluding top and bottom layers)
        ### Note order of self.dpar will
        ### correspond with self.comp
        for n in cnames:
            dpar  = DistParams(name=n)
            totNX = 0.0
            for j in range(nlayer):
                if (j > 0) and (j < nlayer -1):
                    totNX = totNX + self.layer[j].molesX(n)
            dpar.totNX = totNX
            if self.totV > 0.:
                dpar.aveCX = dpar.totNX/self.totV
            self.distpar.append(copy.copy(dpar))

    ########################################################################
    def list_comps(self):
        """
        List of component names
        """
        names = []
        for comp in self.comp:
            names.append(comp.name)
        return names

    ########################################################################
    def list_elements(self):
        """
        List of element names
        """
        return self.elem

    ########################################################################
    def _comp_idx(self,comp):
        """
        Get component index.
        - This should work for both self.comp and self.distpar lists
        - comp maybe a string (comp name) or integer (comp index)
        """
        cidx = -1
        if type(comp) == types.StringType:
            for j in range(len(self.comp)):
                if comp == self.comp[j].name:
                    if self.distpar[j].cname == comp :
                        cidx = j
                        break
                    else:
                        print "Component Name Error!"
                        return
        else:
            cidx = int(comp)
            if cidx not in range(len(self.comp)): cidx = -1
        if cidx == -1:
            print "Component %s not found" % str(comp)
            return -1
        else:
            return cidx
        
    ########################################################################
    def _init_z(self,):
        """
        - Conventions:
          * z values are for the bottom of each slab
          * all d's are positive (layer thicknesses)
        - Therefore:
          self.z + self.d = top of each slab
        - Note self.delta is the target slab size.  however, 
          this is adjusted to ensure that layers maintain
          the correct total thickness
        - If delta <= 0, slabs are same thickness as layers
        - If delta > 0 , we turn off roughness at top of
          any interface layer.
        - To index use: self.z[num.where(self.zidx==layer)]
        """
        self.z      = []    # z-values, len=numz
        self.zidx   = []    # which layer z belongs to, len=numz
        self.d      = []    # thickness of each z segment, len=numz
        self.sig    = []    # roughnesses, len=numz-1 (ignore top)
        #
        delta0 = float(self.delta)
        nlayer = len(self.layer)
        zb     = 0.0
        zt     = 0.0
        #
        for j in range(nlayer):
            mat    = self.layer[j]
            thick  = num.abs(float(mat.thickness))
            rough  = num.abs(float(mat.roughness))
            # sustrate
            if j == 0:
                zz     = num.array([-1.0*thick])
                zzidx  = num.array([0])
                dd     = num.array([thick])
                sig    = num.array([rough])
            # top
            elif j == nlayer - 1:
                zz     = num.concatenate((zz, num.array([zt])) )
                zzidx  = num.concatenate((zzidx, num.array([j])) )
                dd     = num.concatenate((dd, num.array([thick])))
                # ignore top sig
                # sig  = num.concatenate((sig, num.array([rough])))
            # slab using layer thick
            elif delta0 <= 0:
                zz     = num.concatenate((zz, num.array([zt])) )
                zt     = zt + thick
                zzidx  = num.concatenate((zzidx, num.array([j])) )
                dd     = num.concatenate((dd, num.array([thick])))
                sig    = num.concatenate((sig, num.array([rough])))
            # slabify interface layer
            else:
                zt     = zt + thick
                nn     = num.ceil((zt-zb)/delta0)
                delta  = (zt-zb)/nn
                #
                tmp    = num.arange(zb,zt,delta,dtype='float')
                zz     = num.concatenate((zz,tmp))
                #
                tmp    = num.ones(len(tmp),dtype='float')
                zzidx  = num.concatenate((zzidx,tmp*(j)))
                dd     = num.concatenate((dd,tmp*delta))
                #
                tmp    = num.zeros(len(tmp),dtype='float')
                sig    = num.concatenate((sig,tmp))
                # no sigs for slabified!!
                #sig[len(sig)-1] = rough 
                #
                zb = zt
                
        #############
        self.z     = zz
        self.zidx  = zzidx.astype('int')
        self.d     = dd
        self.sig   = sig
        self._get_zrange()

    ########################################################################
    def _get_zrange(self,):
        """
        Get (zmin,zmax),(idxmin,idxmax), nz for layers
        """
        self.layer_range = []
        nlayer = len(self.layer)
        for idx in range(nlayer):
            xx = {}
            idx = num.where(self.zidx==idx)
            zz  = self.z[idx]
            nz  = len(zz)
            xx['nz']     = nz
            xx['zmin']   = zz[0]
            xx['zmax']   = zz[nz-1]
            xx['idxmin'] = idx[0][0]
            xx['idxmax'] = idx[0][-1]
            self.layer_range.append(xx)

    def get_zidx(self,z):
        """
        Get the (first) idx which is closest to z
        """
        if z == None: return -1
        zz  = num.abs(self.z - z)
        idx = num.where(zz == min(zz))
        return idx[0][0]
    
    ########################################################################
    def _init_dist(self):
        """
        Set up the initial distribution model
        """
        ### distribution arrays
        self.rho    = []   # densities, g/cm^3, shape=(numz,)
        self.CX     = []   # component concetrations mole/cm^3, shape=(numX,numz)
        self.CZ     = []   # elements concetrations mole/cm^3, shape=(numEl,numz)
        self.fZ     = []   # elements mole fractions, shape=(numEl,numz) 
        #
        nlayer = len(self.layer)
        numz   = len(self.z)
        numX   = len(self.comp)
        numEl  = len(self.elem)
        #
        self.rho = num.zeros(numz,dtype='float')
        self.CX  = num.zeros((numX,numz),dtype='float')
        self.CZ  = num.zeros((numEl,numz),dtype='float')
        self.fZ  = num.zeros((numEl,numz),dtype='float')

        ### Loop through all components
        ### and all layers and generate initial
        ### box model distro for each 
        for j in range(numX):
            cname = self.comp[j].name
            dpar  = self.distpar[j]
            inter = []
            subs  = {}
            top   = {}
            if cname != dpar.cname:
                print "Error in component/distribution indexing!"
                return
            for k in range(nlayer):
                dd        = {'type':'box'}
                dd['zst'] = self.layer_range[k]['zmin']
                dd['zen'] = self.layer_range[k]['zmax']
                dd['CX']  = self.layer[k].concX(cname)
                if k == 0:
                    subs = copy.copy(dd)
                elif k == nlayer-1:
                    top = copy.copy(dd)
                else:
                    inter.append(copy.copy(dd))
            dpar.init_dist(subs=subs,top=top,inter=inter)
        #
        self.calc_dist()

    ########################################################################
    def add_dpar(self,comp,dist={},norm=None,scale=None,init=False):
        """
        Add a distribution parameter for the component comp
        - comp maybe a string (comp name) or integer (comp index)
        - If init==True, the interface dist is cleared for this comp
          before adding the new dist (top and subs unaffected)
        """
        cidx = self._comp_idx(comp)
        if cidx == -1: return
        
        ### flags
        if norm != None: self.distpar[cidx].norm=norm
        if scale != None: self.distpar[cidx].scale_to_norm=scale

        ### add the distribution
        if init: self.distpar[cidx].init_dist()
        if len(dist) > 0:
            self.distpar[cidx].add_dist(**dist)
    
    ########################################################################
    def set_dpar(self,comp,dist=None,didx=0,norm=None,scale=None,):
        """
        Edit a distribution parameter for the component comp
        - comp maybe a string (comp name) or integer (comp index)
        - didx is the index of the interface distro
        """
        cidx = self._comp_idx(comp)
        if cidx == -1: return
        
        ### flags
        if norm != None: self.distpar[cidx].norm=norm
        if scale != None: self.distpar[cidx].scale_to_norm=scale
        
        ### set params
        if dist!= None:
            self.distpar[cidx].inter[didx].update(dist)

    ########################################################################
    def calc_dist(self,):
        """
        Compute the concentration distribution of all components
        Note density constraints are applied after computing concetrations
          if self.rhoflag == 0 no density normalization
          if self.rhoflag == 1 normalize to substrate density
          if self.rhoflag == 2 normalize to top density
          if self.rhoscale == True the normalization is applied to all
                              component distribution amplitudes 
        Note we are trying to keep the memory locations of 
        self.fZ and self.rho fixed in the calc.  
        """
        
        # reset self.CX, assume array sizes havent changed
        # compute self.CX
        self.CX.fill(0.0)
        ncomp = len(self.distpar)
        for cidx in range(ncomp):
            self._calc_dist(cidx)

        # other arrays (these dont need to be reset)
        #self.CZ.fill(0.0)
        #self.fZ.fill(0.0)
        #self.rho.fill(0.0)
        CZ   = num.zeros(self.CZ.shape)
        fZ   = num.zeros(self.fZ.shape)
        rho  = num.zeros(self.rho.shape)

        # compute elem conc and slab densities
        nelem = len(self.elem)
        for j in range(nelem):
            el = self.elem[j]
            amu = elements.amu(el)
            for k in range(ncomp):
                nu    = self.comp[k].nuZ(el)
                CZ[j] = CZ[j] + self.CX[k] * nu
            rho = rho + CZ[j] * amu
            fZ[j] = CZ[j] * self.d
        denom = fZ.sum(0)
        if num.min(denom) > 0:
            fZ = fZ / denom

        # see if there are density constraints

        # scale to substrate
        if (self.rhoflag == 1) and (rho[1]>0):
            # rescale density
            f = rho[0]/rho[1]
            rho[1:-1] = f*(rho[1:-1])
            # rescale CZ
            for j in range(nelem):
                CZ[j][1:-1] = f*(CZ[j][1:-1])
            # rescale CX
            for k in range(ncomp):
                self.CX[k][1:-1] = f*(self.CX[k][1:-1])
                # adjust ampls
                if self.rhoscale == True:
                    self._scale_dist_ampl(k,scale=f)
        # scale to top
        elif (self.rhoflag == 2) and (rho[-2]>0):
            # rescale density
            f = rho[-1]/rho[-2]
            rho[1:-1] = f*(rho[1:-1])
            # rescale CZ
            for j in range(nelem):
                CZ[j][1:-1] = f*(CZ[j][1:-1])
            # rescale CX
            for k in range(ncomp):
                self.CX[k][1:-1] = f*(self.CX[k][1:-1])
                # adjust ampls
                if self.rhoscale == True:
                    self._scale_dist_ampl(k,scale=f)
 
        # This should keep the original
        # array references valid
        self.CZ.flat[:] = CZ.ravel()[:]
        self.fZ.flat[:] = fZ.ravel()[:]
        self.rho.flat[:] = rho.ravel()[:]
    
    ########################################################################
    def _calc_dist(self,cidx):
        """
        Note all (interface) distributions are assumed to be of form:
           CX[cidx][zrange]= CX[cidx][zrange] + dist['CX']*f(z[zrange])
          
        ie the total component distribution is a sum over indidual dists
        and each distribution is assumed to have a 'scale' = dist['CX']
        that weights to contribution.  This allows us apply normalization
        factors to the distr parameters (ie make normalization stick)

        if self.distpar[cidx].norm == 0  no normalization
        if self.distpar[cidx].norm == 1  mass balance is conserved
        if self.distpar[cidx].norm == 2  normalize to substrate concentration
        if self.distpar[cidx].norm == 3  normalize to top concentration

        """
        # gather all params
        dpar = [self.distpar[cidx].subs]
        for dist in self.distpar[cidx].inter:
            dpar.append(dist)
        dpar.append(self.distpar[cidx].top)

        # loop over dist's and add to self.CX
        # for dist in dpar:
        ndist = len(dpar)
        dist = None
        for j in range(ndist):
            # set interface flag
            if (j == 0):
                interface='s'
            elif (j==ndist-1):
                interface='t'
            else:
                interface='i'
            #
            dist = dpar[j]
            if dist['type'] == 'box':
                self._calc_box(dist,cidx,interface=interface)
            elif dist['type'] in ('erf','erfc'):
                self._calc_erf(dist,cidx,interface=interface)
            elif dist['type'] in ('exp','expc'):
                self._calc_exp(dist,cidx,interface=interface)
            elif dist['type'] == 'gauss':
                self._calc_gauss(dist,cidx,interface=interface)
            elif dist['type'] == 'linear':
                self._calc_linear(dist,cidx,interface=interface)

        # check normalization

        # fixed mass bal
        if self.distpar[cidx].norm == 1:
            f = 1.0
            moles     = self._comp_mb(cidx)
            moles_tot = self.distpar[cidx].totNX
            if moles > 0:
                f = (moles_tot/moles)
            self.CX[cidx][1:-1] = f*(self.CX[cidx][1:-1])
            # adjust amplitudes
            if self.distpar[cidx].scale_to_norm == True:
                self._scale_dist_ampl(cidx,scale=f)

        # fixed first point
        elif (self.distpar[cidx].norm == 2):
            f = 1.0
            if self.CX[cidx][1] != 0:
                f = self.CX[cidx][0]/self.CX[cidx][1]
            self.CX[cidx][1:-1] = f*(self.CX[cidx][1:-1])
            # adjust amplitudes
            if self.distpar[cidx].scale_to_norm == True:
                self._scale_dist_ampl(cidx,scale=f)

        # fixed end point
        elif (self.distpar[cidx].norm == 3):
            f = 1.0
            if self.CX[cidx][-2] != 0:
                f = self.CX[cidx][-1]/self.CX[cidx][-2]
            self.CX[cidx][1:-1] = f*(self.CX[cidx][1:-1])
            # adjust amplitudes
            if self.distpar[cidx].scale_to_norm == True:
                self._scale_dist_ampl(cidx,scale=f)

    ########################################################################
    def _scale_dist_ampl(self,cidx,scale=1.):
        """
        Rescale all dist amplitude factors for component
        distribution cidx by the given scale factor.
        Note this only adjusts interface distributions
        (subs and top unaffected)
        """
        for dist in self.distpar[cidx].inter:
            dist['CX'] = scale*(dist['CX'])

    ########################################################################
    def _dist_range(self,dist,interface='i'):
        """
        Get indicies for distro. These index z array as:
            zz[idxmin:idxmax+1]
        Note
            interface='i' means constrained to between subs/top
            interface='s' means subs
            interface='t' means top
        """
        numz = len(self.z)
        #
        if interface == 's':
            idxmin = 0
            idxmax = 0
            return (idxmin,idxmax)
        #
        if interface == 't':
            idxmin = numz-1
            idxmax = numz-1
            return (idxmin,idxmax)
        #
        if dist.has_key('zst'):
            idxmin = self.get_zidx(dist['zst'])
            if (idxmin < 1): idxmin = 1
            if (idxmin > numz-2): idxmin = 1
        else:
            idxmin = 1
        #
        if dist.has_key('zen'):
            idxmax = self.get_zidx(dist['zen'])
            if (idxmax < 1): idxmax = numz-2
            if (idxmax > numz-2): idxmax = numz-2
        else:
            idxmax = numz-2
        #
        return (idxmin,idxmax)

    ########################################################################
    def _calc_box(self,dist,cidx,interface='i'):
        """
        dist = {'type':'box','zst':10.0,'zen':100,'CX':1.0}
        """
        (idxmin,idxmax) = self._dist_range(dist,interface=interface)
        CX  = dist['CX']
        self.CX[cidx][idxmin:idxmax+1] = self.CX[cidx][idxmin:idxmax+1] + CX

    ########################################################################
    def _calc_linear(self,dist,cidx,interface='i'):
        """
        dist = {'type':'linear','zst':10.0,'zen':100,'CX':1.0,'CXen':1.0}
        """
        CX   = dist['CX']
        CXen = dist['CXen']
        (idxmin,idxmax) = self._dist_range(dist,interface=interface)
        zst   = self.z[idxmin]
        zen   = self.z[idxmax]
        denom = num.fabs(zen - zst)
        if denom == 0.0:
            print "Error, linear model requires z-range!"
            return
        slope = (CXen - CX)/denom

        # compute at center of each increment
        zz  = self.z[idxmin:idxmax+1] + self.d[idxmin:idxmax+1]/2.
        #
        y   = CX + slope*(zz - zst)
        if y.min() < 0:
            y[num.where(y<0)] = 0.0
        #
        self.CX[cidx][idxmin:idxmax+1] = self.CX[cidx][idxmin:idxmax+1] + y
    
    ########################################################################
    def _calc_erf(self,dist,cidx,interface='i'):
        """
        dist = {'type':'erf/c','zst':10.0,'zen':100,
                'cen':50.0,'sig':100.,'CX':1.0}
        """
        CX  = dist['CX']
        cen = dist['cen']
        sig = dist['sig']
        if sig == 0.0: sig = 1.0e-9
        if dist['type'] == 'erfc': compliment = True
        else: compliment = False

        # conc values are computed at the center of each increment
        (idxmin,idxmax) = self._dist_range(dist,interface=interface)
        zz  = self.z[idxmin:idxmax+1] + self.d[idxmin:idxmax+1]/2.
        #
        y   = 0.5*(scipy.special.erf((cen-zz)/(sig/2.))+1.)
        if compliment:
            y = 1. - y
        y  = CX * y
        #
        self.CX[cidx][idxmin:idxmax+1] = self.CX[cidx][idxmin:idxmax+1] + y

    ########################################################################
    def _calc_exp(self,dist,cidx,interface='i'):
        """
        dist = {'type':'exp/c','zst':10.0,'zen':100,
                'cen':50.0,'sig':100.,'CX':1.0}
        """
        CX  = dist['CX']
        cen = dist['cen']
        sig = dist['sig']
        if sig == 0.0: sig = 1.0e-9
        if dist['type'] == 'expc': compliment = True
        else: compliment = False

        # conc values are computed at the center of each increment
        (idxmin,idxmax) = self._dist_range(dist,interface=interface)
        zz  = self.z[idxmin:idxmax+1] + self.d[idxmin:idxmax+1]/2.
        #
        if compliment:
            y   = num.exp((zz-cen)/(sig))
        else:
            y   = num.exp((cen-zz)/(sig))
        y[y>1.] = 1.
        y  = CX * y
        #
        self.CX[cidx][idxmin:idxmax+1] = self.CX[cidx][idxmin:idxmax+1] + y

    ########################################################################
    def _calc_gauss(self,dist,cidx,interface='i'):
        """
        dist = {'type':'gauss','zst':10.0,'zen':100,
                'cen':50.0,'sig':100.,'CX':1.0}
        """
        CX   = dist['CX']
        cen  = dist['cen']
        sig  = dist['sig']
        
        # conc values are computed at the center of each increment
        (idxmin,idxmax) = self._dist_range(dist,interface=interface)
        zz  = self.z[idxmin:idxmax+1] + self.d[idxmin:idxmax+1]/2.
        
        # Note below amp is 'correct' gaussian
        # ampl, ie makes the distribution normalized
        # such that the integral is equal to CX
        #      amp = CX/(sig*num.sqrt(2*num.pi))
        # Use this amp instead so that max of dist = CX
        amp = CX
        #
        y   = amp*num.exp( -1.*(zz-cen)**2. / (2.*(sig**2.)) )
        #
        self.CX[cidx][idxmin:idxmax+1] = self.CX[cidx][idxmin:idxmax+1] + y

    ########################################################################
    def _comp_mb(self,cidx):
        """
        Compute component mass bal (moles)
        Ignore substrate and top
        Note conc are in mole/cm^3, d is in angstroms, 
        and we assume unit area = 1cm^2
        """
        xx = self.CX[cidx] * self.d
        mb = xx[1:-1].sum()*1.e-8
        return mb

    ########################################################################
    def comp_mb(self,show=True):
        """
        Compute component mass (mole) balance change
         relative to original model
        Ignore substrate and top
        """
        ncomp = len(self.comp)
        mb = []
        pmb = []
        for j in range(ncomp):
            m  = self._comp_mb(j)
            mt = self.distpar[j].totNX
            delta = (m - mt)
            if mt != 0:
                pdelta = 100.*delta/mt
            else:
                pdelta = 0.0
            mb.append(delta)
            pmb.append(pdelta)
        if show:
            print "Percent component mass balance change (relative to original model)"
            print "Substrate and top are not included in this calculation"
            for j in range(ncomp):
                print "   %s :  %6.6f" % (self.comp[j].name,pmb[j])
        return mb
    
    ########################################################################
    def plot(self,ty='density',bar=True,hold=False):
        """
        note hold only works with line plots
        """
        from matplotlib import pyplot
        zmin = num.min(self.z)
        zmax = num.max(self.z+self.d)
        if bar == True: pyplot.clf()
        if ty == 'density':
            if bar == True:
                pyplot.bar(self.z[0],    self.rho[0],    width= num.abs(self.d[0]))
                pyplot.bar(self.z[1:-1], self.rho[1:-1], width= num.abs(self.d[1:-1]))
                pyplot.bar(self.z[-1],   self.rho[-1],   width= num.abs(self.d[-1]))
            else:
                pyplot.plot([self.z[0],0.0], [self.rho[0],self.rho[0]],hold=hold)
                pyplot.plot([0.0,self.z[1]], [self.rho[0],self.rho[1]])
                pyplot.plot(self.z[1:], self.rho[1:])
                pyplot.plot([self.z[-1],self.z[-1]+self.d[-1]], [self.rho[-1],self.rho[-1]])
            pyplot.xlim(xmin=zmin,xmax=zmax)
            pyplot.title('Density')
            pyplot.ylabel("rho (g/cm^3)")
            pyplot.xlabel("z (angstroms)")
        if ty == 'comp':
            ncomp = len(self.comp)
            pyplot.subplot(ncomp,1,1)
            for j in range(ncomp):
                pyplot.subplot(ncomp,1,j+1)
                if bar == True:
                    pyplot.bar(self.z[0],    self.CX[j][0],    width= num.abs(self.d[0]))
                    pyplot.bar(self.z[1:-1], self.CX[j][1:-1], width= num.abs(self.d[1:-1]))
                    pyplot.bar(self.z[-1],   self.CX[j][-1],   width= num.abs(self.d[-1]))
                else:
                    pyplot.plot([self.z[0],0.0], [self.CX[j][0],self.CX[j][0]], hold=hold)
                    pyplot.plot([0.0,self.z[1]], [self.CX[j][0],self.CX[j][1]])
                    pyplot.plot(self.z[1:], self.CX[j][1:])
                    pyplot.plot([self.z[-1],self.z[-1]+self.d[-1]], [self.CX[j][-1],self.CX[j][-1]])
                pyplot.title(self.comp[j].name)
                pyplot.xlim(xmin=zmin,xmax=zmax)
                pyplot.ylabel("mole/cm^3")
            pyplot.xlabel("z (angstroms)")
        if ty == 'el':
            nelem = len(self.elem)
            pyplot.subplot(nelem,1,1)
            for j in range(nelem):
                pyplot.subplot(nelem,1,j+1)
                if bar == True:
                    pyplot.bar(self.z[0],    self.CZ[j][0],    width= num.abs(self.d[0]))
                    pyplot.bar(self.z[1:-1], self.CZ[j][1:-1], width= num.abs(self.d[1:-1]))
                    pyplot.bar(self.z[-1],   self.CZ[j][-1],   width= num.abs(self.d[-1]))
                else:
                    pyplot.plot([self.z[0],0.0], [self.CZ[j][0],self.CZ[j][0]], hold=hold)
                    pyplot.plot([0.0,self.z[1]], [self.CZ[j][0],self.CZ[j][1]])
                    pyplot.plot(self.z[1:], self.CZ[j][1:])
                    pyplot.plot([self.z[-1],self.z[-1]+self.d[-1]], [self.CZ[j][-1],self.CZ[j][-1]])
                pyplot.title(self.elem[j])
                pyplot.xlim(xmin=zmin,xmax=zmax)
                pyplot.ylabel("mole/cm^3")
            pyplot.xlabel("z (angstroms)")
        if ty == 'frac':
            nelem = len(self.elem)
            pyplot.subplot(nelem,1,1)
            for j in range(nelem):
                pyplot.subplot(nelem,1,j+1)
                if bar == True:
                    pyplot.bar(self.z[0],    self.fZ[j][0],    width= num.abs(self.d[0]))
                    pyplot.bar(self.z[1:-1], self.fZ[j][1:-1], width= num.abs(self.d[1:-1]))
                    pyplot.bar(self.z[-1],   self.fZ[j][-1],   width= num.abs(self.d[-1]))
                else:
                    pyplot.plot([self.z[0],0.0], [self.fZ[j][0],self.fZ[j][0]], hold=hold)
                    pyplot.plot([0.0,self.z[1]], [self.fZ[j][0],self.fZ[j][1]])
                    pyplot.plot(self.z[1:], self.fZ[j][1:])
                    pyplot.plot([self.z[-1],self.z[-1]+self.d[-1]], [self.fZ[j][-1],self.fZ[j][-1]])
                pyplot.title(self.elem[j])
                pyplot.xlim(xmin=zmin,xmax=zmax)
                pyplot.ylabel("mole frac")
            pyplot.xlabel("z (angstroms)")

###########################################################################
class Model:
    """
    Layer model
    """
    #_TOP = Layer(comp=[(compound.Component(formula={'N':2}),1)],
    #             density=0.001,thickness=1000.)
    #_SUBS = Layer(comp=[(compound.Component(formula={'Si':1,'O':2}),1)],
    #              density=2.65,thickness=1000.)
    #######################################################################
    def __init__(self,substrate=None,layers=[],top=None,theta=[],params={}):
        """
        substrate and top are Layer instances
        layers is a list of Layer instances
        """
        self.ref    = None
        self.slab   = None
        self.layer  = []
        if substrate != None:
            self.layer.append(substrate)
        for m in layers:
            self.layer.append(m)
        if top != None:
            self.layer.append(top)
        
        # Make sure we have some layers...
        if len(self.layer) == 0:
            self.layer.append(_SUBS)
            self.layer.append(_TOP)
        elif len(self.layer) == 1:
            self.layer.append(_TOP)

        # theta
        self.theta  = num.array(theta,dtype=num.double)

        # params
        self.params = {}
        if len(params) > 0:
            self.set_param(**params)
        
        # Flag
        self._initR = True
        
    #######################################################################
    def __repr__(self,):
        lout = "==== Interface Model ====\n"
        idx = range(len(self.layer))
        idx.reverse()
        for j in idx:
            m = self.layer[j]
            lout = "%s** Layer %d  = %s\n"     %  (lout, j, m.__repr__())
            lout = "%s   d=%g, rho=%g, sig=%g\n\n" % (lout,
                                                      m.thickness,
                                                      m.density,
                                                      m.roughness)
        if self.slab:
            lout = lout + self.slab.__repr__()
        return lout
    
    ########################################################################
    def __save__(self,):
        """
        Make object pickleable...
        """
        #print "running model save"
        del self.ref
        self.ref = None
        self._initR = True

    #######################################################################
    def _init_ref(self,):
        """
        Generate arrays for Reflectivity/FY
        """
        if self.slab == None:
            print "Please slabify the model"
            return
        self.ref = RefModel(d=self.slab.d,
                            rho=self.slab.rho,
                            sigma=self.slab.sig,
                            comp=self.slab.fZ,
                            elem_z=self.slab.elem_z, 
                            theta=self.theta, 
                            params=self.params)
        self._initR = False

    #######################################################################
    def slabify(self,delta=10.):
        """
        Generate a 'slab' model from layers
        use delta <= 0 for just slabs...
        """
        self.slab = Slab(self.layer,delta=delta)
        self._initR = True

    #######################################################################
    def calc_dist(self):
        """
        (re)compute the component distributions
        """
        self.slab.calc_dist()
        # May need to use below if losing reference
        #if self.ref:
        #    self.ref.init_model(d=self.slab.d,
        #                        rho=self.slab.rho,
        #                        comp=self.slab.fZ)

    #######################################################################
    def calc_R(self,energy=None):
        """
        Calc reflectivity at given energy
        """
        if energy!=None:
            self.set_param(energy=energy)
        if self._initR:
            self._init_ref()
        self.ref.calc_R()

    #######################################################################
    def calc_FY(self,energy=None,fyel=None,fyenergy=None):
        """
        Calc reflectivity and FY at given energy
        for the given element and given fyenergy
        """
        if energy!=None or fyel!=None or fyenergy!=None:
            self.set_param(energy=energy,fyel=fyel,fyenergy=fyenergy)
        if self._initR:
            self._init_ref()
        self.ref.calc_FY()

    #######################################################################
    #def set_param(self,energy=None,el=None,fyenergy=None):
    def set_theta(self,theta):
        """
        Set theta
        """
        self.theta  = num.array(theta,dtype=num.double)
        self._initR = True
    
    #######################################################################
    #def set_param(self,energy=None,el=None,fyenergy=None):
    def set_param(self,**params):
        """
        Set parameter vals
        """
        #
        if params.has_key('energy'):
            energy = params.pop('energy')
            if energy != None:
                energy = float(energy)
                self.params['energy'] = energy
                if self.ref: self.ref.set_params(energy=energy)
        #
        if params.has_key('fyel'):
            fyel = params.pop('fyel')
            if fyel != None:
                fyidx = -1.
                #
                if type(fyel) == types.StringType:
                    self.params['fyel']  = fyel.strip()
                    fyel = elements.number(self.params['fyel'])
                else:
                    self.params['fyel']  = fyel
                #
                if (fyel > 0) and (self.slab != None):
                    nelem = len(self.slab.elem_z)
                    for j in range(nelem):
                        if fyel == self.slab.elem_z[j]:
                            fyidx = j
                            break
                    if fyidx == -1.:
                        print "Warning element %s not found" % fyel
                    fyidx = float(fyidx)
                else:
                    fyidx = -1.0
                self.params['fyidx'] = fyidx
                if self.ref: self.ref.set_params(fyidx=fyidx)
        #
        if params.has_key('fyenergy'):
            fyenergy = params.pop('fyenergy')
            if fyenergy != None:
                if type(fyenergy) == types.StringType:
                    fyenergy = xrf_lookup.lookup_xrf_line(fyenergy)
                    if fyenergy == None:
                        print "Error getting fyenergy %s " % fyenergy
                        return
                    else:
                        fyenergy = 1000.*fyenergy
                fyenergy = float(fyenergy)
                self.params['fyenergy'] = fyenergy
                if self.ref: self.ref.set_params(fyenergy=fyenergy)
        #
        if len(params) > 0:
            self.params.update(params)
            if self.ref: self.ref.set_params(**params)

    def get_params(self):
        """
        Return a dictionary with all the params
        """
        params = {}
        if self.ref:
            params = self.ref.get_params()
        else:
            params = DEFAULT_PARAMS
        params.update(self.params)
        return params

    #######################################################################
    def plot(self,ty='calc',fignum=True,bar=True,hold=False,FYflag=True,
             Rdat=None,FYdat=None):
        """
        ty = 'calc' for R and/ or FY calc
        ty = 'density' for density
        ty = 'comp' for component concentration
        ty = 'el' for element concentration
        ty = 'frac' for element fractions
        note hold is only useful for 'calc' plots
        """
        from matplotlib import pyplot
        if ty == 'calc':
            if fignum: pyplot.figure(1)
            self._plot_calc(hold=hold,FYflag=FYflag,Rdat=Rdat,FYdat=FYdat)
        else:
            if fignum:
                if ty == 'density': pyplot.figure(2)
                if ty == 'comp': pyplot.figure(3)
                if ty == 'el': pyplot.figure(4)
                if ty == 'frac': pyplot.figure(5)
            self.slab.plot(ty=ty,bar=bar,hold=hold)

    def _plot_calc(self,hold=False,FYflag=True,Rdat=None,FYdat=None):
        from matplotlib import pyplot
        if hold == False:
            pyplot.clf()
        fyidx = int(self.ref.calc_params[6])
        if  (fyidx >= 0) and (FYflag== True):
            pyplot.subplot(2,1,1)
            pyplot.semilogy()
            pyplot.plot(self.ref.theta,self.ref.R)
            pyplot.grid()
            title = 'Reflectivity'
            pyplot.title(title)
            if Rdat != None:
                try:
                    pyplot.plot(self.ref.theta,Rdat,'.')
                except:
                    print "Error plotting R data"
            #
            pyplot.subplot(2,1,2)
            pyplot.plot(self.ref.theta,self.ref.Y)
            pyplot.grid()
            title = '%s Fluorescent Yield' % self.slab.elem[fyidx]
            pyplot.title(title)
            if FYdat != None:
                try:
                    pyplot.plot(self.ref.theta,FYdat,'.')
                except:
                    print "Error plotting FY data"
            #
            pyplot.xlabel("theta (deg)")
        else:
            pyplot.subplot(1,1,1)
            pyplot.semilogy()
            pyplot.plot(self.ref.theta,self.ref.R)
            pyplot.grid()
            title = 'Reflectivity'
            pyplot.title(title)
            if Rdat != None:
                try:
                    pyplot.plot(self.ref.theta,Rdat,'.')
                except:
                    print "Error plotting R data"
            #
            pyplot.xlabel("theta (deg)")
            
############################################################################
############################################################################
############################################################################
def test_model():
    """
    A test case
    """    
    import time
    show_time = True

    ### build some chemical components
    N2     = compound.Component(formula={'N':2})
    qtz    = compound.Component(formula={'Si':1,'O':2})
    fe2o3  = compound.Component(formula={'Fe':2,'O':3})
    al2o3  = compound.Component(formula={'Al':2,'O':3})
    
    ### build some layers
    subs   = Layer(comp=[(qtz,1.),(fe2o3,0.000001)],density=2.65,thickness=1000.,roughness=10.)
    m1     = Layer(comp=[(qtz,1.),(fe2o3,0.0001)],density=2.45,thickness=30.,roughness=10.)
    m2     = Layer(comp=[(N2,1.)],density=0.001,thickness=1000.,roughness=0.)
    top    = Layer(comp=[(N2,1.)],density=0.001,thickness=1000.,roughness=0.)

    ### parameters
    delta = 5.
    theta = num.arange(0.01, 1.0, 0.01)
    calc_params = {'energy':10000.,'wconv':0.01,'slen':20.,'bvert':0.01,
                   'aflag':1.,'fyidx':0,'fyenergy':7000.,'delz':delta,'pdepth':2.0}

    ###################################
    ### create model
    t = time.time()
    model = Model(substrate=subs,layers=[m1,m2],top=top, theta=theta, params=calc_params)
    model.slabify(delta=delta)
    if show_time: print " Elapsed time to slabify = %f seconds\n" % (time.time() - t)

    ###################################    
    ### FY/R on init'd model
    t = time.time()
    #model.set_param(fyel='Fe',fyenergy=7000.)
    #model.calc_FY()
    #model.calc_FY(fyel='Fe',fyenergy=7000.)
    model.calc_FY(fyel='Fe',fyenergy='Fe Ka')
    if show_time: print " Elapsed time to calc FY = %f seconds\n" % (time.time() - t)
    model.plot(ty='calc')
    #model.plot(ty='density')
    #model.plot(ty='comp')
    #model.plot(ty='el')
    #model.plot(ty='frac')
    model.slab.comp_mb()

    #return model

    ###################################
    ### change dist and recompute FY/R
    raw_input('**Enter to continue:')
    #
    dist = {'type':'erf','cen':40.,'sig':20.,'CX':1.}
    model.slab.add_dpar('Fe2_O3',dist=dist,norm=1,init=True)
    #
    dist = {'type':'gauss','cen':70.,'sig':10.,'CX':30.}
    model.slab.add_dpar('Fe2_O3',dist=dist)
    #
    dist = {'type':'erf','cen':40.,'sig':20.,'CX':1.}
    model.slab.add_dpar('Si1_O2',dist=dist,norm=1,init=True)
    #
    t = time.time()
    model.slab.calc_dist()
    if show_time: print " Elapsed time to calc dist = %f seconds\n" % (time.time() - t)
    t = time.time()
    model.calc_FY()
    if show_time: print " Elapsed time to calc FY = %f seconds\n" % (time.time() - t)
    model.plot(ty='calc',hold=True)
    model.plot(ty='density')
    model.plot(ty='comp')
    model.slab.comp_mb()

    #################################
    ### mod dist and recompute FY/R
    raw_input('** Enter to continue:')
    dist = {'cen':60.,'sig':40.,'CX':1.}
    model.slab.set_dpar('Fe2_O3',dist=dist,didx=0,norm=1)
    #
    t = time.time()
    model.slab.calc_dist()
    if show_time: print " Elapsed time to calc dist = %f seconds\n" % (time.time() - t)
    t = time.time()
    model.calc_FY()
    if show_time: print " Elapsed time to calc FY = %f seconds\n" % (time.time() - t)
    model.plot(ty='calc',hold=True)
    model.plot(ty='density')
    model.plot(ty='comp')
    model.slab.comp_mb()
    
    return model


############################################################################
if __name__ == "__main__":
    model = test_model()
    
