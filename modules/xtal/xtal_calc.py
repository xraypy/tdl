##########################################################################
"""
Tom Trainor (tptrainor@alaska.edu) 
Xtal calcs

Modifications:
--------------
- convert matlab to python, Kunal Tanwar
- Make object oriented, TPT

"""
##########################################################################
"""
Todo
 - Lattice:
    - add cell volume calcs
    - add real_to_recip and recip_to_real (ie get vector indicies
      in different basis)
    - calculate angles between real and recip vectors "angle_rr"

 - LatticeTransform class 

 - Reading cif files and others (maybe seperate module)
   - e.g. get xyz file from fit and par files

 - Structure analysis and bond valence calcs (seperate module)

 - Make another module with space groups and typical lattice transformations
   from the international tables...
 
"""
##########################################################################

import numpy as num
import sys
import types

##########################################################################
class Lattice:
    """
    Class that defines a lattice and various
    operations within the lattice
    """
    def __init__(self,a=1.,b=1.,c=1.,alpha=90.,beta=90.,gamma=90.):
        self.update(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma)

    def __repr__(self,):
        lout = "a=%6.5f, b=%6.5f, c=%6.5f" % (self.a, self.b, self.c)
        lout = "%s, alpha=%6.5f,beta=%6.5f,gamma=%6.5f\n" % (lout,
                                                             self.alpha,
                                                             self.beta,
                                                             self.gamma)
        lout = "%sar=%6.5f, br=%6.5f, cr=%6.5f" % (lout,self.ar, self.br, self.cr)
        lout = "%s, alphar=%6.5f,betar=%6.5f,gammar=%6.5f" % (lout,
                                                              self.alphar,
                                                              self.betar,
                                                              self.gammar)
        return lout

    def update(self,a=None,b=None,c=None,alpha=None,beta=None,gamma=None):
        """
        update lattice parameters
        """
        if a != None: self.a = float(a)
        if b != None: self.b = float(b)
        if c != None: self.c = float(c)
        if alpha != None: self.alpha = float(alpha)
        if beta  != None: self.beta  = float(beta)
        if gamma != None: self.gamma = float(gamma)
        # update calc quantities
        self._calc_g()

    def _calc_g(self):
        """
        calculate the metric tensors and recip lattice params
           self.g  = real space metric tensor
           self.gr = recip space metric tensor
        """
        a = self.a
        b = self.b
        c = self.c
        alp = num.radians(self.alpha)
        bet = num.radians(self.beta)
        gam = num.radians(self.gamma)
        # real metric tensor
        self.g = num.array([ [ a*a, a*b*num.cos(gam), a*c*num.cos(bet) ],
                             [ b*a*num.cos(gam), b*b, b*c*num.cos(alp) ],
                             [ c*a*num.cos(bet), c*b*num.cos(alp), c*c ] ])
        # recip lattice metric tensor
        # and recip lattice params
        self.gr   = num.linalg.inv(self.g)
        self.ar   = num.sqrt(self.gr[0,0])
        self.br   = num.sqrt(self.gr[1,1])
        self.cr   = num.sqrt(self.gr[2,2])
        self.alphar = num.degrees(num.arccos(self.gr[1,2]/(self.br*self.cr)))
        self.betar  = num.degrees(num.arccos(self.gr[0,2]/(self.ar*self.cr)))
        self.gammar = num.degrees(num.arccos(self.gr[0,1]/(self.ar*self.br)))
        
    def cell(self):
        """
        return array of real lattice cell parameters
        """
        return num.array([self.a,self.b,self.c,
                          self.alpha,self.beta,self.gamma],dtype=float)

    def rcell(self):
        """
        return array of reciprocal lattice cell parameters
        """
        return num.array([self.ar,self.br,self.cr,
                          self.alphar,self.betar,self.gammar],dtype=float)
        
    def V(self,recip=False):
        """
        Calculate the cell volume.
          If recip == False, this is the real space vol in ang**3
          If recip == True, this is the recip space vol in ang**-3
        """
        if recip == True: g = self.gr
        else: g = self.g
        det = num.linalg.det(g)
        if det > 0:
            return num.sqrt(det)
        else:
            return 0.0
        
    def dot(self,u,v,recip=False):
        """
        calculate dot product of two vectors
        if u and v are real space vectors, recip = False
        if u and v are recip space vectors, recip = True
        Note u and v are assumed to be normal numpy arrays
        (ie not matrix objects)
        """
        if recip == True: g = self.gr
        else: g = self.g
        u = num.array(u,dtype=float)
        v = num.array(v,dtype=float)
        # dot product = u*g*v, from Sands.
        return num.dot(u,num.dot(g,v))

    def mag(self,v,recip=False):
        """
        calculate the norm of a vector
        if v is real space vector, recip = False
        if v is recip space vector, recip = True
        Note v is assumed to be normal numpy array
        (ie not matrix objects)
        """
        m = num.sqrt(self.dot(v,v,recip=recip))
        return m

    def angle(self,u,v,recip=False):
        """
        calculate dot product of two vectors
        if u and v are real space vectors, recip = False
        if u and v are recip space vectors, recip = True
        Note u and v are assumed to be normal numpy arrays
        (ie not matrix objects)
        """
        uv = self.dot(u,v,recip=recip)
        um = self.mag(u,recip=recip)
        vm = self.mag(v,recip=recip)
        arg = uv/(um*vm)
        if num.fabs(arg) > 1.0:
            arg = arg / num.fabs(arg)
        alpha = num.arccos(arg)
        return num.degrees(alpha)

    def d(self,hkl):
        """
        calculate d space for given [h,k,l]
        """
        if len(hkl)!=3:
            print "Error need an array of [h,k,l]"
            return 0.
        H = self.mag(hkl,recip=True)
        if H == 0.:
            print "Error [h,k,l] magnitude is zero:", hkl
            return 0.0
        d = 1./H
        return d

    def tth(self,hkl,lam=1.5406):
        """
        calculate 2Theta for given [h,k,l] and lambda
        default lambda = Cu Ka1
        """
        d = self.d(hkl)
        if d == 0.: return 0.
        r = lam/(2.*d)
        if num.fabs(r) > 1.0:
            r = r/num.fabs(r)
        tth = 2.*num.arcsin(r)
        tth = num.degrees(tth)
        return tth

    def dvec(self,hkl):
        """
        calculate the real space vector d
        which has a magnitude of the d spacing
        and is normal to the plane HKL
        """
        dspc = self.d(hkl)
        if dspc == 0: return num.array([0.,0.,0.])
        # convert hkl vector to real space indicies
        hkl  = num.array(hkl,dtype=float)
        dvec = num.dot(self.gr,hkl)
        dvec = (dspc**2.)*dvec
        return dvec



##########################################################################
"""
Note make a new class:
class LatticeTransform:
    def __init__(lattice,)
    def cartesian()
      - calculates new cartesian basis
    def new_basis(Va,Vb,Vc,shift)
      - takes set of vectors in lattice basis that
        defines a new basis
    def real_indicies(uvw):
      - takes uvw in lattice and gives same vector in new basis
    def recip_indicies
      - same as above except in recip space
    need same as above 2 functions but going in opposite direction,
    ie transform from new indicies to old indicies
"""
def basis_transform_cart(cell):
    """
    * basis_transform_cart(cell) :
      calculate basis and coordinate transformation matrix i.e. Gc and Fc to
      transform to cartesian space from a basis defined by cell
    """
    #calculate Gc and Fc for transformations to cartesian basis
    a = cell[0]
    b = cell[1]
    c = cell[2]
    alp = cell[3]
    bet = cell[4]
    gam = cell[5]
    Pi = num.pi
    rcell = rlat(cell)
    a_r = rcell[0]
    b_r = rcell[1]
    c_r = rcell[2]
    alp_r = rcell[3]
    bet_r = rcell[4]
    gam_r = rcell[5]
    # Gc transforms basis to cartesian
    Gc = num.array([[1/a, 0 , 0],
                     [-1/(a*num.tan(gam*Pi/180)), 1/(b*num.sin(gam*Pi/180)), 0],
                     [a_r*num.cos(bet_r*Pi/180), b_r*num.cos(alp_r*Pi/180), c_r]])
    # F transforms coordinates to cartesian
    # F = inverse(G')
    Fc = num.array([[a, b*num.cos(gam*Pi/180), c*num.cos(bet*Pi/180)],
                     [0, b*num.sin(gam*Pi/180), -1*c*num.sin(bet*Pi/180)*num.cos(alp_r*Pi/180)],
                     [0, 0, 1/c_r]])
    return Fc

##########################################################################
def transform_to_cart(u,v,w,cell):
    """
    
    * transform_to_cart(u,v,w,cell) :
      transform given u,v,w coordinates from basis defined by cell to cartesian basis

    """
    # transform u,v,w to cartesian
    Fc = basis_transform_cart(cell) 
    uvw = [u, v, w]
    uvw_xyz = num.matrixmultiply(Fc, num.transpose(uvw))
    xyz = num.transpose(uvw_xyz1)
    return xyz

##########################################################################
def cart_to_cell(x,y,z,cell):
    """
    * cart_to_cell(x,y,z,cell) :
      transform cartesian x,y,z coordinates to a basis defined by given cell 
    """
    # transform x,y,z to basis defined by cell
    Fc = basis_transform_cart(cell)
    xyz = [x, y, z]
    xyz_uvw = num.matrixmultiply((num.linalg.inverse(Fc)),num.transpose(xyz))
    uvw = num.transpose(xyz_uvw)
    return uvw

##########################################################################
# some lattice transformations
"""
should make a dictionary/ data structure that defines these in terms
of set of lattice vectors / tranlations that defines new in terms of old
then use to calc a new Lattice...
See international tables...
"""
##########################################################################
def trans_hexa_to_rhombo(Va):
    """
    * trans_hexa_to_rhombo(Va):
      transfrom any vector Va from hexa to rhombo system
    """
    # transfrom any vector from hexa to rhombo system
    MM = num.array([[0.6667, 0.3333, 0.3333],
                      [-0.3333, 0.3333, 0.3333],
                      [-0.3333, -0.6667, 0.3333]])
    NN = num.linalg.inv(num.transpose(MM))
    V_rh = num.matrixmultiply(NN,num.transpose(Va))
    return V_rh

##########################################################################
def trans_rhombo_to_hexa(Va):
    """
    * trans_rhombo_to_hexa(Va):
      transfrom any vector from rhombo to hexa system
    """
    # transfrom any vector from rhombo to hexa system
    MM = num.array([[0.6667, 0.3333, 0.3333],
                      [-0.3333, 0.3333, 0.3333],
                      [-0.3333, -0.6667, 0.3333]])
    NN = num.linalg.inv(num.transpose(MM))
    NN_inv = num.linalg.inv(NN)
    V_hx = num.matrixmultiply(NN_inv,num.transpose(Va))
    return V_hx
    
##########################################################################


##########################################################################
# generators
##########################################################################
class PositionGenerator:
    """
    class to generate atomic positions given symmetry operations
    """
    def __init__(self):
        self.ops = []

    ###########################################################    
    def add_op(self,sym='x,y,z',shift=''):
        """
        add a new symmetry operator for generating positions
        sym and shift are strings with comma delimeted set of
        characters that defines the opertions.  e.g.
          sym = "x,y,z", shift = "0,0,0"
          sym = "-y,z,x+y", shift = "1/2,1/2,0"
        """
        #check shift
        if len(shift) > 0:
            shifts = shift.split(',')
            if len(shifts) != 3:
                print "Error parsing shift, should have 3 components: ", shift
                return None
        else:
            shifts = ['','','']
        # break up sym into x,y,z parts
        syms = sym.split(',')
        if len(syms) != 3:
            print "Error parsing sym, should have 3 components: ", shift
            return None

        x = syms[0] + shifts[0]
        y = syms[1] + shifts[1]
        z = syms[2] + shifts[2]
        #print 'x=',x,'y=',y,'z=',z
        
        m = self._make_seitz_matrix(x,y,z)
        if m != None:
            self.ops.append(m)

    ###########################################################    
    def _make_seitz_matrix(self,x,y,z):
        """
        generate augmented (seitz) matrix given
        string symbols for x,y,z, coordinates
        """
        def _vec(sym):
            v = num.array([0.,0.,0.,0.])
            if type(sym) != types.StringType:
                print "Error, passed a non-string symbol"
                return None
            sym = sym.replace('+','')
            sym = sym.replace(' ','')
            if '-x' in sym:
                v[0] = -1
                sym = sym.replace('-x','')
            elif 'x' in sym:
                v[0] = 1
                sym = sym.replace('x','')
            if '-y' in sym:
                v[1] = -1
                sym = sym.replace('-y','')
            elif 'y' in sym:
                v[1] = 1
                sym = sym.replace('y','')
            if '-z' in sym:
                v[2] = -1
                sym = sym.replace('-z','')
            elif 'z' in sym:
                v[2] = 1
                sym = sym.replace('z','')
            if len(sym) > 0:
                sym = sym + '.'
                v[3] = eval(sym)
            return v
        #
        v1 = _vec(x)
        if v1 == None: return None
        v2 = _vec(y)
        if v2 == None: return None
        v3 = _vec(z)
        if v3 == None: return None
        v4 = [0.,0.,0.,1.]
        m = num.array([v1,v2,v3,v4])
        #print m
        return m

    ###########################################################    
    def copy(self,x,y,z,reduce=True,rem_dups=True):
        """
        given x,y,z coordinates, calc all sym copies
        """
        v0 = num.array([float(x),float(y),float(z),1.0])
        vectors = []
        for m in self.ops:
            vc = num.dot(m,v0)
            vectors.append(vc[:3])
        
        # reduce values
        def _reduce(x):
            while 1:
                if x >= 1.0:
                    x = x-1.0
                if x < 0.0:
                    x = x + 1.0
                if num.fabs(x) < 1.0:
                    return x
        #
        if reduce == True:
            for j in range(len(vectors)):
                if vectors[j][0] >= 1.0 or vectors[j][0] < 0.0:
                    vectors[j][0] = _reduce(vectors[j][0])
                if vectors[j][1] >= 1.0 or vectors[j][1] < 0.0:
                    vectors[j][1] = _reduce(vectors[j][1])
                if vectors[j][2] >= 1.0 or vectors[j][2] < 0.0:
                    vectors[j][2] = _reduce(vectors[j][2])

        # remove duplicates
        def _rem_dups(vectors):
            unique = []
            while len(vectors) > 0:
                v = vectors.pop(0)
                #print "-->new vector=", v
                add = True
                for u in unique:
                    #print '    ', u
                    tmp = num.equal(v,u)
                    if tmp.all():
                        add = False
                        break
                if add:
                    unique.append(v)
            return unique

        if rem_dups == True:
            vectors = _rem_dups(vectors)

        return vectors

##########################################################################
##########################################################################
def _test1():
    """
    C2/m
    """
    sym1 = "x,y,z"
    sym2 = "-x,y,-z"
    sym3 = "-x,-y,-z"
    sym4 = "x,-y,z"
    shift0 = "0,0,0"
    shift1 = "1/2, 1/2, 0"
    p = PositionGenerator()
    p.add_op(sym=sym1,shift=shift0)
    p.add_op(sym=sym2,shift=shift0)
    p.add_op(sym=sym3,shift=shift0)
    p.add_op(sym=sym4,shift=shift0)
    p.add_op(sym=sym1,shift=shift1)
    p.add_op(sym=sym2,shift=shift1)
    p.add_op(sym=sym3,shift=shift1)
    p.add_op(sym=sym4,shift=shift1)
    return p

##########################################################################
if __name__ == "__main__":
    """
    test PositionGenerator
    for C2/m
    """
    p = _test1()
    print "0.15,0.0,0.33"
    vecs = p.copy(0.15,0.0,0.33,reduce=True,rem_dups=True)
    for v in vecs: print v

    print "0.5,0.11,0.5"
    vecs = p.copy(0.5,0.11,0.5,reduce=True,rem_dups=True)
    for v in vecs: print v

    print "0.25,0.25,0.25"
    vecs = p.copy(0.25,0.25,0.25,reduce=True,rem_dups=True)
    for v in vecs: print v
