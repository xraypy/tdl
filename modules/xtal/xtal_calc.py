##########################################################################
"""
Tom Trainor (tptrainor@alaska.edu) 
Xtal calcs

Modifications:
--------------
- convert matlab to python, Kunal Tanwar
- Make into class, TPT

"""
##########################################################################
"""
Todo

 - LatticeTransform class 

 - Reading cif files and others (maybe seperate module)
   - e.g. get xyz file from fit and par files

 - Structure analysis and bond valence calcs (seperate module)

 - Make another module with space groups and typical lattice transformations
   from the international tables...
 
"""
##########################################################################
"""
Notes on Convention and Generalized Transforms
----------------------------------------------
The set of basis vectors that defines the real lattice are
covarient quantities, while the set of indicies that define
real space vectors are contravarient:
                 |a|
    v = (x,y,z)* |b| = x*a + y*b + z*c
                 |c|
If A is an arbitrary matrix that transforms the contravarient
vector indicies, its assumed that the indicies left multiply 
    (x',y',z') = (x,y,z) * A

Note that the usual convention is to express covarient quantities
as column vectors and contravarient quantities as row vectors
(subscripts and superscripts respectivley in tensor notation)
However, in the notes below we'll just take 1-D arrays to be
column vectors to simplify the notation.  Therefore, if A is defined
as above, the matrix vector multiplications will need to be changed to:
    |x'|                 |x|
    |y'| = transpose(A)* |y|
    |z'|                 |z|
We can express the same in numpy given a vector v and matrix A 
    dot(v,A) = dot(transpose(A),v)
where v is multiplied as a row vector on the lhs and a column
vector on the rhs; ie the transpose of v bewteen column and
row vecor is implied by its position.

Assume there is a matrix that transforms the covarient basis
(via rotation and/or dilation) to a new set of basis vectors:
    |a'|      |a|        a' = (F11*a + F12*b + F13*c)
    |b'| = F* |b|   or   b' = (F21*a + F22*b + F23*c)
    |c'|      |c|        c' = (F31*a + F32*b + F33*c)
with a well defined inverse:
    |a|           |a'|   
    |b| = inv(F)* |b'|   
    |c|           |c'|   

The covarient/contravarient relationship bewteen basis vectors
and the vector indicies implies that the indicies of a (stationary)
vector in the new (') basis are obtained via the inverse of the
F matrix:
    |x'|      |x|
    |y'| = M* |y|
    |z'|      |z|
where
    M = transpose(inv(F)) = inv(transpose(F))

The reciprocal relation is:
    |x|      |x'|
    |y| = N* |y'|
    |z|      |z'|
where
    N = inv(M) = transpose(F)

With respect to the real lattice the recip lattice basis vectors
are contravarient and the vector indicies are covarient:
                     |h|
    vr = (ar,br,cr)* |k| = h*ar + k*br + l*cr
                     |l|

Given the generalized basis transform F above, the covarient recip
lattice indicies transform in the same way as the lattice.
Therefore:
     |h'|      |h|
     |k'| = F* |k|
     |l'|      |l|
     
The reciprocal relation is:
    |h|      |k'|
    |k| = G* |k'|
    |l|      |l'|
where
    G = inv(F) = transpose(M)


The metric tensor
------------------
The metric tensor is defined in terms of the vector dot product
for a general real space system
  u*v = u*g*v
where
       |a*a, a*b, a*c|
  g =  |b*a, b*b, b*c|
       |c*a, c*b, c*c|
The metric tensor may also be defined in terms of a lattice transform
ie the real space basis vectors can be expressed in terms of the
recip space basis vectors via
      |a|      |ar|
      |b| = g* |br|
      |c|      |cr|
therefore
      |ar|       |a|
      |br| = gr* |b|
      |cr|       |c|
where gr = inv(g).  Hence the recip lattice vectors
can be plotted in the real lattice via
      ar = gr11*a + gr12*b + gr13*c
      br = gr21*a + gr22*b + gr23*c
      cr = gr31*a + gr32*b + gr33*c
and visa versa

Using the above relations for lattice transforms the
indicies of real lattice vectors can be calc in the recip
lattice and visa-versa via:
   |h|                             |x|      |x|
   |k| = (x,y,z)*g = transpose(g)* |y| = M* |y|
   |l|                             |z|      |z|
and
   |x|                      |h|      |h|
   |y| = inv(transpose(g))* |k| = N* |k|
   |z|                      |l|      |l|

e.g. in above we take recip lattice as primed basis and
real lattice as the unprimed basis, therefore equate:
  F = gr
  G = inv(F) = inv(gr) = g 
  M = transpose(inv(F)) = transpose(g)
  N = inv(M) = transpose(F) = inv(transpose(g)) = transpose(gr)

Note that g and gr are always symmetric, therefore
  g = transpose(g) and gr = transpose(gr).
  F = N = gr
  G = M = g
  
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
        Update lattice parameters
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
        Calculate the metric tensors and recip lattice params
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
        Return array of real lattice cell parameters
        """
        return num.array([self.a,self.b,self.c,
                          self.alpha,self.beta,self.gamma],dtype=float)

    def rcell(self):
        """
        Return array of reciprocal lattice cell parameters
        """
        return num.array([self.ar,self.br,self.cr,
                          self.alphar,self.betar,self.gammar],dtype=float)
        
    def vol(self,recip=False):
        """
        Calculate the cell volume.
          If recip == False, this is the real space vol in ang**3
          If recip == True, this is the recip space vol in ang**-3

        V_real  = sqrt(determiant(g))
        V_recip = sqrt(determiant(inv(g)))

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
        Calculate dot product of two vectors
          If u and v are real space vectors, use recip = False
          If u and v are recip space vectors, use recip = True
        Note u and v are assumed to be normal numpy arrays
        (ie not matrix objects)

        The generalized dot product for real space vectors
        may be defined in terms of the metric tensor (g):
           u*v = u*g*v
        If u and v are recip lattice vectors replace g with
        gr = inv(g)
        """
        if recip == True: g = self.gr
        else: g = self.g
        u = num.array(u,dtype=float)
        v = num.array(v,dtype=float)
        return num.dot(u,num.dot(g,v))

    def mag(self,v,recip=False):
        """
        Calculate the norm of a vector
          If v is real space vector, use recip = False
          If v is recip space vector, use recip = True
        Note v is assumed to be normal numpy array
        (ie not matrix objects)
        """
        m = num.sqrt(self.dot(v,v,recip=recip))
        return m

    def angle(self,u,v,recip=False):
        """
        Calculate dot product of two vectors
          If u and v are real space vectors, use recip = False
          If u and v are recip space vectors, use recip = True
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

    def angle_rr(self,x,h):
        """
        Calculate the angle between a real vector x = [x,y,z]
        and a recip vector h = [h,k,l]
        """
        x = num.array(x,dtype=float)
        h = num.array(h,dtype=float)
        hx = x*h
        xm = self.mag(x,recip=False)
        hm = self.mag(h,recip=True)
        arg = hx/(hm*xm)
        if num.fabs(arg) > 1.0:
            arg = arg / num.fabs(arg)
        alpha = num.arccos(arg)
        return num.degrees(alpha)

    def d(self,hkl):
        """
        Calculate d space for given [h,k,l]
        """
        if len(hkl)!=3:
            print "Error need an array of [h,k,l]"
            return 0.
        h = self.mag(hkl,recip=True)
        if h == 0.:
            print "Error [h,k,l] magnitude is zero:", hkl
            return 0.0
        d = 1./h
        return d

    def tth(self,hkl,lam=1.5406):
        """
        Calculate 2Theta for given [h,k,l] and lambda
        Default lambda = Cu Ka1
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
        Calculate the real space vector d
        which has a magnitude of the d spacing
        and is normal to the plane hkl = [h,k,l]
        """
        # convert hkl vector to real space indicies
        #hkl  = num.array(hkl,dtype=float)
        #dvec = num.dot(self.gr,hkl)
        dvec = self.recip_to_real(hkl)
        dspc = self.d(hkl)
        if dspc == 0: return num.array([0.,0.,0.])
        dvec = (dspc**2.)*dvec
        return dvec

    def recip_to_real(self,hkl):
        """
        Given recip vector hkl = [h,k,l] calculate the
        vectors indicies in the real lattice
        """
        hkl = num.array(hkl,dtype=float)
        #v   = num.dot(self.gr.transpose(),hkl)
        v   = num.dot(self.gr,hkl)
        return v
    
    def real_to_recip(self,v):
        """
        Given real vector v = [x,y,z] calculate the
        vectors indicies in the reciprocal lattice
        """
        v   = num.array(v,dtype=float)
        hkl = num.dot(v,self.g)
        return hkl

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
of set of lattice vectors / translations that defines new in terms of old
then use to calc a new Lattice...
See international tables for usual / typical transforms
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
if __name__ == "__main__":
    """
    test 
    """
