##########################################################################
"""
Tom Trainor (tptrainor@alaska.edu) 
lattice calcs

Modifications:
--------------
- convert matlab to python, Kunal Tanwar
- Convert classes, TPT

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

from mathutil import cosd, sind, tand
from mathutil import arccosd, arcsind, arctand

##########################################################################
class Lattice:
    """
    Class that defines a lattice and various
    operations within the lattice

    """
    def __init__(self,a=10.,b=10.,c=10.,alpha=90.,beta=90.,gamma=90.,lam=1.5406):
        """
        Initialize by passing a,b,c in angstroms 
        alpha, beta, gamma in degrees,
        and lambda in angstroms (default lambda is Cu Ka1)
        """
        self.update(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,lam=lam)

    def __repr__(self,):
        lout = "a=%6.5f, b=%6.5f, c=%6.5f" % (self.a, self.b, self.c)
        lout = "%s, alpha=%6.5f,beta=%6.5f,gamma=%6.5f\n" % (lout,
                                                             self.alpha,
                                                             self.beta,
                                                             self.gamma)
        lout = "%sar=%6.5f, br=%6.5f, cr=%6.5f" % (lout,self.ar, self.br, self.cr)
        lout = "%s, alphar=%6.5f,betar=%6.5f,gammar=%6.5f\n" % (lout,
                                                                self.alphar,
                                                                self.betar,
                                                                self.gammar)
        lout = "%sDefault wavelength for angle calculations=%6.5f\n" % (lout,
                                                                      self.lam)
        return lout

    def update(self,a=None,b=None,c=None,alpha=None,beta=None,gamma=None,lam=None):
        """
        Update lattice parameters
        """
        if a != None:     self.a = float(a)
        if b != None:     self.b = float(b)
        if c != None:     self.c = float(c)
        if alpha != None: self.alpha = float(alpha)
        if beta  != None: self.beta  = float(beta)
        if gamma != None: self.gamma = float(gamma)
        if lam != None:   self.lam = float(lam)
        # update calc quantities
        self._calc_g()

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

    def _calc_g(self):
        """
        Calculate the metric tensors and recip lattice params
           self.g  = real space metric tensor
           self.gr = recip space metric tensor
        """
        (a,b,c,alp,bet,gam) = self.cell()
        # real metric tensor
        self.g = num.array([ [ a*a, a*b*cosd(gam), a*c*cosd(bet) ],
                             [ b*a*cosd(gam), b*b, b*c*cosd(alp) ],
                             [ c*a*cosd(bet), c*b*cosd(alp), c*c ] ])
        # recip lattice metric tensor
        # and recip lattice params
        self.gr     = num.linalg.inv(self.g)
        self.ar     = num.sqrt(self.gr[0,0])
        self.br     = num.sqrt(self.gr[1,1])
        self.cr     = num.sqrt(self.gr[2,2])
        self.alphar = arccosd(self.gr[1,2]/(self.br*self.cr))
        self.betar  = arccosd(self.gr[0,2]/(self.ar*self.cr))
        self.gammar = arccosd(self.gr[0,1]/(self.ar*self.br))
        
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
        Calculate the angle between two vectors
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
        alpha = arccosd(arg)
        return alpha

    def angle_rr(self,x,h):
        """
        Calculate the angle between a real vector x = [x,y,z]
        and a recip vector h = [h,k,l]
        """
        x = num.array(x,dtype=float)
        h = num.array(h,dtype=float)
        hx = num.sum(x*h)
        xm = self.mag(x,recip=False)
        hm = self.mag(h,recip=True)
        arg = hx/(hm*xm)
        if num.fabs(arg) > 1.0:
            arg = arg / num.fabs(arg)
        alpha = arccosd(arg)
        return alpha

    def d(self,hkl):
        """
        Calculate d space for given [h,k,l]
        """
        if len(hkl)!=3:
            print "Error need an array of [h,k,l]"
            return 0.
        h = self.mag(hkl,recip=True)
        if h == 0.:
            #print "Error [h,k,l] magnitude is zero:", hkl
            return 0.0
        d = 1./h
        return d

    def tth(self,hkl,lam=None):
        """
        Calculate 2Theta for given [h,k,l] and wavelength
        
        If lam is None the default lambda will be used (Cu Ka1)
        If you pass in lam, this will change the default for
        subsequent calls
        """
        if lam != None: self.lam = float(lam)
        d = self.d(hkl)
        if d == 0.: return 0.
        r = self.lam/(2.*d)
        if num.fabs(r) > 1.0:
            r = r/num.fabs(r)
        tth = 2.*arcsind(r)
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
class LatticeTransform:
    """
    Generalized lattice transformations
    """
    def __init__(self,lattice,Va=None,Vb=None,Vc=None,shift=None):
        """
        Initialize with a Lattice instance.
        All transforms are defined with respect to this
        lattice, ie this is the original unprimed lattice
        """
        self.lattice = lattice
        self.Va    = num.array([1.,0.,0.])
        self.Vb    = num.array([0.,1.,0.])
        self.Vc    = num.array([0.,0.,1.])
        self.shift = num.array([0.,0.,0.])
        self._update(Va=Va,Vb=Vb,Vc=Vc,shift=shift)
    
    def basis(self,Va=None,Vb=None,Vc=None,shift=None):
        """
        Define new basis vectors.
        Va, Vb and Vc should define the a',b',c' lattice
        vectors of the new basis (rotation/dialation part).
        Shift describes an origin shift of the new lattice
        (ie take the new basis defined by Va,Vb,Vc then
        apply the shift vector)
        
        All the vectors should be defined in terms of fractional
        coordinate indicies in the original basis, for example:
           a' = x1*a + y1*b + z1*c 
           b' = x2*a + y2*b + z2*c 
           c' = x3*a + y3*b + z3*c
           shift = x4*a + y4*b + z4*c
        and
           Va = [x1,y1,y1]
           Vb = [x2,y2,y2]
           Vc = [x3,y3,y3]
           shift = [x4,y4,z4]
        """
        self._update(Va=Va,Vb=Vb,Vc=Vc,shift=shift)

    def _update(self,Va=None,Vb=None,Vc=None,shift=None):
        if Va != None: self.Va = num.array(Va,dtype=float)
        if Vb != None: self.Vb = num.array(Vb,dtype=float)
        if Vc != None: self.Vc = num.array(Vc,dtype=float)
        if shift != None: self.shift = num.array(shift,dtype=float)
        F = [self.Va,self.Vb,self.Vc]
        self.F = num.array(F, dtype=float)
        self.G = num.linalg.inv(self.F)
        self.M = self.G.transpose()
        self.N = self.F.transpose()
        
    def cartesian(self,shift=[0.,0.,0.]):
        """
        Calculates a cartesian basis using:
          Va = a' is parrallel to a
          Vb = b' is perpendicular to a' and in the a/b plane
          Vc = c' is perpendicular to the a'/c' plane
        A shift vector may be specified to shift the origin
        of the cartesian lattice relative to the original lattice
        origin (specify shift in fractional coordinates of
        the original lattice)
        """
        (a,b,c,alp,bet,gam) = self.lattice.cell()
        (ar,br,cr,alpr,betr,gamr) = self.lattice.rcell()
        Va = [1./a,                 0. ,             0.]
        Vb = [-1./(a*tand(gam)), 1./(b*sind(gam)),   0.]
        Vc = [ar*cosd(betr),     br*cosd(alpr),      cr]
        self.basis(Va=Va,Vb=Vb,Vc=Vc,shift=shift)

    def xp(self,x):
        """
        Given x = [x,y,z] in the original lattice
        compute the indicies of the vector in the new basis
           xp = M*(x-shift)
        """
        x  = num.array(x,dtype=float)
        if self.shift.sum() != 0.:
            x = x - self.shift
        xp = num.dot(self.M,x)
        return xp
        
    def x(self,xp):
        """
        Given xp = [x',y',z'] in the primed basis
        compute the indicies of the vector in the original basis
            x = N*(xp+shiftp) = N*xp + shift
        """
        xp  = num.array(xp,dtype=float)
        x = num.dot(self.N,xp)
        if self.shift.sum() != 0.:
            x = x + self.shift
        return x

    def hp(self,h):
        """
        Given h = [h,k,l] in the original (recip) lattice
        compute the indicies of the vector in the new (recip) basis
        """
        h  = num.array(h,dtype=float)
        return num.dot(self.F,h)
        
    def h(self,hp):
        """
        Given hp = [h',k',l'] in the primed (recip) basis
        compute the indicies of the vector in the original (recip) basis
        """
        hp  = num.array(hp,dtype=float)
        return num.dot(self.G,hp)

    def plat(self):
        """
        Return a Lattice instance for the primed basis 
        """
        a = self.lattice.mag(self.Va)
        b = self.lattice.mag(self.Vb)
        c = self.lattice.mag(self.Vc)
        alp = self.lattice.angle(self.Vb,self.Vc)
        bet = self.lattice.angle(self.Va,self.Vc)
        gam = self.lattice.angle(self.Va,self.Vb)
        return Lattice(a,b,c,alp,bet,gam)
    
##########################################################################
##########################################################################
def test_lattice():
    # create a new lattice instance
    cell = Lattice(5.6,5.6,13,90,90,120)

    # show lattice parameters and real
    # recip lattice cell volumes
    print cell
    print cell.vol()
    print cell.vol(recip=True)

    # what the hkl = [0,0,1] dspace
    # and 2 theta value for lambda = 1 ang.
    print cell.d([0,0,1])
    print cell.tth([0,0,1],lam=1.)
    
    # create a real space vector perpendicular
    # to hkl = [0,0,1], with magnitude = dspace(001)
    dv = cell.dvec([0,0,1])
    print dv
    print cell.mag(dv)

    # compute the angle between some vectors
    # recip lattice [1,0,0] and [0,1,0]
    print cell.angle([1,0,0],[0,1,0],recip=True)
    # recip lattice [1,0,0] and [0,0,1]
    print cell.angle([1,0,0],[0,0,1],recip=True)
    # recip lattice [0,0,1] and [1,1,1]
    print cell.angle([0,0,1],[1,1,1],recip=True)

    # real lattice [1,0,0] and [0,1,0]
    print cell.angle([1,0,0],[0,1,0])
    # real lattice [1,0,0] and [0,0,1]
    print cell.angle([1,0,0],[0,0,1])
    # real lattice [0,0,1] and [1,1,1]
    print cell.angle([0,0,1],[1,1,1])

    # real [0,0,1] and recip [0,0,1]
    print cell.angle_rr([0,0,1],[0,0,1])
    # real [1,1,1] and recip [0,0,1]
    print cell.angle_rr([1,1,1],[0,0,1])
    # real dv and recip [0,0,1]
    print cell.angle_rr(dv,[0,0,1])
    # real dv and recip [1,1,0]
    print cell.angle_rr(dv,[1,1,0])

    return cell

def test_transform():
    # create a new hexagonal lattice instance
    cell = Lattice(5.6,5.6,13,90,90,120)
    
    # test lattice transform
    # below defines rhombohedral basis vectors
    # in terms of the heaxagonal basis vectors
    Va_rhom = [0.6667, 0.3333, 0.3333]
    Vb_rhom = [-0.3333, 0.3333, 0.3333]
    Vc_rhom = [-0.3333, -0.6667, 0.3333]
    t = LatticeTransform(cell,Va=Va_rhom,Vb=Vb_rhom,Vc=Vc_rhom)

    # create an instance of the rhombohedral lattice
    # shows the rhomb lattice params
    print t.plat()

    # given the hkl = [001] in hex calc hkl in rhom
    print t.hp([0,0,1])
    
    # given the hkl = [111] in rhom calc hkl in hex
    print t.h([1,1,1])

    # create a cartesian representation
    t.cartesian()
    cart = t.plat()
    print cart
    
    # convert a [1,1,0] vector in hex lattice to cartesian
    vc = t.xp([1,1,0])
    print vc

    # check that the vector is the same length
    print cell.mag([1,1,0])
    print cart.mag(vc)

    return t


##########################################################################
##########################################################################
if __name__ == "__main__":
    """
    test 
    """
    test_lattice()
    test_transform()
    