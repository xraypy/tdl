"""
Crystal structure generator

Authors/Modifications:
----------------------
* Tom Trainor (tptrainor@alaska.edu) 

Todo:
-----
* Include dictionary of space groups from the international tables
* Handle 2D plane group operations (include dictionary of plane groups) 
* Reading cif files and others (maybe seperate module)
   - e.g. get xyz file from fit and par files
* Structure analysis and bond valence calcs (seperate module...)
 
"""
##########################################################################

import numpy as num
import sys
import types

##########################################################################
class PositionGenerator:
    """
    Class to generate equivalent positions given symmetry operations
    """
    def __init__(self):
    """ init """
        self.ops = []

    ###########################################################    
    def add_op(self,sym='x,y,z',shift=''):
        """
        Add a new symmetry operator for generating positions

        Parameters:
        ----------
        * sym and shift are strings with comma delimeted set of
          characters that defines the opertions.  e.g.
          sym = "x,y,z", shift = "0,0,0"
          sym = "-y,z,x+y", shift = "1/2,1/2,0"

        Example:
        --------
        >>sym1 = "x,y,z"
        >>shift1 = "1/2, 1/2, 0"
        >>p.add_op(sym=sym1,shift=shift1)
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
        Generate augmented (seitz) matrix given
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
        Calc all sym copies of a position

        Parameters:
        -----------
        * x,y,z are fractional coordinates
        * reduce is flag to indicate that all
          positions must be in bounds  0 to 1
        * rem_dups is a flag to indicate if duplicates
          should be removed

        Outputs:
        --------
        * list of vectors of symmetry copy positions
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
    sym1   = "x,y,z"
    sym2   = "-x,y,-z"
    sym3   = "-x,-y,-z"
    sym4   = "x,-y,z"
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
