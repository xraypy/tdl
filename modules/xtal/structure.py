"""
Crystal structure generator

Authors/Modifications:
----------------------
* Tom Trainor (tptrainor@alaska.edu) 

Todo:
-----
* UnitCell class
* Reading cif files and others (maybe seperate module)
   - e.g. get xyz file from fit and par files
* Structure analysis and bond valence calcs (seperate module...)
* Include dictionary of space groups (a seperate module)
* Handle 2D plane group operations (include dictionary of plane groups) 
 
"""
##########################################################################

import numpy as num
import sys
import sets
import types
import string
import lattice

# try importing the aussie pycif module
# if it wont load we should have our own simple backup
# for reading/writing structures as cif files
try:
    
    import CifFile
    import StarFile
except:
    print "No Pycif module available"

##########################################################################
class UnitCell:
    """
    Need to include assymetric unit and thermal factors
    Also include slab repeat vector as an option (or maybe we
    should make a new SurfaceCell class that subclasses this?)
    """
    def __init__(self):
        # assymetric unit, debye wallers, lattice (Lattice instance)
        # and sym_ops (PositionGenerator instance)
        self.assym = []
        self.dw = []
        self.lat = None
        self.PositionGenerator = []
        pass
    
    def read_cif(self,fname):
        cf = CifFile.ReadCif(fname)
        #Converts CifFile to CifBlock which is readable
        cb = cf['global']
        positions = cb.GetLoop('_atom_site_label')
        positions = positions.values()
        atoms = positions[0]
        num_atoms = len(atoms)
        x = positions[1]
        y = positions[2]
        z = positions[3]
        counter = 0
        n = 0
        assym_unit = []
        elem_positions = {}
        
        while (counter < num_atoms): 
            elem_positions = {'name':atoms[n], 'x':x[n], 'y':y[n], 'z':z[n]}
            assym_unit.append(elem_positions)
            counter = counter + 1
            n = n + 1

        self.assym_unit = assym_unit

        a = cb['_cell_length_a']
        self.a = a
        b = cb['_cell_length_b']
        self.b = b
        c = cb['_cell_length_c']
        self.c = c

        alpha = cb['_cell_angle_alpha']
        self.alpha = alpha
        beta = cb['_cell_angle_beta']
        self.beta = beta
        gamma = cb['_cell_angle_gamma']
        self.gamma = gamma

         # store values
        self.lat = lattice.Lattice(a,b,c,alpha,beta,gamma)
        print ""
        print ""
        print "Lattice Information is as follows:"
        print self.lat
        
        sym_operations = cb.GetLoop('_space_group_symop_operation_xyz')
        sym_values = sym_operations.values()
        # This converts all the operations into a list of strings
        sym_values = sym_values[0] 

        num_opers = len(sym_values)
        p = PositionGenerator()

        dict_length = len(assym_unit)
        counter = 0
        position = 0
        atom_name = 1
        
        while (counter < dict_length):
            print "Fractional coordinates for ", atom_name
            while (position < num_opers):
                split = sym_values[position]
                xyz = string.split(split, ',')
                x=xyz[0]
                y=xyz[1]
                z=xyz[2]
                matrix = p._make_seitz_matrix(x,y,z)
                p.ops.append(matrix)
                position = position + 1
            vecs = p.copy(0.11,0.5,0.5,reduce=True,rem_dups=True)
            for v in vecs: print v   
            counter = counter + 1
            position = 0
            atom_name = atom_name + 1

       #Store positions
        self.symm = PositionGenerator()
        # loop through and add sym ops to self.symm
        
    def write_cif(self):
        pass
   
    def generate_p1(self,na=1,nb=1,nc=1):
        cell = self.lattice(a,b,c,alpha,beta,gamma) 
        # cell is a new lattice instance

        # Now generate positions

        
            
    def transform(self):
        """
        this computes/returns a new UnitCell
        given basis transform vectors.
        Note this needs to generate a new Lattice,
        new coordinates for assymetric unit,
        correclty transform thermal tensor,
        and transform the symmetry operators.
        This will be a fun code to work on!!!
        """
        pass
    def bond_valence(self):
        """
        compute bond valence sums (and coordination chem ie
        coordation sphere, bond lenghts and angles)
        """
        pass
    def visualize(self,na=1,nb=1,nc=1):
        """
        output a jmol script to view the structure
        add more arguments for controlling jmol attributes 
        """
        pass
    
##########################################################################
class PositionGenerator:
    """
    Class to generate equivalent positions given symmetry operations
    """
    ###########################################################    
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

def test2():
    uc = UnictCell()
    uc.read_cif('pyrite.cif')
    p1_coords = uc.generate_p1()
    print p1_coords

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
