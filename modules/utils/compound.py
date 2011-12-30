"""
Utilities for handling the stiochiometry of chemical compounds

Authors/Modifications:
----------------------
* T. Trainor (tptrainor@alaska.edu)

Notes:
------
This module defines a set of routines for working with
chemical compounds.  The following are defined:
* element = a pure chemical element
* component = a chemical entity described in terms of element stiochiometry
* species = a chemical entity made up from one or more chemical components in
            definite proportions.  Species has a well define composition
            and structure
            
* gas
* liquid
* solution 
* material = a (solid) species with a fixed density (we assume homogeneous) 
   - amorphous
   - xtaline
* mixture

The following shorthand/nomenclature is used in this module:
* element = elem/El/Z/z, conc = CZ, moles = nZ, mass = mZ, mole frac = fZ
* component = comp/X, conc = CX, moles = nX, mass = mX, mole frac = fX
* species = spec/S/s, conc = CS, moles = nS, mass = mS, mole frac = fS
* material = mat/M/m
* stiochiometric coefficient = nu


Examples:
---------
# create some chemical components, these can be
# made up from an arbirtrary stiochiometry
>>c1 = Component(formula={'Si':1.,'O':2})

# add some more elements/change element stiochiometry
>>c1.set({'O':4,'H':2})
>>print c1

# print the stiochiometric coefficient of O
>>print c1.nuZ('O')

# make another component
>>c2 = Component()
>>c2.set({"H":2.,"O":1})

# make a species by combining the components in a given proportion
>>spec1 = Species(comp=[(c1,2),(c2,4)])
>>print spec1
>>print spec1.elements()
>>print spec1[1]
>>print spec1.nuZ('O'), spec1.nuZ('Si'), spec1.nuZ('H')

# this computes a species from an expression.  
>>expr = "{(1.0)[Si_O2]}{(1.232)[Fe2.3454_O11]}{(.12)[H_C.2_Al0.0001203]}"
>>comps = parse_fmt_formula(expr)
>>spec = Species(comp=comps)

To Do
------
* We should use the appropriate terminology, ie molar mass
  rather than formula or molecular weight.  or use wt instead of mw

* Note see the IUCR description of formulas in cif files....
  http://ww1.iucr.org/cif/cifdic_html/1/cif_core.dic/Cchemical_formula.html

* Use element class for elements...  and define a ChemicalFormula base class
  that holds the basic info related to chemical formula's.  a Component
  is therefor just a formula class, and a species is a collection
  of formula classes.  Should be able to maka a simple formula instance using:
  >>f = formula('Al2 O3')

* Add charge stiochiometry

* Add more units/conversions.  E.g. allow assigning stioch coeff in terms
  of wt% or ppm and reporting components in mass and vol based units.

* Add methods __add__, __sub__, __mul__ etc... so we can do stiochiometric
  arithmetic.  ie would be cool if:

  Al = element('Al',oxid=3,isotope=None)
  Fe = element('Fe',oxid=3)
  O  = element('O')

  alumina = 2*Al + 3*O
  hem = 2*Fe + 3*O
  mix = 0.98*alumina + 0.02*hem

  all the above expr should evaluate to a species class

* Should have a method
     species.as_component() 
  which returns a component object having the net formula etc...

* then add a method to species to change the basis set:
     species.change_basis([comps])
  ie this recomputes stioch coefficients using new set of comps

"""
#########################################################################

import numpy as num
import types, exceptions

#import element_data
#import physcon as pc
from tdl.modules.utils import elements
from tdl.modules.utils import physcon as pc

#########################################################################
class Element:
    """
    An element
    """
    def __init__(self, symbol, name, Z, wgt):
        self.sym  = symbol # short symbol
        self.name = name   # long name
        self.Z    = Z      # atomic number
        self.wgt  = wgt    # amu g/mol
        self.ox   = 0.     # oxidation state

#########################################################################
class Component:
    """
    A component is a entity built in terms of a formula of elements

    Notes:
    ------
    * __getitem__ and __setitem__ operate on element
      stiochiometric coefficients.  The stiochiometric
      coefficients (nu) are moles_of_element / mole_of_component.

    Examples:
    ---------
    # create a component
    >>c = Component(formula={'H':2,'O':1})
    or
    >>c = Component(formula='H2 O')
    >>c['H'] = 10
    >>x = c['H']
    # we can also modify it by reference to element Z value
    >>c[8]=12  

    Todo:
    -----
    * we arent currently using the element class, but should add...
    
    """
    #####################################################################
    def __init__(self,**kw):
        """
        Initialize.  See init
        """
        self.init(**kw)

    def init(self,formula=None):
        """
        (re)init the object
        
        Parameters:
        -----------
        * formula is a dictionary of {'El',nu} or {Z,nu}
          or a string that can be parsed by the parse_formula
          function
        """
        self.name      = ''  # computed string name 
        self.formula   = {}  # dict of {'El':nu}
        self.mw        = 0.0 # computed wgt
        if formula: self.set(formula)

    #####################################################################
    def __repr__(self,):
        """
        """
        self.update()
        lout = self.name + ': MW = ' + str(self.mw) + ' g/mol'
        return lout

    ####################################################################
    def __setitem__(self,el,nu):
        """
        set/add stioch of element
        """
        if type(el) != types.StringType:
            el = elements.symbol(el)
            if el == None: return
        el = el.title()
        if el in self.formula.keys():
            self.formula[el] = float(nu)
        else:
            self.formula.update({el:float(nu)})
        self.update()
        
    ####################################################################
    def __getitem__(self,el):
        """
        return stiochiometric coefficient of element
        """
        if type(el) != types.StringType:
            el = elements.symbol(el)
            if el == None: 0.0
        el = el.title()
        if el in self.formula.keys():
            return self.formula[el]
        else:
            return 0.0
        
    ####################################################################
    def _calc_name(self,):
        """
        compute a string repr
        """
        name=''
        for key in self.formula.keys():
            name="%s_%s%g"% (name,key,self.formula[key])
        name = name[1:]
        self.name = name
    
    ####################################################################
    def _calc_mw(self,):
        """
        compute molecular wieght (g/mole)
        """
        mw = 0.0
        for el in self.formula.keys():
            nu = self.formula[el]
            aw = elements.amu(el)
            mw = mw + nu*aw
        self.mw = mw
        
    ####################################################################
    def update(self,):
        """
        recompute molec weight and name
        """
        self._calc_name()
        self._calc_mw()

    ####################################################################
    def clear(self,):
        """
        clear formula
        """
        self.formula = {}
        self.mw      = 0.0
        self.name    = ''

    ####################################################################
    def set(self,formula):
        """
        This methods appends to or modifies the existing formula

        Paramters:
        ----------
        * formula is a dictionary: {'El':nu}
          where 'El' is the atomic symbol or Z (integer)
          and nu is the stoichiometric coefficient

          Or formula is a string: 'Si O2' (see parse_formula)

          Note nu = moles_of_El / mole_of_component
          therefore nu may be fractional and are converted to
          float even if passed as int
        
        """
        #if type(formula) == types.StringType:
        #    formula = {formula:1.}
        if type(formula) == types.StringType:
            formula = parse_formula(formula)
            if formula == None: return
            
        for el,nu in formula.items():
            if type(el) != types.StringType:
                el = elements.symbol(el)
                if el == None:
                    print "Error setting component formula"
                    return
            el = el.title()
            if el in self.formula.keys():
                self.formula[el] = float(nu)
            else:
                self.formula.update({el:float(nu)})
        self.update()
        
    ####################################################################
    def nuZ(self,el):
        """
        return the stoichiometric coefficient
        
        Paramters:
        ----------
        * el is a string 'El' or integer Z
        """
        if type(el) != types.StringType:
            el = elements.symbol(el)
            if el == None: return 0.0
        el = el.title()
        if el in self.formula.keys():
            return self.formula[el]
        else:
            return 0.0

    ####################################################################
    def elements(self,):
        """
        return list of elements
        """
        return self.formula.keys()

########################################################################
class Species:
    """
    A species is a chemical compound made up of one or more components

    Notes:
    ------
    * __getitem__ and __setitem__ operate on component
      stiochiometric coefficients.  The stiochiometric
      coefficients (nu) are moles_of_component / mole_of_species.

    Examples:
    ---------
    # create components
    >>c1 = Component(formula={'H':2,'O':1})
    >>c2 = Component(formula={'Al':2,'O':3})

    # combine components,
    # e.g. 4 moles of c1 and 1 mole of c2 per mole of s
    >>s = Species(comp=[(c1,4),(c2,1)])

    # change the stiochiometry of the 2nd component
    >>s[1] = 8.2

    # You can also pass the formula directly.
    # this just creates each element as a component.  
    >>s = Species(formula={'H':2,'O':1.})
    or
    >>s = Species(formula='{(1)[H2 O]}')
    """
    ####################################################################
    def __init__(self,**kw):
        """
        Initialize.  See init
        """
        self.init(**kw)

    def init(self,comp=[],formula=None):
        """
        (re)init the object
        
        Parameters:
        -----------
        * comp is a list of components [(c,nu),(c,nu)]
        * formula is a dictionary of {'Element':nu} or
          a string that can be parsed by the parse_fmt_formula
          function

        Note use one or the other of the above methods to init
        in the case you use 'formula' then all elements are taken
        as components.  
        """
        self.name   = ''   # computed name string
        self.comp   = []   # list of (comp,nu)
        self.mw     = 0.0  # computed wieght
        if len(comp) > 0:
            self.set(comp)
        elif formula != None:
            if type(formula) == types.StringType:
                comp = parse_fmt_formula(formula)
                if comp == None:
                    print "Error parsing formula string"
                    return
            elif type(formula) == types.DictType:
                comp = []
                for el in formula.keys():
                    c = Component(formula={el:1.})
                    comp.append((c,formula[el]))
            #print comp
            self.set(comp)

    ####################################################################
    def __repr__(self,):
        """
        """
        self.update()
        lout = self.name + ': MW = ' + str(self.mw) + 'g/mol'
        return lout

    ####################################################################
    def __setitem__(self,comp,nu):
        """
        set the component stoichiometric coefficent
        """
        idx = -1
        nu = float(nu)
        if type(comp) != types.StringType:
            idx = int(comp)
            if comp < 0 or comp > len(self.comp):
                raise exceptions.IndexError
        else:
            for j in range(len(self.comp)):
                c = self.comp[j][0]
                if c.name == comp:
                    idx = j
                    break
        if idx >= 0:
            c  = self.comp[idx][0]
            self.comp[idx] = (c,nu)
        else:
            # this should convert string to a comp object
            #comp = parse_comp(comp)
            #self.comp.append[(comp,nu)]
            pass
        self.update()
        
    ####################################################################
    def __getitem__(self,comp):
        """
        return the stoich coeff of the component

        Paramters:
        ---------
        * comp - component class instance, or string or the component index
        """
        idx = -1
        if hasattr(comp,'name'):
            comp = comp.name
        if type(comp) != types.StringType:
            idx = int(comp)
            if comp < 0 or comp > len(self.comp):
                raise exceptions.IndexError
        else:
            for j in range(len(self.comp)):
                c = self.comp[j][0]
                if c.name == comp:
                    idx = j
                    break
        #print idx
        if idx in range(len(self.comp)):
            return self.comp[idx][1]
        else:
            return 0.0

    ####################################################################
    def _calc_name(self,):
        """
        compute name
        """
        name=''
        for (comp,nu) in self.comp:
            name="%s{(%s)[%s]}" % (name,str(nu),comp.name)
        self.name = name
    
    ####################################################################
    def _calc_mw(self,):
        """
        compute molecular wieght (g/mol)
        """
        mw = 0.0
        for (comp,nu) in self.comp:
            mw = mw + nu*comp.mw
        self.mw = mw
        
    ####################################################################
    def update(self,):
        """
        update 
        """
        for (comp,nu) in self.comp:
            comp.update()
        self._calc_name()
        self._calc_mw()

    ####################################################################
    def clear(self,):
        """
        clear formula
        """
        self.comp = []
        self.mw   = 0.0
        self.name = ''

    ####################################################################
    def set(self,comp):
        """
        This appends new component units to existing species

        Parameters:
        -----------
        * comp is a list of [(c,nu)]
          where c is a Component instance and
          nu is the stoichiometric coefficient
          note nu may be fractional (and are converted to
          float even if passed as int)

          Note these are passed by reference, so if you
          change the components it will change the compound!
        """
        if type(comp) != types.ListType:
            comp = [comp]
        for xx in comp:
            if len(xx) == 1:
                c = xx
                nu = 1.0
            elif len(xx) == 2:
                c  = xx[0]
                nu = float(xx[1])
            #if type(c) == types.StringType:
            #   c = parse_comp(comp)
            self.comp.append((c,nu))
        self.update()
        
    ####################################################################
    def nuX(self,comp):
        """
        return the stoichiometric coefficient
        of the component (int or str)
        """
        return self.__getitem__(comp)

    ####################################################################
    def nuZ(self,el):
        """
        return the stoichiometric coefficient
        of the element given by Z (int or str)
        """
        nu_tot = 0.0
        if type(el) != types.StringType:
            el = elements.symbol(el)
            if el == None: return 0.0
        el = el.title()
        for (comp,nu) in self.comp:
            if el in comp.formula.keys():
                nu_tot = nu_tot + nu*comp.formula[el]
        return nu_tot

    ####################################################################
    def elements(self):
        """
        return list of elements
        """
        elem = []
        for (comp,nu) in self.comp:
            for el in comp.formula.keys():
                if el not in elem:
                    elem.append(el)
        return elem

    ####################################################################
    def net_formula(self):
        """
        Return the net elemental formula
        """
        elem = self.elements()
        formula = {}
        for c,nu in self.comp:
            nu = float(nu)
            f = c.formula
            for el,enu in f.items():
                if formula.has_key(el):
                    formula[el] = formula[el] + float(enu)*nu
                else:
                    formula.update({el:float(enu)*nu})
        return formula
    
    ####################################################################
    def components(self):
        """
        return list of components
        """
        names = []
        for (comp,nu) in self.comp:
            names.append(comp.name)
        return names

########################################################################
class Material(Species):
    """
    A material is the same as species class except
    - we assume its homogeneous with a fixed (and specified) density.
    - given a (unit) volume we may also compute total moles etc.
    
    """
    def __init__(self,**kw):
        """
        Initialize.  See init
        """
        self.init_mat(**kw)
        
    def init_mat(self,comp=[],density=1.0,volume=1.0):
        """
        (re)init the object
        
        Parameters:
        -----------
        * comp is a list of components [(c,nu),(c,nu)]
        * density is the density in g/cm^3
        * volume is the total volume of material in cm^3
        """
        self.density     = density     # density, g/cm^3
        self.vol         = volume      # cm^3
        self.init(comp=comp)

    #######################################################################
    def concZ(self,z):
        """
        compute element concetration moles/cm^3
        """
        self.update()
        nu    = self.nuZ(z)
        conc = (self.density/self.mw)*nu
        return conc

    #######################################################################
    def concX(self,comp):
        """
        compute component concetration moles/cm^3
        """
        self.update()
        nu    = self.nuX(comp)
        conc  = (self.density/self.mw)*nu
        return conc

    #######################################################################
    def molesZ(self,z):
        """
        compute total moles of given element in the volume
        """
        self.update()
        moles = self.concZ(z)*self.vol
        return moles

    #######################################################################
    def molesX(self,comp):
        """
        compute total moles of given component in the volume
        """
        self.update()
        moles = self.concX(comp)*self.vol
        return moles
    
    #######################################################################
    def fracZ(self,z):
        """
        compute mole fraction of z
        """
        nZ = self.molesZ(z)
        nt = 0.0
        for z in self.elements():
            nT = nT + self.molesZ(z)
        return nZ/nT
    
    #######################################################################
    def fracX(self,comp):
        """
        compute mole fraction of comp
        """
        nX = self.molesX(z)
        nT = 0.0
        for c in self.components():
            nT = nT + self.molesC(z)
        return nX/nT

    #######################################################################
    def numZ(self,z):
        """
        compute total number of atoms of given element in the volume
        """
        self.update()
        moles = self.molesZ(z)
        return moles*pc.N_A

    #######################################################################
    def numX(self,comp):
        """
        compute total number of component units in the volume
        """
        self.update()
        moles = self.molesX(comp)
        return moles*pc.N_A    

###################################################################
def parse_fmt_formula(expr):
    """
    Convert string expr to a list of component objects

    The return value is of the form: [(comp,nu),(comp,nu)]

    Parameters:
    -----------
    * expr is a string expression of the form:
         '{(2)[A1_B2]}{(1.1)[D_F0.1]}...'
       where
       - {} deliniates components.
       - inside the curly brackets the (2) specifies the
         number of [A1_B2] (=component) formula units

    Todo:
    -----
    put this and similiar in chem_parse module, and improve
    """
    expr = expr.strip()
    if len(expr) == 0: return
    x = expr
    formula = []
    while 1:
        st = x.find("{")
        if st == -1: break
        en = x.find("}")
        compstr = x[st+1:en]
        #
        x = x[en+1:].strip()
        #
        st = compstr.find("(")
        en = compstr.find(")")
        numstr = compstr[st+1:en]
        st = compstr.find("[")
        en = compstr.find("]")
        compstr = compstr[st+1:en]
        #
        nu   = float(numstr)
        comp = parse_comp(compstr)
        formula.append((comp,nu))
        #
        if len(x) == 0: break
    #
    return formula

def parse_comp(compstr):
    """
    Convert a string to a component object

    Parameters:
    ----------
    * compstr is a string with a formula unit representing a chemical component
      the string should be formatted as e.g. 'Si O2' or 'Si_O2'.  ie the space
      or underscore is needed to seperate elements in the formula.  The element
      stoichiometry needs to be adjacent to the element (no space).  If there is
      no coefficient it is assumed to be 1.  
    """
    formula = parse_formula(compstr)

    return Component(formula=formula)

def parse_formula(fstr):
    """
    Convert a string to a formula list 

    Parameters:
    ----------
    * fstr is a chemical formula string 
      The string should be formatted as e.g. 'Si O2' or 'Si_O2'.  ie the space
      or underscore is needed to seperate elements in the formula.  The element
      stoichiometry needs to be adjacent to the element (no space).  If there is
      no coefficient it is assumed to be 1.
    
    Output:
    --------
    * formula: a list of [(element,nu),(element,nu)...]
    
    """
    if len(fstr) == 0: return None
    fstr = fstr.strip()
    fstr = fstr.replace(' ','_')
    items = fstr.split("_")

    formula = {}
    for item in items:
        if len(item) == 0:
            pass
        elif len(item) == 1:  
            el = item[0].upper()
            nu = 1.0
        elif len(item) == 2:
            el = item[0].upper()
            if item[1].isalpha():
                el = el + item[1].lower()
                nu = 1.0
            else:
                nu = float(item[1])
        else:
            el = item[0].upper()
            item = item[1:]
            if item[0].isalpha():
                el = el + item[0].lower()
                item= item[1:]
            nu = float(item[:])
        #
        formula.update({el:nu})
    return formula

###################################################################
def ideal_gas(p=1.,t=298.15,mw=0.0):
    """
    compute the gas phase concentration - mole/cm^3 and density

    Parameters:
    -----------
    * p is pressure in atm
    * t is temperature in K
    * mw is molecular weight in g/mole
    """
    m2cm = 10.**-6        # m^-3 to cm^-3
    a2p = 101.325*(10**3) # pa/atm
    p = p*a2p
    c = p/(pc.R*t)        # conc in moles/m^3
    c = c*m2cm            # conc in moles/cm^3
    rho = c*mw            # g/cm3
                          # x 1e3 = kg/m^3
    return (c,rho)

###################################################################
def element_mass_fraction(f,el):
    """
    Compute the mass fraction of an element in a species
    (multiply by 100 for wt%, 1e6 for ppm and 1e9 for ppb)

    Parameters
    ----------
    * f is a species or component instance (or string - see parse_formula)
    * el is a a string, e.g. 'Al', 'As' etc..

    """
    if type(f) == types.StringType:
        f = parse_formula(f)
        f = Species(formula=f)
    nu = f.nuZ(el)
    if nu == 0: return 0.
    aw = elements.amu(el)
    frac = nu*aw/f.mw
    return frac

def mix_weight_percent(mixture,el):
    """
    Compute wt% of element in mixture

    Parameters:
    -----------
    * mixture is a list of [(f,wt%),(f,wt%)..]
      where f is a formula string (see parse_formula)
      or a species or component instance
    * el is an element string
    """
    frac = 0.
    for f,wtp in mixture:
        xx = wtp*element_mass_fraction(f,el)
        frac = frac + xx
    return frac

def ppm_to_molar(spec,ppm):
    """
    given a ppm concentration of a species
    compute the molar concentration (assuming aqueous
    system and 1kg/L water density)
    """
    ppm = float(ppm)
    C = ppm * (1.0e-3)
    C = C / spec.mw
    return C

def molar_to_ppm(spec,molar):
    """
    given a molar concentration of a species
    compute the ppm concentration (assuming aqueous
    system and 1kg/L water density)
    """
    molar = float(molar)
    ppm = molar * (1.0e3)
    ppm = ppm * spec.mw
    return ppm

#############################################################
#############################################################
def test_1():
    """ test """
    c1 = Component(formula={'Si':1.00001,'O':2})
    c1.set({'O':2,'H':2})
    print "****"
    print c1
    print c1.nuZ('O')
    print "****"
    c2 = Component()
    c2.set({"H":2.,"O":1})
    print c2
    print "****"
    spec1 = Species(comp=[(c1,2),(c2,4)])
    print spec1
    print spec1.elements()
    print spec1[1]
    print spec1.nuZ('O'), spec1.nuZ('Si'), spec1.nuZ('H')
    print "****"

def test_2():
    """ test """
    expr = "{(1.0)[Si_O2]}{(1.232)[Fe2.3454_O11]}{(.12)[H_C.2_Al0.0001203]}"
    #spec = parse_fmt_formula(expr)
    spec = Species(formula=expr)
    print spec.net_formula()
    return spec

def test_3():
    """ test """
    s = Species(formula={'Na':2.,'S':1.,'O':4.})
    print s
    C = ppm_to_molar(s,10000.)
    print C
    ppm = molar_to_ppm(s,C)
    print ppm

def test_4():
    """test"""
    fstr = "Si3.9 Al.86 Li.08 Fe.1 Mg.14 O10 H O2"
    x = parse_formula(fstr)
    sp = Species(formula=x) 
    print x,sp
    return x,sp

#############################################################
if __name__ == "__main__":
    spec=test_2()
    #x,sp = test_4()
    