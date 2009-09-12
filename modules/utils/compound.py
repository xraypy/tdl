########################################################################
"""
T. Trainor (fftpt@uaf.edu)
Utilities for handling the stiochiometry of chemical compounds

Modifications:
--------------


"""
########################################################################
"""
To Do
------
Note we should use the appropritate terminology, ie molar mass
rather than formula or molecular weight.  or use wt instead of mw

Note see the IUCR description of formulas in cif files....
http://ww1.iucr.org/cif/cifdic_html/1/cif_core.dic/Cchemical_formula.html

Use element class in component

Add charge stiochiometry

Add more units/conversions.  E.g. allow assigning stioch coeff in terms
of wt% or ppm and reporting components in mass and vol based units.

--------------------------------------------------------------------
Add methods __add__, __sub__, __mul__ etc... so we can do stiochiometric
arithmetic.  ie would be cool if:

Al = element('Al',oxid=3,isotope=None)
Fe = element('Fe',oxid=3)
O = element('O')

alumina = 2*Al + 3*O
hem = 2*Fe + 3*O
mix = 0.98*alumina + 0.02*hem

- notice that these all appear to be formation rxns...

- all the above expr should evaluate to a species class

- but we should have a method
     species.as_component() 
  which returns a component object having the net formula etc...
- then add a method to species to change the basis set:
     species.change_basis([comps])
  ie this recomputes stioch coefficients using new set of comps

"""
#########################################################################

import numpy as num
import types, exceptions

import element_data
import physcon as pc

#########################################################################
# element = elem/El/Z/z, conc = CZ, moles = nZ, mass = mZ, mole frac = fZ
class Element:
    def __init__(self, symbol, name, Z, wgt):
        self.sym  = symbol # short symbol
        self.name = name   # long name
        self.Z    = Z      # atomic number
        self.wgt  = wgt    # amu g/mol
        self.ox   = 0.     # oxidation state

#########################################################################
# component = comp/X, conc = CX, moles = nX, mass = mX, mole frac = fX
class Component:
    """
    * component = formula of elements
    
    * __getitem__ and __setitem__ operate on element
      stiochiometric coefficients.

    * note we arent currently using the element class, but should add...
    
    """
    #####################################################################
    def __init__(self,**kw):
        self.init(**kw)

    def init(self,formula=None):
        self.name      = ''  # computed string name 
        self.formula   = {}  # dict of {'El':nu}
        self.mw        = 0.0 # computed wgt
        if formula: self.set(formula)

    #####################################################################
    def __repr__(self,):
        self.update()
        lout = self.name + ': MW = ' + str(self.mw) + ' g/mol'
        return lout

    ####################################################################
    def __setitem__(self,el,nu):
        """
        set/add stioch of element
        el = string symbol or z
        note nu = moles_of_el / mole_of_component
        """
        if type(el) != types.StringType:
            el = element_data.symbol(el)
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
        return stioch of element
        el = string symbol or z 
        note nu = moles_of_el / mole_of_component
        """
        if type(el) != types.StringType:
            el = element_data.symbol(el)
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
            aw = element_data.amu(el)
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
        formula is a dictionary: {'El':nu}
        where 'El' is the atomic symbol
        and nu is the stoichiometric coefficient

        Note nu = moles_of_El / mole_of_component
        therefore nu may be fractional and are converted to
        float even if passed as int
        
        This methods appends new formula units to existing
        """
        if type(formula) == types.StringType:
            formula = {formula:1.}
            
        for el,nu in formula.items():
            if type(el) != types.StringType:
                el = element_data.symbol(el)
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
        of the element given by Z (int or str)
        """
        if type(el) != types.StringType:
            el = element_data.symbol(el)
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
# species = spec/S/s, conc = CS, moles = nS, mass = mS, mole frac = fS
class Species:
    """
    * species = chemical compound made up of one or more components
    * __getitem__ and __setitem__ operate on component
      stiochiometric coefficients.  
    """
    ####################################################################
    def __init__(self,**kw):
        self.init(**kw)

    def init(self,comp=[]):
        self.name   = ''   # computed name string
        self.comp   = []   # list of (comp,nu)
        self.mw     = 0.0  # computed wieght
        if len(comp) > 0:
            self.set(comp)

    ####################################################################
    def __repr__(self,):
        self.update()
        lout = self.name + ': MW = ' + str(self.mw) + 'g/mol'
        return lout

    ####################################################################
    def __setitem__(self,comp,nu):
        """
        set the component stoichiometric coefficent
        comp name or int idx
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
            self.comp[idx] = (f,nu)
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
        comp name or int idx
        """
        idx = -1
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
        compute molecular wieght
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
        self.comp = []
        self.mw   = 0.0
        self.name = ''

    ####################################################################
    def set(self,comp):
        """
        comp is a list of [(formula,nu)]
        where formula is a formula instance and
        nu is the stoichiometric coefficient
        note nu may be fractional (and are converted to
        float even if passed as int)
        
        This appends new comp units to existing

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
            el = element_data.symbol(el)
            if el == None: return 0.0
        el = el.title()
        for (comp,nu) in self.comp:
            if el in comp.formula.keys():
                nu_tot = nu_tot + nu*comp.formula[el]
        return nu_tot

    ####################################################################
    def elements(self):
        elem = []
        for (comp,nu) in self.comp:
            for el in comp.formula.keys():
                if el not in elem:
                    elem.append(el)
        return elem
    
    ####################################################################
    def components(self):
        names = []
        for (comp,nu) in self.comp:
            name.append(comp.name)
        return names

########################################################################
# material = mat/M/m
class Material(Species):
    """
    same as species class except
    - we assume its a homogeneous material of fixed density.
    - given a (unit) volume we may also compute total moles etc.
    """
    def __init__(self,**kw):
        self.init_mat(**kw)
        
    def init_mat(self,comp=[],density=1.0,volume=1.0):
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
###################################################################
def parse_species_formula(expr):
    """
    convert string expr to a component or species object, e.g.
        {(2)[A1_B2]}
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
    convert string (e.g. Si_O2)
    to a component object
    """
    if len(compstr) == 0: return None
    items = compstr.split("_")

    formula = {}
    for item in items:
        if len(item) == 1:  
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
    return Component(formula=formula)

def ideal_gas(p=1.,t=298.15,mw=0.0):
    """
    compute the gas phase concentration - mole/cm^3
    and density
    p in atm
    t in K
    mw in g/mole
    """
    m2cm = 10.**-6 #m^-3 to cm^-3
    a2p = 101.325*(10**3) #pa/atm
    p = p*a2p
    c = p/(pc.R*t) # conc in moles/m^3
    c = c*m2cm     # conc in moles/cm^3
    rho = c*mw     # g/cm3
                   # x 1e3 = kg/m^3
    return (c,rho)

#############################################################
#############################################################
def test_1():
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
    expr = "{(1.0)[Si_O2]}{(1.232)[Fe2.3454_O11]}{(.12)[H_C.2_Al0.0001203]}"
    spec = parse_species_formula(expr)
    print spec

#############################################################
if __name__ == "__main__":
    #test_1()
    test_2()
    