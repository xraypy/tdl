#   This module provides a number of utility functions for
#   X-ray fluorescence (XRF)
#
#   Author:  Mark Rivers
#   Created: Sept. 16, 2002
#   Modifications:

nels = 100   # Number of elements in table
nlines = 14  # Number of x-ray lines in table

atomic_symbols = (
   'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', 'Na', 'Mg', 'Al',
   'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe',
   'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 
   'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
   'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
   'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt',
   'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
   'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm')

# The xrf_lines database is a dictionary in the form {'FE': {'KA':6.400}}
from xrf_lines import xrf_lines

# lines for common radioactive sources
gamma_lines = {'CO57':{'G1':14.413,'G2':122.0614,'G3':136.4743},
              'CD109':{'G1':88.04}}

def atomic_number(symbol):
      """
      Returns the atomic number of an element, given its atomic symbol 
  
      Inputs:
         symbol: The atomic symbol of the element whose atomic number is being
                 requested.  This is a 1 or 2 character string, e.g. 'H', 'Si',
                 etc.  It is case insensitive and leading or trailing blanks
                 are ignored.

      Outputs:
         This function returns the atomic number of the input element.  If an
         invalid atomic symbol is input then the function returns 0.
      """
      s = symbol.split()[0].upper()
      for i in range(len(atomic_symbols)):
         if (s == atomic_symbols[i].upper()): return i+1
      return None


def atomic_symbol(z):
      """
      This function returns the atomic symbol of an element, given its atomic
      number.
      
      Inputs:
         z: The atomic number of the element whose atomic symbol is being
            requested.  

       Outputs:
          This function returns the atomic symbol of the input element as a
          string.  If Z is an invalid atomic number then the function returns 
          None.
      """
      if (z > 0) and (z <= len(atomic_symbols)): return(atomic_symbols[z-1])
      return None


def lookup_xrf_line(xrf_line):
   """
   This function returns the energy in keV for a particular x-ray
   fluorescence line.

   Inputs:
      xrf_line: A string of the form 'Element Line', where Element is an
                atomic symbol, and Line is an acronym for a fluorescence line.
                Both Element and Line are case insensitive.  There must be a space
                between Element and Line.
                The valid lines are
                 ka  - K-alpha (weighted average of ka1 and ka2)
                 ka1 - K-alpha 1
                 ka2 - K-alpha 2
                 kb  - K-beta (weighted average of kb1 and kb2)
                 kb1 - K-beta 1
                 kb2 - K-beta 2
                 la1 - L-alpha 1
                 lb1 - L-beta 1
                 lb2 - L-beta 2
                 lg1 - L-gamma 1
                 lg2 - L-gamma 2
                 lg3 - L-gamma 3
                 lg4 - L-gamma 4
                 ll  - L-eta

      Examples of XRF_Line:
          'Fe Ka'  - Fe k-alpha
          'sr kb2' - Sr K-beta 2
          'pb lg2' - Pb L-gamma 2

   Outputs:
      This function returns the fluoresence energy of the specified line.
      If the input is invalid, e.g. non-existent element or line, then the
      function returns None

   Examples:
      energy = lookup_xrf_line('Fe Ka')  ; Look up iron k-alpha line
      energy = lookup_xrf_line('Pb lb1') ; Look up lead l-beta 1 line
   """

   try:
      words = xrf_line.split()
      element = words[0].strip().upper()
      line    = words[1].strip().upper()
      return xrf_lines[element][line]
   except:
      return None


def lookup_gamma_line(gamma_line):
   """
   Returns the energy in keV for a particular gamma emmission line.

   Inputs:
      gamma_Line: A string of the form 'Isotope Line', where Isotope is a
                  the symbol for a radioactive isotope, and Line is an index of the form
                  g1, g2, ... gn.
                  Both Isotope and Line are case insensitive.  There must be a space
                  between Isotope and Line.

         Examples of Gamma_Line:
            'Cd109 g1'  - 88 keV line of Cd-109
            'co57 g2'   - 122 keV line of Co-57

   Outputs:
      This function returns the gamma energy of the specified line.
      If the input is invalid, e.g. non-existent isotope or line, then the
      function returns 0.

   Restrictions:
       This function only knows about a few isotopes at present.  It is
       intended for use with common radioactive check sources.  It is easy to
       add additional isotopes and gamma lines to this function as needed.
       The current library is:
           'CO57 G1' = 14.413
           'CO57 G2' = 122.0614
           'CO57 G3' = 136.4743
           'CD109 G1'= 88.04

   Example:
       energy = lookup_gamma_line('Co57 g1')  ; Look up 14 keV line of Co-57
   """
   try:
      words = gamma_line.split()
      isotope = words[0].strip().upper()
      line    = words[1].strip().upper()
      return gamma_lines[isotope][line]
   except:
      return None
