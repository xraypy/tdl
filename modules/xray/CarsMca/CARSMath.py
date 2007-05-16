"""
This module contains miscellaneous math functions

Author:
   Mark Rivers

Created: 
   Sept. 16, 2002

Modifications:
   Mar 28 2007  MN: we're only using compress_array and expand_array --
                got rid of cruft and Numeric dependencies.

   Sept. 26, 2002  MLR
      - Added newton function from scipy.optimize.  Put here so users don't
        need scipy
"""

import numpy as N

############################################################
def compress_array(array, compress):
   """
   Compresses an 1-D array by the integer factor "compress".  
   Temporary fix until the equivalent of IDL's 'rebin' is found.
   """

   l = len(array)
   if ((l % compress) != 0):
      print 'Compression must be integer divisor of array length'
      return array

   temp = N.resize(array, (l/compress, compress))
   return N.sum(temp, 1)/compress

############################################################
def expand_array(array, expand, sample=0):
   """
   Expands an 1-D array by the integer factor "expand".  
   if 'sample' is 1 the new array is created with sampling, if 1 then
   the new array is created via interpolation (default)
   Temporary fix until the equivalent of IDL's 'rebin' is found.
   """

   l = len(array)
   if (expand == 1): return array
   if (sample == 1): return N.repeat(array, expand)

   kernel = N.ones(expand)/expand
   # The following mimic the behavior of IDL's rebin when expanding
   temp = N.convolve(N.repeat(array, expand), kernel, mode=2)
   # Discard the first "expand-1" entries
   temp = temp[expand-1:]
   # Replace the last "expand" entries with the last entry of original
   for i in range(1,expand): temp[-i]=array[-1]
   return temp

