#!/usr/bin/python
# M. Newville Univ of Chicago (2005)
#
# --------------
# Modifications
# --------------
#  12-Mar-2007  MN  simplified by *requiring* numpy 1.0 or higher, and
#                   not using scipy
#  23-Apr-2006  MN  added test for scipy 0.4.8 or higher, required for
#                   some numeric functionality
#
#  21-Mar-2006  MN  switched to requiring NumPy 0.9.6 or higher
#
##########################################################################
from __future__ import division

try:
    import numpy as Num
    factor = 100.
    version = 0.0
    for v in [int(i) for i in Num.__version__.split('.')]:
        version = version + v * factor
        factor = factor / 10.

    if version >= 100:   num_version = 'numpy %s' % Num.__version__
except:
    raise ImportError, 'Need numpy version 1.0 or higher'

from numpy import *
ArrayType = ndarray

if __name__ == '__main__':
    print 'tdl will use numpy %s ' % num_version
    
