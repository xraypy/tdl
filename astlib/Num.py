#!/usr/bin/python
# M. Newville Univ of Chicago (2005)
#
# --------------
# Modifications
# --------------
#  12-Mar-2007  MN  simplified by *requiring* numpy 1.0 or higher, and
#                   not using scipy here.
#  23-Apr-2006  MN  added test for scipy 0.4.8 or higher, required for
#                   some numeric functionality
#
#  21-Mar-2006  MN  switched to requiring NumPy 0.9.6 or higher
#
##########################################################################
from __future__ import division

def ParseVersion(s):
    """ parse a version string to an integer '1.0.1' -> 101"""
    factor = 100.
    version = 0.0
    for v in [int(i) for i in s.split('.')]:
        try:
            version = version + v * factor
        except:
            pass
        factor = factor * 0.1000
    return version
    
def NoNumpy():
    raise ImportError, 'Need numpy version 1.1 or higher'    

def NoScipy():
    has_scipy = False

try:
    import numpy as Num
    version = ParseVersion(Num.__version__)
    if version < 110:  NoNumpy()
    num_version = 'numpy %s' % Num.__version__
except:
    NoNumpy()

try:
    import scipy
    has_scipy = True
    version = ParseVersion(scipy.__version__)
    if version < 50: NoScipy()
    num_version = '%s with scipy %s ' % (num_version,scipy.__version__)
except:
    NoScipy()

from numpy import *
ArrayType = ndarray

if __name__ == '__main__':
    print 'tdl will use %s ' % num_version
