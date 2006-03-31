#!/usr/bin/python
# M. Newville Univ of Chicago (2005)
#
# --------------
# Modifications
# --------------
#  21-Mar-2006  MN  switched to requiring NumPy 0.9.6 or higher
#
##########################################################################
from __future__ import division

version_needed = 'Need NumPy version 0.9.6 or higher'
try:
    import numpy as Num
    v = [int(i) for i in Num.__version__.split('.')]
    if v[1]<9 or v[2]<6: raise ImportError, version_needed
except:
    raise ImportError, version_needed
    

if __name__ == '__main__':
    print 'using numpy version ', Num.__version__
