#!/usr/bin/python
# M. Newville Univ of Chicago (2005)
#
# --------------
# Modifications
# --------------
#  23-Apr-2006  MN  added test for scipy 0.4.8 or higher, required for
#                   some numeric functionality
#
#  21-Mar-2006  MN  switched to requiring NumPy 0.9.6 or higher
#
##########################################################################
from __future__ import division

numpy_needed = 'Need SciPy 0.4.8 or (at least) NumPy version 0.9.6 or higher'
num_version = None
try:
    import scipy as Num
    v = [int(i) for i in Num.__version__.split('.')]
    v = 100*v[0] + 10*v[1] + v[2]
    if v >= 47:
        num_version = 'scipy %s' % Num.__version__
except:
    try:
        import numpy as Num
        v = [int(i) for i in Num.__version__.split('.')]
        v = 100*v[0] + 10*v[1] + v[2]
        if v >= 96:
            num_version = 'numpy %s' % Num.__version__            
    except:
        pass

if num_version is None:
    raise ImportError, numpy_needed

if __name__ == '__main__':
    print 'tdl will use %s ' % num_version
    
