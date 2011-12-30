"""
Wrapper for the the test fcns in _hello.dll

Authors/Modifications:
----------------------
* T. Trainor (tptrainor@alaska.edu)

Notes:
------

"""
#######################################################################

import os, sys
import numpy as num
import ctypes as C

#######################################################################

#Get the name of the 'lib' directory
libspath = os.path.dirname(__file__)
libspath = os.path.abspath(libspath)

# import the dll
if sys.platform == 'win32':
    hellodll = num.ctypeslib.load_library('_hello.dll',libspath)
else:
    hellodll = num.ctypeslib.load_library('_hello.so',libspath)

#######################################################################
def test_hello():
    x = num.array([[1.1,2.1,3.1],[4.1,5.1,6.1]])
    y = num.array([-1.,-20.,-30.])
    dptr = C.POINTER(C.c_double)
    xptr = (dptr*len(x))(*[row.ctypes.data_as(dptr) for row in x])
    yptr = y.ctypes.data_as(dptr)
    nr   = len(x)
    nc   = len(x[0])
    hellodll.wrhello(xptr,yptr,nr,nc)

#######################################################################
if __name__ == "__main__":
    test_hello()
    print dir()
    
