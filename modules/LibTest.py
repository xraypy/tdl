import os
import sys

def test(a,b):
    print a, type(a)
    print b, type(b)
    print a + b
    return


def test_tdl(a,b,tdl=None):
    print a, type(a)
    print b, type(b)
    return


_func_ = {"test.test":(test,None)}
_var_  = {"x":[1,2,3],"test.dat":[1.3,6.5]}


