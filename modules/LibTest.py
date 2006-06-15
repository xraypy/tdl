import os
import sys

def test(file_path=None):
    if hasattr(open,'file_path'):
        open.file_path = file_path
    print open.file_path
    return


def test_tdl(a,b,tdl=None):
    print a, type(a)
    print b, type(b)
    return


_func_ = {"test.test":(test,None)}
_var_  = {"x":[1,2,3],"test.dat":[1.3,6.5]}


