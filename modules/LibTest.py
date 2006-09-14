import os
import sys
import time

def test2(file_path=None):
    if hasattr(open,'file_path'):
        open.file_path = file_path
    print open.file_path
    return

def testfunc(x):
    print 'This is python function testfunc ', x
    return 22*x


def libtest_init(tdl=None):
    print '  This is Libtest: ', time.ctime()
    x = tdl.symbolTable.getVariable('_sys.searchGroups').value
    x.append('test')
    tdl.symbolTable.setVariable('_sys.searchGroups',x)
    return

# tell tdl to create these groups:
_groups_ = [('test',True)]

_func_ = {"test.test2":(test2,None),
          "test.testfunc": (testfunc,None)}

_var_  = {"x":[1,2,3],
          "test.dat":[1.3,6.5]}

# code to run on initialization (no args, but will get a 'tdl reference')
_init_ = libtest_init


