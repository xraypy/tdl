#!/usr/bin/env python
"""utilities for larch
"""
from __future__ import print_function
import re
import sys
def PrintExceptErr(err_str, print_trace=True):
    " print error on exceptions"
    print('\n***********************************')
    print(err_str)
    #print 'PrintExceptErr', err_str
    try:
        print('Error: %s' % sys.exc_type)
        etype, evalue, tback = sys.exc_info()
        if print_trace == False:
            tback = ''
        sys.excepthook(etype, evalue, tback)
    except:
        print('Error printing exception error!!')
        raise
    print('***********************************\n')

class Closure:
    """Give a reference to a function with arguments so that it 
    can be called later, optionally changing the argument list.  

    The class provids a simple callback function which is then
    executed when called as a function. It can be defined as:

       >>>def my_action(x=None):
       ...        print 'my action: x = ', x
       >>>c = Closure(my_action,x=1)
  
    and used as:
       >>>c()
       my action: x = 1
       >>>c(x=2)
        my action: x = 2

    The code is based on the Command class from
    J. Grayson's Tkinter book.
    """
    def __init__(self, func=None, **kwds):
        self.func = func
        self.kwds = kwds

    def __repr__(self):
        return "<function %s>" % (self.func.__name__)
    __str__ = __repr__

    def __call__(self, *args, **c_kwds):
        if self.func is None:
            return None
        # avoid overwriting self.kwds here!!
        kwds = {}
        for key, val in list(self.kwds.items()):
            kwds[key] = val
        kwds.update(c_kwds)
        return self.func(*args, **kwds)

def strip_comments(sinp, char='#'):
    "find character in a string, skipping over quoted text"
    if sinp.find(char) < 0:
        return sinp    
    i = 0    
    while i < len(sinp):
        tchar = sinp[i]
        if tchar in ('"',"'"):
            eoc = sinp[i+1:].find(tchar)
            if eoc > 0:
                i = i + eoc
        elif tchar == char:
            return sinp[:i].rstrip()
        i = i + 1
    return sinp


RESERVED_WORDS = ('and', 'as', 'assert', 'break', 'continue', 'def',
                  'del', 'elif', 'else', 'except', 'finally', 'for',
                  'from', 'if', 'import', 'in', 'is', 'not', 'or',
                  'pass', 'print', 'raise', 'return', 'try', 'while',
                  'group', 'end', 'endwhile', 'endif', 'endfor',
                  'endtry', 'enddef', 'True', 'False', 'None')

NAME_MATCH = re.compile(r"[a-z_][a-z0-9_]*(.[a-z_][a-z0-9_]*)*$").match

def isValidName(name):
    "input is a valid name"
    tnam = name[:].lower()
    if tnam in RESERVED_WORDS:
        return False
    return NAME_MATCH(tnam) is not None

def isNumber(num):
    "input is a number"
    try:
        cnum = complex(num)
        return True
    except ValueError:
        return False

def isLiteralStr(inp):
    "is a literal string"
    return ((inp.startswith("'") and inp.endswith("'")) or
            (inp.startswith('"') and inp.endswith('"')))
