#!/usr/bin/env python
"""
"""

import re

class EvalError(Exception):
    def __init__(self,error,descr = None,node = None):
        self.error = error
        self.descr = descr
    def __repr__(self):
        return "%s: %s" % (self.error, self.descr)
    __str__ = __repr__


class closure:
    """Give a reference to a function with arguments so that it 
    can be called later, optionally changing the argument list.  

    The class provids a simple callback function which is then
    executed when called as a function. It can be defined as:

       >>>def my_action(x=None):
       ...        print 'my action: x = ', x
       >>>c = closure(my_action,x=1)
  
    and used as:
       >>>c()
       my action: x = 1
       >>>c(x=2)
        my action: x = 2

    The code is based on the Command class from
    J. Grayson's Tkinter book.
    """
    def __init__(self,func=None, **kw):
        self.func  = func
        self.kw    = kw
        
    def __repr__(self):
        return "<function %s>" % (self.func.__name__)
    __str__ = __repr__

    def __call__(self,  *args, **kw):
        if self.func is None: return None
        kwds = {}
        # avoid overwriting self.kw here!!
        for k,v in self.kw.items():  kwds[k] = v
        kwds.update(kw)
        return self.func(*args,**kwds)
    


def strip_comments(s,char='#'):
    "find character in a string, skipping over quoted text"
    if s.find(char) < 0: return s
    
    i = 0    
    while i < len(s):
        t = s[i]
        if t in ('"',"'"):
            x = s[i+1:].find(t)
            if x > 0:
                i = i + x
        elif t == char:
            return s[:i].rstrip()
        i = i + 1
    return s


reserved_words = ('and','as','assert','break','continue','def','del',
                  'elif','else','except','finally', 'for', 'from',
                  'if', 'import', 'in', 'is','not', 'or','pass',
                  'print', 'raise', 'return','try', 'while', 'group',
                  'end','endwhile','endif','endfor','endtry','enddef',
                  'True','False','None')
py_words = ('class', 'yield','with','lambda','exec')

name_match = re.compile(r"[a-z_][a-z0-9_]*(.[a-z_][a-z0-9_]*)*$").match

def isValidName(name):
    n = name[:].lower()
    if n in reserved_words: return False
    return name_match(n) is not None


def isNumber(x):
    try:
        y = complex(x)
        return True
    except:
        return False

def isLiteralStr(x):
    return ( (x.startswith("'") and x.endswith("'")) or
             (x.startswith('"') and x.endswith('"')) )
