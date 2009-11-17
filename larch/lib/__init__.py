#!/usr/bin/python2.6


"""
  Larch: a data processing macro language for python
  
"""

__version__  = '0.8.0'
__date__     = '16-Nov-2009'

##
## require that numpy be available right away!!
import numpy

import interpreter
from symbolTable import Group, symbolTable
from shell import shell

larch = interpreter.Interpreter

