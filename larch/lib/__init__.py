#!/usr/bin/env python2.6

"""
  Larch: a data processing macro language for python
"""

#
import sys
major, minor = sys.version_info[0], sys.version_info[1]
if major == 2 and minor < 6:
    raise EnvironmentError('requires python 2.6 or higher')

## require that numpy be available right away!!
import numpy

from . import interpreter
from .symbolTable import Group, symbolTable
from .shell import shell
from .interpreter import Interpreter
from .inputText import InputText

interp= Interpreter
input = InputText

__version__ = interpreter.__version__
__date__    = '22-Nov-2009'
