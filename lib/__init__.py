#/usr/bin/python
# Name:      tdl: tiny data language 
# Purpose:   small, embeddable data processing language written in python
# Author:    Matthew Newville
# Copyright: Matthew Newville, The University of Chicago, 2005
# Licence:   Python
# Created:   2005-Nov-15
#-----------------------------------------------------------------------------

"""
Tiny Data Language
"""
import version
__version__ = version.version

import os
import sys
import types

import Num
import Eval
import Expression
import Help
import Util
import Symbol
import Plotter
import TdlBuiltins
import TdlNumLib
import Shell
shell = Shell.shell
TdlBuiltins = TdlBuiltins
TdlNumLib = TdlNumLib
Plotter = Plotter
