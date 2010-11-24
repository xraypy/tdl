#!/usr/bin/env python
#
# this attempts to run a tdl script with the python interpreter.
# it will definitely fail for many cases (including 'x.y' names!!)
# but is useful for testing syntax and results.
# 

from numpy import *
import sys
import os
import glob

enddef = endif = endfor = endwhile = None
ls = glob.glob

if len(sys.argv)>1: execfile(sys.argv[1])
    

