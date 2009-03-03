#!/usr/bin/env python2.6

import numpy
import compiler
import inputText
from util import EvalError
import sys
compiler = compiler.Compiler()
input    = inputText.InputText(prompt='>',interactive=False)

text = '''x=3
x =x+1
print x'''
input.put(text)                
while len(input) >0:
    block,fname,lineo = input.get()
    if block is not None:
        compiler.eval(block)
        

