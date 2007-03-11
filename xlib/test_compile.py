#!/usr/bin/env python
# test expression parsing

from Num import Num
from Eval       import Evaluator
from Expression import ExpressionParser
from Symbol     import SymbolTable

import types
import sys
expr = sys.argv[1]


floattypes   = Num.sctypes['float']    + [types.FloatType]
complextypes = Num.sctypes['complex']  + [types.ComplexType]
inttypes     = Num.sctypes['int']      + Num.sctypes['uint'] + [types.IntType]


symtable  = SymbolTable()
parser    = ExpressionParser()
evaluator = Evaluator(symbolTable=symtable)

symtable.setVariable('a',Num.arange(7.))
symtable.setVariable('pi',Num.pi)

print 'parse ', expr, 
stack = parser.compile(expr)
stack.reverse()
print ' :: got     = ', stack

