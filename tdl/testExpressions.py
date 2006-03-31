#!/usr/bin/env python
# test expression parsing

from Num import Num
from Eval       import Evaluator
from Expression import ExpressionParser
from Symbol     import SymbolTable

import unittest
import types

def myf1(a,b,c=10):
    return (a + b)/ c

def myf2(c=100,**kw):
    return c/2.0

def myf3(x,b):
    return x * b

parse_file = '../tests/core_parsing.txt'
f = open(parse_file,'r')
parse_tests = f.readlines()
f.close()

eval_file = '../tests/core_eval.txt'
f = open(eval_file,'r')
eval_tests = f.readlines()
f.close()


floattypes   = Num.sctypes['float']    + [types.FloatType]
complextypes = Num.sctypes['complex']  + [types.ComplexType]
inttypes     = Num.sctypes['int']      + Num.sctypes['uint'] + [types.IntType]

class ParserTest(unittest.TestCase):
    def setUp(self):
        self.symtable  = SymbolTable()
        self.parser    = ExpressionParser()
        self.evaluator = Evaluator(symbolTable=self.symtable)
        s = self.symtable

        b = Num.arange(64)
        b.shape=(8,8)

        s.addVariable('a',Num.arange(7.))
        s.addVariable('b',b)
        s.addVariable('c',Num.arange(10.)/10)
        s.addVariable('u',0.5)
        s.addVariable('x',3)
        s.addVariable('y',2)
        s.addVariable('z',7.4)
        s.addVariable('alist',[[1.,2.,3.,4.],[6.,7.,8.,9.],['a','b','c','d']])
        s.addVariable('adict',{'a':1.0, 'b':2.0, 'c':3.0})
        s.addVariable('yes',True)
        s.addVariable('no', False)
        s.addVariable('nottrue',False)
        s.addVariable('orx',2)
        s.addVariable('andy',99)
        s.addFunction('f1', myf1)
        s.addFunction('f2', myf2)
        s.addFunction('f3', myf3)
        self.nlines = 0

    def __testA_Load(self):
        out = self.evaluator.eval('1')
        self.assertEqual(out,1)
        
    def testB_ExpressionCompile(self):
        print 'Parsing Test:' 
        self.nlines = 0
        for line in parse_tests:
            line = line[:-1].strip()
            if line.startswith('#') or len(line)<2: continue
            if line.find('=>') > 3:
                self.nlines = self.nlines + 1
                expr, val= line.split('=>')
                val.strip()
                test  = val.split(',')
                print 'parsing ', expr, 
                stack = self.parser.compile(expr)
                stack.reverse()
                # print ' :: got     = ', stack
                # print ' :: expected= ', test
                for i in range(len(stack)):
                    out = stack[i]
                    val  = test[i]
                    if type(out) in floattypes:
                        self.assertAlmostEqual(out,float(val.strip()),6)

                    elif type(out) in inttypes:
                        self.assertEqual(str(out),str(val.strip()))
                        
                    elif type(out) in complextypes:
                        val = complex(val.strip())
                        self.assertAlmostEqual(out.real,val.real,6)
                        self.assertAlmostEqual(out.imag,val.imag,6)

                    elif type(out) == types.StringType:
                        self.assertEqual(str(out.strip()),str(val.strip()))                        
                    else:
                        print ' Unknown type for element ', out, ' in ', expr, type(out)
                        self.assertEqual(out,val)
                print ' OK.'
        print ' Compiled %i expressions successfully.'  % self.nlines

    def testC_ExpressionEval(self):
        print 'Evaluation Test:'
        self.nlines = 0
        for line in eval_tests:
            line = line[:-1].strip()
            if line.startswith('#') or len(line)<2: continue            
            if line.find('=>') > 3:
                self.nlines = self.nlines + 1
                expr, val= line.split('=>')
                print 'eval ', expr, 
                out = self.evaluator.eval(expr)
                if type(out) in floattypes:
                    self.assertAlmostEqual(out,float(val.strip()),6)
                elif type(out) in inttypes:
                    self.assertEqual(str(out),str(val.strip()))
                elif type(out) in complextypes:
                    val = val.strip()
                    xr,xi = val[1:-1].split(',')
                    self.assertAlmostEqual(out.real,float(xr),6)
                    self.assertAlmostEqual(out.imag,float(xi[:-1]),6)
                elif type(out) == types.StringType:
                    self.assertEqual("'%s'"%str(out),str(val.strip()))
                elif type(out) in (Num.ndarray,):
                    oshape = out.shape
                    out = out.ravel()
                    val = val.replace('[',' ').replace(']', ' ').strip().split()
                    for v,o in zip(val,out):
                        self.assertAlmostEqual(o,float(v),6)
                elif type(out) == types.BooleanType: 
                    self.assertEqual(repr(out),(val.strip()))
                else:
                    self.assertEqual(repr(out),(val.strip()))                    
                print ' OK.'
        print ' Evaluated %i expressions successfully.'  % self.nlines                

if __name__ == '__main__':
    unittest.main()

    
