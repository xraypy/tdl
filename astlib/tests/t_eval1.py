#!/usr/bin/env python
# test expression parsing
from __future__ import print_function

import sys
sys.path.insert(0,'../')

import compiler
import numpy
import unittest
import types

def myf1(a,b,c=10):  return (a + b)/ c

def myf2(c=100,**kw): return c/2.0

def myf3(x,b):        return x * b

eval_tests  = open('eval1.txt','r').readlines()

floattypes   = numpy.sctypes['float']    + [types.FloatType]
complextypes = numpy.sctypes['complex']  + [types.ComplexType]
inttypes     = numpy.sctypes['int']      + numpy.sctypes['uint'] + [types.IntType]


class EvalTest(unittest.TestCase):
    def setUp(self):
        self.compiler = compiler.Compiler()
        self.eval = self.compiler.eval
        s  = self.compiler.symtable

        b = numpy.arange(64)
        b.shape=(8,8)

        s.setSymbol('a',numpy.arange(7.))
        s.setSymbol('b',b)
        s.setSymbol('c',numpy.arange(10.)/10)
        s.setSymbol('u',0.5)
        s.setSymbol('x',3)
        s.setSymbol('y',2)
        s.setSymbol('z',7.4)
        s.setSymbol('alist',[[1.,2.,3.,4.],[6.,7.,8.,9.],['a','b','c','d']])
        s.setSymbol('adict',{'a':1.0, 'b':2.0, 'c':3.0})
        s.setSymbol('yes',True)
        s.setSymbol('no', False)
        s.setSymbol('nottrue',False)
        s.setSymbol('orx',2)
        s.setSymbol('andy',99)
        s.setSymbol('f1', myf1)
        s.setSymbol('f2', myf2)
        s.setSymbol('f3', myf3)
        self.nlines = 0
        print("Setup done.")

    def testA(self):
        self.nlines = 0
        for line in eval_tests:
            line = line[:-1].strip()
            if line.startswith('#') or len(line)<2: continue            
            if line.find('=>') > 3:
                self.nlines = self.nlines + 1
                expr, val= line.split('=>')
                print( 'eval ', expr,end='')
                ast = self.compiler.compile(expr)
                out = self.compiler.interp(ast)

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
                elif type(out) in (numpy.ndarray,):
                    oshape = out.shape
                    out = out.ravel()
                    val = val.replace('[',' ').replace(']', ' ').strip().split()
                    for v,o in zip(val,out):
                        self.assertAlmostEqual(o,float(v),6)
                elif type(out) == types.BooleanType: 
                    self.assertEqual(repr(out),(val.strip()))
                elif type(out) == types.DictType:
                    # print 'Dict: ', repr(out)
                    # print 'Dict : ', val.strip()
                    self.assertEqual(repr(out),(val.strip()))
                else:
                    self.assertEqual(repr(out),(val.strip()))

                print(' OK.')
        print(' Evaluated %i expressions successfully.'  % self.nlines )

if __name__ == '__main__':
    unittest.main()

    
