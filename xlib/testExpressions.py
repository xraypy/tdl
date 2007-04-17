#!/usr/bin/env python
# test expression parsing

import Num
import Eval

import unittest
import types

def myf1(a,b,c=10):
    return (a + b)/ c

def myf2(c=100,**kw):
    return c/2.0

def myf3(x,b):
    return x * b

parse_tests = open('../tests/core_parsing.txt','r').readlines()
eval_tests  = open('../tests/core_eval.txt','r').readlines()

floattypes   = Num.sctypes['float']    + [types.FloatType]
complextypes = Num.sctypes['complex']  + [types.ComplexType]
inttypes     = Num.sctypes['int']      + Num.sctypes['uint'] + [types.IntType]

class ParserTest(unittest.TestCase):
    def setUp(self):
        self.tdl = Eval.Evaluator()
        self.compiler = self.tdl.expr_compile
        s  = self.tdl.symbolTable

        b = Num.arange(64)
        b.shape=(8,8)

        s.setVariable('a',Num.arange(7.))
        s.setVariable('b',b)
        s.setVariable('c',Num.arange(10.)/10)
        s.setVariable('u',0.5)
        s.setVariable('x',3)
        s.setVariable('y',2)
        s.setVariable('z',7.4)
        s.setVariable('alist',[[1.,2.,3.,4.],[6.,7.,8.,9.],['a','b','c','d']])
        s.setVariable('adict',{'a':1.0, 'b':2.0, 'c':3.0})
        s.setVariable('yes',True)
        s.setVariable('no', False)
        s.setVariable('nottrue',False)
        s.setVariable('orx',2)
        s.setVariable('andy',99)
        s.setFunction('f1', myf1)
        s.setFunction('f2', myf2)
        s.setFunction('f3', myf3)
        self.nlines = 0

    def testA_Load(self):
        print 'Load Test:' 
        out = self.tdl.eval('1')
        self.assertEqual(out,1)
        
    def testB_ExpressionCompile(self):
        print 'Parsing Test:' 
        self.nlines = 0
        for line in parse_tests:
            line = line[:-1].strip()
            if line.startswith('#') or len(line)<2: continue
            if line.find('=>') > 3:
                self.nlines = self.nlines + 1
                expr, result = line.split('=>')
                out  = [i.strip() for i in result.split(',')]
                test = [repr(i) for i in self.compiler(expr)]

                print 'parse ', expr, 
                if test == out:
                    print 'OK.'
                else:
                    print 'failed: ',
                    for t, o in zip(test,out):
                        if t!=o: print 'expected element: ', o, ', got ', t

        print ' Compiled %i expressions successfully.'  % self.nlines

    def testC_ExpressionEval(self):
        print 'Evaluation Test ==== :'
        self.nlines = 0
        for line in eval_tests:
            line = line[:-1].strip()
            if line.startswith('#') or len(line)<2: continue            
            if line.find('=>') > 3:
                self.nlines = self.nlines + 1
                expr, val= line.split('=>')
                print 'eval ', expr,
                out = self.tdl.eval(expr)
                print ' = ', repr(out), 
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
                elif type(out) == types.DictType:
                    # print 'Dict: ', repr(out)
                    # print 'Dict : ', val.strip()
                    self.assertEqual(repr(out),(val.strip()))
                else:
                    self.assertEqual(repr(out),(val.strip()))

                print ' OK.'
        print ' Evaluated %i expressions successfully.'  % self.nlines                

if __name__ == '__main__':
    unittest.main()

    
