#!/usr/bin/python
# test expression parsing  & evaluation with a SymbolTable
#
from tdl.Evaluator import Evaluator
from tdl.Symbol import SymbolTable
import numpy as num

def myf1(a,b,c=10):
    return (a + b)/ c

def myf2(c=100,**kw):
    return c/2.0
def myf3(x,b):
    return x * b

s = SymbolTable()
p = Evaluator(symbolTable=s)

b = num.arange(64)
b.shape=(8,8)

s.addVariable('a',num.arange(7.))
s.addVariable('b',b)
s.addVariable('c',num.arange(10.)/10)
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

array = num.array

t = """
1+2*3                    =>  7                => int
-pi/3                    => -1.0471975512     => float
2*(4/3)                  =>  2.6666666667     => float
2*4/3                    =>  2.6666666667     => float
2*(4+3)                  =>  14               => int
2*4+3                    =>  11               => int
5/2+3                    =>  5.5              => float
5/2*3                    =>  7.5              => float
6/2/3                    =>    1              => int
2^3+3                    =>    11             => int
2^(3+3)                  =>   64              => int
2^(3*3)                  =>  512              => int
+3                       =>   3               => int
3 * -2                   =>  -6               => int
-3^2                     =>  -9               => int
2^3^2                    =>  512              => int
2^4 + 1                  =>  17               => int
8*3^2 + 2                =>  74               => int
-1**3                    =>  -1.00            => float
cos(0.1)                 => 0.9950041653      => float
sin(0.1)                 => 0.0998334166      => float
sin(pi/2 - tan(0.1))     => 0.9949706981      => float
tan(0.2 + 1.0j)          => (0.0831511848,0.7744312669j) => complex
e^(1j*pi) + 1            => (0.,0.j)          => complex
a[2] + 3                 =>   5               => float
b[2][4]                  => 20                => float
2j                       =>  (0.,2.j)         => complex
sqrt(33.2j / cos(88.0))  => (4.0755870761,4.0755870761j) => complex
x == 2                   =>  False    =>  bool
x == 3                   =>  True     =>  bool
x != 3                   =>  False    =>  bool
x >= 1                   =>  True     =>  bool
x >= 3                   =>  True     =>  bool
x >= 4                   =>  False    =>  bool
x < 99                   =>  True     =>  bool
(x>u) or x*u>100         =>  True     =>  bool
yes and no or nottrue     =>  False   => bool
yes and (no or nottrue)   =>  False   => bool
(yes and no) or nottrue   =>  False   => bool
yes or no and nottrue     =>  True    => bool
yes or (no and nottrue)   =>  True    => bool
(yes or no) and nottrue   =>  False   => bool
yes or not no             =>  True    => bool
not (no or yes)           =>  False   => bool
not no or yes             =>  True    => bool
not yes                   =>  False   => bool
not no                    =>  True    => bool
! yes                     =>  False   => bool
! no                      =>  True    => bool
a[0]                      =>  0.      => float
a[1]                      =>  1.      => float
a[0] and a[1]             =>  0.      => float
not a[0]                  =>  True    => bool
not a[1]                  =>  False   => bool
not a[2]                  =>  False   => bool
f1(2,8.0)                 =>  1.      => float
f1(2,2*7,c=11.0)          =>  1.4545454545 => float
f1(2,8.0,c=sin(0.3))      =>  33.8386336182 => float
f2()                      =>   50.0   => float
f2(c=99)                  =>  49.5    => float
f2(d=88)                  =>  50.0    => float
max(a[3]*01,9,3)          =>   9.0    => float
adict['b']                =>    2     => float
adict['c'] /adict['b']    =>  1.5     => float
0.2                       =>  0.2     => float
.2                        =>   .2     => float
0.223e+4                  => 0.223e+4 => float
.32101e+3j                => (0,.32101e+3j) => complex
f3(2.5,b)[3][4]           =>  70.0    =>float
alist[2][3]               =>  'd'     => str
'a simple string '        => ' a simple string ' => str
' %s = %g '  % ('sqrt(22.3)',sqrt(22.3))   => ' sqrt(22.3) = 4.72229 ' => str
adict                     =>  {'a': 1.0, 'c': 3.0, 'b': 2.0}  => dict
[1,2,['a', 'b','c'], 4.0] => [1.0, 2.0, ['a', 'b', 'c'], 4.0] => list
alist[2]                  =>  ['a', 'b', 'c', 'd']  => list
{'a':1, 'b':2.0, 'c':[1,2,3]} => {'a': 1.0, 'c': [1.0, 2.0, 3.0], 'b': 2.0} => list
alist                     =>   [[1.0, 2.0, 3.0, 4.0], [6.0, 7.0, 8.0, 9.0], ['a', 'b', 'c', 'd']]  => list
[0,1,2,4,8,16,32]         =>  [0. 1. 2. 4. 8. 16. 32.] => array
alist[1]                  =>  [6. 7. 8. 9.]  => array
f3(2.5,b)[3]              =>  [ 60.   62.5  65.   67.5  70.   72.5  75.   77.5]  => array
b[2:4]                    =>  [[16 17 18 19 20 21 22 23] [24 25 26 27 28 29 30 31]]  => array
-(a[3:5] + 2)             =>  [-5. -6.]  => array
-a[3:7] + 2.0             =>  [-1. -2. -3. -4.]  => array
sin(a*pi/20)              =>  [ 0.          0.15643447  0.30901699  0.4539905   0.58778525  0.70710678  0.80901699]  => array
"""

total, passed = 0,0
for s in t.split('\n'):
    s = s.strip()
    if len(s)<1 or s.startswith('#'): continue
    w = s.split('=>')
    expr,expval = w[0].strip(),w[1].strip()
    vtype = 'float'
    if w[2]: vtype = w[2].strip()
    outval = p.do_eval(expr)
    total = total + 1
    if vtype =='int':
        outval  = "%i" % int(outval)        
    elif vtype =='float':
        outval = "%.7f" % float(outval)
        expval  = "%.7f" % float(expval)
    elif vtype =='str':
        outval = outval.strip()
        expval  = expval.strip()[1:-1].strip()  # strip off quotes
    elif vtype =='list':
        outval = repr(outval).strip()
    elif vtype =='array':
        for ix in ('\r','\n','[',']'):
            outval = str(outval).replace(ix,'')
            expval = str(expval).replace(ix,'')
        outval = outval.split()
        expval = expval.split()
    elif vtype =='dict':
        outval = str(outval)
    elif vtype =='bool':
        outval = outval in (1,'1',True,'True')
        expval = expval  in (1,'1',True,'True')
    elif vtype =='complex':
        x = expval[1:-1].split(',')
        outval = "(%.7f,%.7fj)" % (outval.real, outval.imag)
        expval = "(%.7f,%.7fj)" % (float(x[0]),float(x[1][:-1]))

    if expval != outval:
        print 'Warning evaluating expression "%s"  type=%s' % (expr, vtype)
        print '    expected: ', expval
        print '    got:      ', outval
        print '-------------------------------------'
    else:
        passed = passed + 1


print "Passed %i of %i tests." % (passed , total)

    
newtests = """
# adict
# f3(2.5,b)[3]
# b[2:4]
# -(a[3:5] + 2)
# -a[3:7] + 2.0
# sin(a*pi/20)
"""
# 
# for s in newtests.split('\n'):
#     if len(s)>0 and not s.startswith('#'):
#         expr = s.strip()
#         val = p.do_eval(expr)
#         print expr, ' => ', str(val), ' => array'
# 

