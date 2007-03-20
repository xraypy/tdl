
import Eval

import numpy

tdl = Eval.Evaluator()
s = tdl.symbolTable

x1 = numpy.arange(60)*1.20
x1.shape = (5,3,4)

s.setSymbol('_main.a', value = 1)
s.addGroup('g1')
s.setSymbol('_sys.var', value=1)
s.setSymbol('_main.b',  value=2)

s.addGroup('mymod',toplevel=True)

xdict = {'a': 3, 'b':4}
s.addGroup('subgroup')
s.addGroup('subgroup.s2')
s.setSymbol('subgroup.s2.x',value='x1')
s.setSymbol('subgroup.bb',value=22.312)
s.setSymbol('subgroup.string1',value='another string')
s.setSymbol('subgroup.dict1',value=xdict)

s.LocalGroup='mymod'
s.setSymbol('mymod.a' , value=33.)
s.setSymbol('ooo' ,     value='aa' )
s.ModuleGroup='_main'

s.setSymbol('subgroup.foo', value=x1)

s.setSymbol('y', value = 'value for y')

s.setSymbol('boo', value = 12)

tdl.help.show_group('subgroup')
    


# print '========================'
# for g in s.listGroups():
#     if g not in ('_plot', '_math','_builtin'):
#         tdl.help.show_group(g)


# print s.getSymbol('subgroup')

# s.delGroup('subgroup.s2')
# s.delSymbol('subgroup')
