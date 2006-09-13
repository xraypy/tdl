from Symbol import SymbolTable, symTypes, symGroup, Symbol

import numpy
x1 = numpy.arange(60)*1.20
x1.shape = (5,3,4)


s = SymbolTable()

s.addSymbol('_main.a', value = 1)
s.addGroup('g1')
s.addSymbol('_sys.var', value=1)
s.addSymbol('_main.b',  value=2)

s.addGroup('mymod',toplevel=True)
s.addGroup('_sys.subgroup')
s.addGroup('_sys.subgroup.s2')
s.addSymbol('_sys.subgroup.s2.x',value='x1')
s.addSymbol('_sys.subgroup.bb',value=22)
s.LocalGroup='mymod'
s.addSymbol('mymod.a' , value=33.)
s.addSymbol('ooo' ,     value='aa' )
s.ModuleGroup='_main'

s.addSymbol('_sys.subgroup.foo', value=x1)


s.showTable()

s.addSymbol('y', value = 'value for y')

s.addSymbol('boo', value = 12)

s.showTable()


print s.getSymbol('subgroup')

# s.delGroup('subgroup.s2')
# s.delSymbol('subgroup')
print 'Groups '
print s.listGroups()
