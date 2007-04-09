from Symbol import SymbolTable, symTypes, Group, Symbol
import numpy

main = Group(name='g1')
main['u'] = Symbol(name='u',value=15)
main['subgroup1'] = Group(name='subgroup1')
main['subgroup1']['x'] = Symbol(name='x', value='a string')
main['subgroup1']['y'] = Symbol(name='y', value=99)
main['v']   = Symbol(name='v', value=22.3)

s = SymbolTable()
s.setSymbol('_main.a', value = 1)
s.addGroup('g1')

sym  = s.getSymbolLocalGroup('g1.y')

s.showTable()


