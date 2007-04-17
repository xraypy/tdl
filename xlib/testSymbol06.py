from Symbol import SymbolTable, symTypes, Group, Symbol, isSymbol, isGroup
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
s.addGroup('g1.g2')
s.setSymbol('g1.x', value = 22)
s.setSymbol('g1.y', value = 'a string')
s.setSymbol('g1.g2.doc', value = 'nested string')

print 'Groups: '
print s.listGroups()

g1 = s.getGroup('_main')

print g1, g1.stats()
           
def show_group(g,indent=0):
    if not isGroup(g): return None
    tab = '   '*indent
    for k,v in g.items():
        print "%s %s %s" %(tab, k,  v.getinfo())
        if isGroup(v): show_group(v,indent=indent+1)


show_group(s.getGroup('_main'),indent=1)
