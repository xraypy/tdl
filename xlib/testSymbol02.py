
import Eval

tdl = Eval.Evaluator()
s = tdl.symbolTable

s.addSymbol('_main.a', value = 1)
s.addGroup('g1')

s.addSymbol('_sys.var', value=1)
s.addSymbol('_main.b',  value=2)
s.addSymbol('_main.g1.s', value='a string')
s.addSymbol('_main.g1.t', value='another string')

s.addGroup('g2')
s.addSymbol('g2.data',value=[1,2,3])
s.addTempGroup()
s.addTempGroup()

print 'Groups: ' 
for g in s.listGroups():
    tdl.help.show_group(g)


