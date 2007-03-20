
import Eval

tdl = Eval.Evaluator()
s = tdl.symbolTable

s.setSymbol('_main.a', value = 1)
s.addGroup('g1')

s.setSymbol('_sys.var', value=1)
s.setSymbol('_main.b',  value=2)
s.setSymbol('_main.g1.s', value='a string')
s.setSymbol('_main.g1.t', value='another string')

s.addGroup('g2')
s.addGroup('_main.g3')
s.setSymbol('g2.data',value=[1,2,3])
s.addTempGroup()
s.addTempGroup()

print '>>>>>>>Groups: ' 
# for g in s.listGroups():
#    tdl.help.show_group(g)
s.showTable(skip=['_builtin','_math'])


