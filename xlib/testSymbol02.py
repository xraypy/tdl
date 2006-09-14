from Symbol import SymbolTable, symTypes, symGroup, Symbol

import numpy
s = SymbolTable()

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


s.showTable(skip=('_math','_builtin'))

s.clearTempGroups()

s.showTable(skip=('_math','_builtin'))
# s.addGroup('_sys.subgroup')
# s.addGroup('_sys.subgroup.subsubgroup')
# 
# s.addSymbol('_main.c', value = 44)
# s.addSymbol('g1.y', value = 'a')
# s.addSymbol('g2.y', value = 'b')
# s.addSymbol('g2.boo', value = 1222)
# 
# s.addSymbol('_sys.y', value = 'sys y')
# s.addSymbol('_main.b', value =6000)
# s.addSymbol('_main.x', value ='x')
# 
# # s.showTable()
# 
# x1 = numpy.arange(60)*1.20 ; x1.shape = (5,3,4)
# s.addSymbol('foo', value=x1)
# 
# print ' ============== ShowTable' 
# s.showTable()
# 
# print '*************************'
# 
# print s.getSymbol('foo')
# print s.getSymbol('b')
# print s.getSymbol('_main.c')
# 
# print s.getSymbol('foo.shape')
# print s.getSymbol('_main.foo.size')
# 
# s.addGroup('q')
# s.addGroup('mymodule', toplevel=True)
# 
# 
# # print 'TopLevel Groups: '
# # for k,v in  s.data.items():
# #     if type(v)==symGroup:    print k
# # 
# # 
# # print 'Groups under _main: '
# # for k,v in  s.data['_main'].items():
# #     if type(v)==symGroup: print k
# # 
# u =  s.getSymbol('mymodule')
# print u
# 
# print s.hasSymbol('_main.x')
# print s.hasSymbol('x')
# print s.hasSymbol('mymodule')
# print s.hasSymbol('mymodule.o')
# 
# s.addSymbol('mymodule.o', value =22)
# print s.getSymbol('mymodule.o')
# 
# print '======================'
# class bob:
#     foo = 22
#     x   = x1
# 
# s.addSymbol('_main.bob',value=bob)
# print s.getSymbol('_main.bob')
# 
# print '---'
# 
# print s.getSymbol('_main.bob.foo')
# print s.getSymbol('_main.bob.x')
# print s.getSymbol('_main.bob.x.shape')
# 
# print s.getSymbol('_main.bob.x.size')
# 
# # 
# # s.delGroup('_main.g1')
# # s.delGroup('mymodule')
# 
# # s.showTable()
# 
