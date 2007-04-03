from Symbol import SymbolTable, symTypes, Group, Symbol

import numpy
main = Group(name='mg')
main['u'] = Symbol(name='u',value=15)
main['subgroup1'] = Group(name='subgroup1')
main['subgroup1']['x'] = Symbol(name='x', value='a string')
main['subgroup1']['y'] = Symbol(name='y', value=99)
main['v']   = Symbol(name='v', value=22.3)

s = SymbolTable()
s.setSymbol('_main.a', value = 1)
s.addGroup('g1')

sym  = s.getSymbolLocalGroup('g1.y')
sym  = s.getSymbolLocalGroup('g1.v')

s.showTable()

        
# s.setSymbol('_sys.var', value=1)
# s.setSymbol('_main.b',  value=2)
# s.setSymbol('_main.g1.s', value='a string')
# s.setSymbol('_main.g1.t', value='another string')
# 
# s.addGroup('g2')
# s.addGroup('_sys.subgroup')
# s.addGroup('_sys.subgroup.subsubgroup')
# 
# s.setSymbol('_main.c', value = 44)
# s.setSymbol('g1.y', value = 'a')
# s.setSymbol('g2.y', value = 'b')
# s.setSymbol('g2.boo', value = 1222)
# 
# s.setSymbol('_sys.y', value = 'this is _sys.y')
# s.setSymbol('_main.b', value =6000)
# 
# # s.showTable()
# 
# print '==================='
# 
# print s.getSymbol('a')
# print s.getSymbol('c')
# print s.getSymbol('g1')
# 
# # print s.searchGroups
# 
# print s.getSymbol('g1.y')
# print s.getSymbol('_main.g1.t')
# 
# 
# # s.searchGroups.append('_main.g2')
#    
# # 
# print s.setSymbol('kexp', value='90.3')
# 
# s.setSymbol('_main.x',value = 'main x')
# s.setSymbol('_main.g1.x',value = 'main g1 x')
# s.setSymbol('_main.g2.x',value = 'main g2 x')
# 
# n = 0
# for i in ('_main', '_main.g1', '_main.g2'):
#     n = n + 1
#     s.LocalGroup = i
#     print ':: ', i, s.getSymbol('y').value
#     s.setSymbol('x', n)
# 
# 
# print s.searchGroups
# 
# x1 = numpy.arange(60)*1.20
# x1.shape = (5,3,4)
# s.LocalGroup = '_main'
# s.setSymbol('foo', value=x1)
# # 
 
# s.initialize()
# print ' ============== ' 

# 
# print ' ============== '
# 
# print s.getSymbol('foo')
# 
# print s.getSymbol('foo.shape')
# print s.getSymbol('foo.size')
# 
# s.addGroup('q')
# 
# s.addGroup('ww', toplevel=True)
# 
# print 'hasgroup ww' , s.hasGroup('ww')
# print 'hasgroup q' , s.hasGroup('q')
# print 'hasgroup _sys.q' , s.hasGroup('_sys.q')
# print 'hasgroup bob' , s.hasGroup('bob')
# 
# 
# 
# 
# s.showTable()
# 
# # print 'TopLevel Groups: '
# # for k,v in  s.data.items():
# #     if type(v)==Group:    print k
# # 
# # 
# # print 'Groups under _main: '
# # for k,v in  s.data['_main'].items():
# #     if type(v)==Group: print k
# # 
# 
# # 
# # 
# # print s.hasSymbol('_main.x')
# # print s.hasSymbol('x')
# # print s.hasSymbol('ww')
# # print s.hasFunction('sort')
