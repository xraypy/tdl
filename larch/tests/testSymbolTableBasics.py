print '== Symbol Test 01'

import larch
import larch.symboltable
print '== Import OK... create SymbolTable....'

s = larch.SymbolTable()
print '== OK!'

print '== Initialized OK.  Search Groups:'

for i in s._sys.searchGroups:
    print '   %s ' % i

print '== Test getVariable, value '
print '  getSymbol(_sys.searchGroups) = ', s.getSymbol('_sys.searchGroups')

print '== Test addTempGroup '
for   i in range(5):
    g = s.createGroup(name="g%i" % i, x=i, y=2+i/3.0)    
    print 'add group: ', s.setSymbol('g%i' % i, g)

print '== Groups: Name      elements    SearchGroup?'

allgroups = s._subgroups()
for i in s._sys.searchGroups:
    print '     %s         %i        yes' % ((i+' '*15)[:15], len(s.getSymbol(i)))
    allgroups.remove(i)

for i in allgroups:
    sym = s.getSymbol(i)
    # print i, larch.symboltable.isgroup(sym)
    if larch.symboltable.isgroup(sym):        
        print '     %s         %i        no' % ((i+' '*15)[:15], len(sym._members()))
    else:
        print '     %s         no' % ((i+' '*15)[:15])

print s.show_group('g3')

print '== Passed all tests for Symbol Test 01 '
#s.showTable()
