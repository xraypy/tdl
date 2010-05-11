print '== Symbol Test 01'

import larch
print '== Import OK... create SymbolTable....'

s = larch.SymbolTable()
print '== OK!'

print '== Initialized OK.  Search Groups:'

for i in s._sys.searchGroups:
    print '   %s ' % i

print '== Test getVariable, value '
print '  getSymbol(_sys.searchGroups) = ', s.getSymbol('_sys.searchGroups')

print '== Test addTempGroup '
for   i in range(3):
    print 'add temp group: ', s.createGroup()

print '== Groups: Name      elements    SearchGroup?'

allgroups = s._subgroups()
for i in s._sys.searchGroups:
    print '     %s         %i        yes' % ((i+' '*10)[:10], len(s.getSymbol(i)))
    allgroups.remove(i)

for i in allgroups:
    print '     %s         %i        no' % ((i+' '*10)[:10], len(s.getSymbol(i)))    
    
print '== Passed all tests for Symbol Test 01 '
#s.showTable()
