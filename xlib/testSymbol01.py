print '== Symbol Test 01'

import Symbol
print '== Import OK... create Table....'

s = Symbol.SymbolTable(init=False)
print '== OK!'

print '== Initialization.... '
s.initialize()
print '== Initialized OK.  Search Groups:'

for i in s.searchGroups:
    print '   %s ' % i

print '== Test getVariable, value '
print '  getVariable(_sys.searchGroups) = ', s.getVariable('_sys.searchGroups')
print '  value of _sys.searchGroups     = ', s.getVariable('_sys.searchGroups').value

print '== Test addTempGroup '
for   i in range(3):
    print 'add temp group: ', s.addTempGroup()

print '== Groups: Name      elements    SearchGroup?'

allgroups = s.listGroups()
for i in s.searchGroups:
    print '     %s         %i        yes' % ((i+' '*10)[:10], len(s.getSymbol(i)))
    allgroups.remove(i)

for i in allgroups:
    print '     %s         %i        no' % ((i+' '*10)[:10], len(s.getSymbol(i)))    
    
print '== Passed all tests for Symbol Test 01 '
#s.showTable()
