print 'Symbol Tests'

import Symbol
print 'imported OK'
from Symbol import SymbolTable


s = SymbolTable()
s.initialize()

print "Groups: "
print s.groups.keys()
print s.searchGroups

for i in s.searchGroups:
    print '== ', i
    print s.groups[i].keys()

print s.getVariable('_sys.searchGroups')
print s.getVariable('_sys.searchGroups').value

for   i in range(10): s.addRandomGroup()


print s.groups.keys()
