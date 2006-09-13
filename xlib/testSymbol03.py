import Symbol, time
s = Symbol.SymbolTable()

t0 = time.time()

for i in range(400):  s.addTempGroup()

t=time.time()-t0
sl =  s.data.keys()

print len(sl), t

t=time.time()-t0
s.clearTempGroups()
print len(s.data.keys()), t
