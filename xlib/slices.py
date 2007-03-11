from numpy import arange, ndarray
from types import FloatType, StringType, ListType
class opcodes:
    empty = '@EMP'
    
def test_eval(s):
    r = []
    r.append("==== %s:\n" % s)
    r.append(str(eval(s)))
    return ''.join(r)


def take_subarray(val,elems):
    # takes subarrays / slices of lists and numeric arrays,
    # including syntax like: x[:2], x[3:], x[2:5,3], x[2:40:2,4:8]
    # print 'This is take_subarray ', type(val),  elems, type(elems)
    
    el = []
    for e in elems:
        if type(e) == FloatType: e = int(e)
        el.append(e)
    if type(val) == ndarray:
        if len(el) > len(val.shape): el = el[:len(val.shape)]
        val = val[tuple(el)]
    else:
        for e in el: val = val[e]
    return val



x = arange(60,dtype='float64') ; x.shape=(5,4,3)
s = 'a long string here'
l = [1,2,'abcdefg',4,[5.,7.,9.,8,10.,'a']]

tests = [
    (l[2],         l, (2, ) ),
    (l[2][1:4],    l, (2,slice(1,4))),
    (l[4][4],      l, (4,4)),
    (l[4][1:3],    l, (4,slice(1,3))),     
    (s[1:4],       s, (slice(1,4),)),
    (s[:-1],       s, (slice(None,-1),)),
    (x[0:3],       x, (slice(0,3),)),
    (x[3,1,2],     x, (3,1,2)), 
    (x[3,2,1],     x, (3,2,1)), 
    (x[3,2,-1],    x, (3,2,-1)), 
    (x[3,-1,-1],   x, (3,-1,-1)),
    (x[3,2],       x, (3,2)),
    (x[-1,2,1],    x, (-1,2,1)), 
    (x[:,2],       x, (slice(None,None), 2)), 
    (x[0:4,2],     x, (slice(0,4), 2)),  
    (x[:,2,2],     x, (slice(None,None), 2,2)),
    (x[:,1:3],     x, (slice(None,None), slice(1,3))),
    (x[:,1:3,1],   x, (slice(None,None), slice(1,3),1)),
    (x[:,1:3,0:1], x, (slice(None,None), slice(1,3),slice(0,1))),
    (x[:,:,0:1],   x, (slice(None,None), slice(None,None),slice(0,1))),
    (x[:,:,1],     x, (slice(None,None), slice(None,None),1)),
    (x[0:3,2:3,],  x, (slice(0,3),    slice(2,3))),
    (x[0:3,0:2,2], x, (slice(0,3),    slice(0,2), 2)),
    (x[0:3,2,:],   x, (slice(0,3),    2, slice(None,None))),
    (x[0:3,2,],    x, (slice(0,3),    2))
]
# 
for t,v,s in tests:
    print '-> ', str(type(v))[6:-1], s, 
    assert(repr(t) == repr(take_subarray(v,s)))
    print 'OK. ' 
