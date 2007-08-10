from Parser import ExprParser
p = ExprParser()

t = """
a[0]                      
a['b']

"""

exprs = t.split('\n')

for i in exprs:
    i = i.strip()
    if len(i)>2:
        x = p.compile(i)
        # n = i + 30*' '
        print " %s => " % (i),
        for l in x[:-1]: print '%s, '% l,
        print '%s'% x[-1]        
