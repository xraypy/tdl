import types

class TNode:
    def __init__(self,catcode=None,exprstack=None,codeline=None):
        self.catcode   = catcode
        self.exprstack = exprstack
        self.codeline  = codeline
        self.arg   = []
        self.block = []
    def __repr__(self):
        return "<SNode: %s: '%s' (%s)>" % (self.catcode,self.exprstack,self.codeline)
    def look(self):
        print 'look: ', self.catcode
        return self

class Fifo:
    __slots__ = ('data',)
    def __init__(self,dat=None):
        self.data = [[], []]
        if dat:
            for i in dat: self.append(i)
    def append(self, value):
        self.data[1].append(value)
    def pop(self,n=-1):
        d = self.data
        if not d[0]:
            d.reverse()
            d[0].reverse()
        return d[0].pop()
    def __len__(self):
        return len(self.data[0]) + len(self.data[1])
   
    def __repr__(self):
        return "<Fifo len=%i id=0x%0x>" % (len(self),id(self))

class CodeLine:
    def __init__(self,text,fname=None,lineno=-1):
        self.text  = text
        self.fname = fname
        self.lineno= lineno
    def __repr__(self):
        return "<CodeLine %s (%s / %i)>" % (self.text,self.fname,self.lineno)
            
class Ev:
    EOF = '@EOF@'
    def __init__(self):
        self.x = 1
        self.text = Fifo()

    def load_file(self,fname):
        try:
            f = open(fname)
            l = f.readlines()
            f.close()
        except IOError:
            raise IOError, 'cannot load file %s '% fname
        self.load_statements(l,fname=fname)
        self.text.append(CodeLine(self.EOF,fname))

    def load_statements(self,t,fname='stdin'):
        " load a list of text lines to be parsed & executed"
        if type(t) == types.StringType:    t = [t]
        self.infile = fname
        for x,n in zip(t,range(len(t))):
            self.text.append(CodeLine(x[:-1],fname=fname,lineno=n))

def test():
    x = TNode(catcode='if',exprstack='code here',codeline=None)
    e  = Ev()
    e.load_file('test.tdl')
    # print e.text.data
    while True:
        try:
            u = e.text.pop()
            print u
        except IndexError:
            print 'done'
            break
