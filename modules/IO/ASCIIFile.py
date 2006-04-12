#!/usr/bin/python
# M. Newville Univ of Chicago (2005)
#
# -------------
# Modifications
# -------------
#
# Should this get moved to pyMods?  
#
##########################################################################

import os
import types, re
import numpy as Num

class ASCIIFile:
    def __init__(self,fname=None,label=None,delim=None):
        self.titles  = []
        self.fname= None
        self.data    = {}
        self.delim   = r'\s+' # delimiter for labels...
        self.title_match = re.compile(r"[#a-d,h-z\{\}\[\]\(\):;\"'\?<>=_%]+",re.IGNORECASE)
        self.labline_match = re.compile(r"#?\-{3,}")
        self.commentchars  = re.compile(r"[#%;]+")
        if (fname):
            self.read(fname=fname,label=label)
        
    def write(self,fname=None,label=None):
        print 'ASCIIFile.write fname = ', fname
        
    def read(self,fname=None,label=None):
        if fname: self.fname =fname
        if self.fname == None :  return -1
        if not os.path.exists(self.fname): return -1
        f = open(fname)
        lines = f.readlines()
        f.close()
        lab_line = False
        tmp_dat = []
        tmp_label = ''
        for i in lines:
            i = i[:-1].strip()
            if self.title_match.match(i):
                if self.commentchars.match(i[0:1]):
                    i = i[1:].strip()
                if lab_line:
                    lab_line = False
                    tmp_label = i
                elif self.labline_match.match(i):
                    lab_line = True
                elif len(i)>1:
                    self.titles.append(i)
            else:
                words = i.split()
                all_numbers = True
                xw   = []
                for w in words:
                    try:
                        xw.append(float(w))
                    except:
                        all_numbers = False
                if all_numbers:
                    tmp_dat.append(xw)
                else:
                    if (self.commentchars.match(i[0:1])): i = i[1:].strip()
                    self.titles.append(i)
        #
        self.dat   = Num.transpose(Num.array(tmp_dat))
        self.ncols, self.npts = self.dat.shape
        print 'data shape: ', self.ncols, self.npts
        
        # column labels
        if label != None: tmp_label = label
        tmp_label = tmp_label.strip()
        try:
            lx = re.split(self.delim,tmp_label)
        except:
            lx = []
        # print 'label: ', len(lx), lx
        for i in range(len(lx),self.ncols+1):
            lx.append("col%i" % (i+1))
        self.labels = [l.lower() for l in lx[:self.ncols]]
        
    def get_titles(self):
        return self.titles
    
    def get_array(self,name):
        if type(name) == types.IntType:
            return self.dat[name]
        else:
            n = name.lower()
            if (n in self.labels):
                return self.dat[self.labels.index(n)]

    
if (__name__ == '__main__'):
    import sys
    u = ASCIIFile(sys.argv[1]) # ,label='c1 c2 c3')
    print u.get_titles()
    for i in u.labels:
        d = u.get_array(i)
        print i, ': ',  d[0:2]
