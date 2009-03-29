##########################################################################
"""
M. Newville

- module is out of date

Modifications:
--------------


"""
##########################################################################

import os
import sys
import types

from util import datalen, verify_tdl
import ASCIIFile.ASCIIFile as ASCIIFile

##########################################################################
def tdl_read_ascii(fname, tdl=None,debug=False, **kw):
    " read ascii file of tdl code"
    verify_tdl(tdl,'read_ascii')

    if debug: print 'reading.... ', fname
    if not os.path.exists(fname):
        print 'read_ascii: cannot find file %s ' % fname
        return None        
    try:
        f = ASCIIFile(fname)
    except:
        print 'read_ascii: error reading file %s ' % fname        
        return None

    group = tdl.symbolTable.addTempGroup(toplevel=True,**kw)
    gname = group.name
    def setvar(name,value):
        tdl.symbolTable.setVariable("%s.%s" % (gname,name), value)
        
    setvar('titles', f.get_titles())
    setvar('column_labels', f.labels)

    for i in f.labels: setvar(i, f.get_array(i))
    if debug: print 'read done.'
    return group
#
def tdl_write_ascii(fname,  *arr,**kw):
    " write ascii file of tdl code -- could be a lot better! "
    print 'write ascii ', fname, kw
    tdl    = kw['tdl']
    labels = kw.get('label', '')
    if tdl is None: return None
    arr_out=[]
    sca_val = []
    sca_nam = []
    group = tdl.symbolTable.getDataGroup()
    npts = []
    for item in arr:
        n     = datalen(item)
        dtype = type(item)
        itype = 'literal'
        iname = ' '        
        ival  = item
        if dtype is types.StringType:
            if True:# try:
                sym   = tdl.symbolTable.getVariable(item)
                dtype = sym.type
                ival  = sym.value
                n     = datalen(sym.value)
                iname = sym.name
                itype = 'symbol'
            else: # except:
                itype = 'unknown'
                print ' cannot handle ', item
        if itype == 'unknown': continue
        print '>> ', n, itype, dtype, iname
        if n == 1 or dtype == 'string':
            sca_nam.append(iname)
            sca_val.append(ival)            
            print '-> scalar ', iname, ival
        else:
            labels = "%s %s" % (labels,iname)
            arr_out.append(list(ival))
            npts.append(len(ival))
            
                 
    f = open(fname,'w')
    f.write("# file written by tdl write_ascii() \n")
    for n,v in zip(sca_nam,sca_val):
        f.write("# %s = %s \n" %  (n,v))
    f.write("#-----------------------------------------\n")
    f.write("# %s\n" % labels)
    nout = len(npts)
    mout = npts[0]
    for i in npts:
        if i > mout: mout = i
        
    for i in range(mout):
        t = ""
        for s in range(nout):
            t = "%s    %g " % (t,arr_out[s][i])
        f.write("%s\n" % t)
        
    f.close()

##########################
title = 'File IO routine'

HelpIO = """
  File Handling in tdl:
   
"""

_help_ = {'file_io': HelpIO}
_func_ = {'_builtin.read_ascii':(tdl_read_ascii, None),
          '_builtin.write_ascii':(tdl_write_ascii, None)}
