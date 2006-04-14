#!/usr/bin/env python
##
## IO functionality

import os
import sys
import types
import FileIO.ASCIIFile as ASCIIFile

##
def tdl_read_ascii(fname, group=None, tdl=None,debug=False, **kw):
    " read ascii file of tdl code"
    if tdl == None:
        if debug: print 'cannot read file %s ' % fname
        return None
    if debug: print 'reading.... ', fname
    if not os.path.exists(fname):
        print 'read_ascii: cannot find file %s ' % fname
        return None        
    try:
        f = ASCIIFile(fname)
    except:
        print 'read_ascii: error reading file %s ' % fname        
        return None
    # save current group name
    savegroup = tdl.symbolTable.getDataGroup()
    if group == None:
        group = savegrup
    else:
        group = tdl.symbolTable.setDataGroup(group)
        
    tdl.symbolTable.setVariable('titles', f.get_titles(), group=group)
    tdl.symbolTable.setVariable('column_labels', f.labels, group=group)        
    ret = [group]
    for i in f.labels:
        ret.append(i)
        tdl.symbolTable.setVariable(i, f.get_array(i), group=group)
    if debug: print 'read done.'
    # return default group to original
    tdl.symbolTable.setDataGroup(savegroup)
    return ret
#
def tdl_write_ascii(fname,  *arr,**kw):
    " write ascii file of tdl code -- could be a lot better! "
    print 'write ascii ', fname
    tdl = kw['tdl']
    if tdl is None: return None
    arr_out=[]
    sca_val = []
    sca_nam = []
    labels =''
    group = tdl.symbolTable.getDataGroup()
    npts = []
    for i in arr:
        o = tdl.symbolTable.getVariable(i).value
        try:
            n = len(o)
        except TypeError:
            n = 0
        if n<=1:
            sca_nam.append(i)
            sca_val.append(o)
        else:
            labels = "%s %s" % (labels,i)
            arr_out.append(list(o))
            npts.append(len(o))
                 
    f = open(fname,'w')
    f.write("# file written by tdl write_ascii() \n")
    for n,v in zip(sca_nam,sca_val):
        f.write("# %s = %s \n" %  (n,v))
    f.write("-----------------------------------------\n")
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
# todo:
#  -  control of legend, titles, etc  (tdl functions??)
#  -  draw box on zooming
#  -  GUI config

title = 'File IO routine'

HelpIO = """
  File Handling in tdl:
   
"""

_help_ = {'file_io': HelpIO}
_func_ = {'_builtin.read_ascii':(tdl_read_ascii, None),
          '_builtin.write_ascii':(tdl_write_ascii, None)}
