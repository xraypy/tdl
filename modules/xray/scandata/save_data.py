"""
Saving and restoring data objects

Authors/Modifications:
----------------------
Tom Trainor (tptrainor@alaska.edu)

Notes:
------

Structure hdf files as follows:

root/scan
root/image
root/ctr
etc...


Examples:
--------

>>scandata.save(ctr,file='ctr.h5')
>>ctr = scandata.restore(file='ctr.h5',name='ctr')


Todo:
-----
- work in progress

"""
#######################################################################

import types
import os
import copy
import numpy as num
import tables

################################################################################
def get_file(file,path=None):
    """
    Open / create file
    """
    if path != None:
        fname = os.path.join(path,file)
        fname = os.path.abspath(fname)
    else:
        fname = os.path.abspath(file)
    if os.path.exists(fname):
        h = tables.openFile(fname,mode="a")
    else:
        h = tables.openFile(fname,mode="w",title="Scan Data Archive")
    return h

def calc_next_setname(names):
    if len(names) == 0:
        return 'S000'
    idx = []
    try:
        for n in names:
            idx.append(int(n[1:]))
        m = num.max(idx)
        s = "S%03d" % (m)
        return s
    except:
        return 'S000'
    
################################################################################
def write_scandata(file,data,path=None,setname=None):
    """
    """
    if type(data) != types.ListType:
        data = [data]
    h = get_file(file,path)

    if not hasattr(h.root,'scans'):
        h.createGroup(h.root,'scans',"Scan Data")
    if setname == None:
        names = h.root.scans._v_children.keys()
        setname = calc_next_setname(names)
    if hasattr(h.root.scans,setname):
        print "   Archive File %s: Setname %s already exists" % (fname,setname)
    else:
        h.createArray('/scans',setname,data,'Scan Data')

    # images
    if not hasattr(h.root,'images'):
        h.createGroup(h.root,'images',"Image Data")
    #
    h.close()

def _scan_data(data):
    """
    turn scan data object into dictionary of arrays
    """
    d = {}
    d['name'] = data.name
    d['dims'] = data.dims
    d['primary_axis'] = data.primary_axis
    d['primary_det'] = data.primary_det
    # scalers
    d['scnames'] = []
    d['scvals']  = []
    for (n,v) in data.scalers.items():
        d['scnames'].append(n)
        d['scvals'].append(v)
    # positioners
    d['ponames'] = []
    d['povals']  = []
    for (n,v) in data.positioners.items():
        d['ponames'].append(n)
        d['povals'].append(v)
    # state
    d['stnames'] = []
    d['stvals']  = []
    for (n,v) in data.state.items():
        d['stnames'].append(n)
        d['stvals'].append(v)
    # images
    if hasattr(data,'image'):
        d['imparams'] = {}
        imdata = []
        for j in range(len(d.image.image)):
            # tedious
            #d['imparams']['aparam'][j] = d.image.aparam[j]
            #s[0].image._is_integrated    
            #s[0].image.bgrpar
            #s[0].image.im_max
            #s[0].image.peaks             
            #s[0].image.rois
            #s[0].image.rotangle          
            imdata.append(d.image.image[j])
    else:
        d['imparams'] = None
        imdata = None
    # need same for med and xrf
    return (data,image,med,xrf)

################################################################################
def read_scandata(path,file,name):
    h = get_file(file,path)
    n  = h.getNode('/scans',name)
    data = n.read()
    h.close()
    return data

################################################################################
################################################################################
if __name__ == '__main__':
    pass
