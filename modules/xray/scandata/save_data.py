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
root/med
root/xrf
etc

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
def get_file(fname,path=None):
    """
    Open / create file
    """
    try:
        if path != None:
            fname = os.path.join(path,fname)
            fname = os.path.abspath(fname)
        else:
            fname = os.path.abspath(fname)
        if os.path.exists(fname):
            h = tables.openFile(fname,mode="a")
        else:
            h = tables.openFile(fname,mode="w",title="Scan Data Archive")
        return h
    except:
        _cleanup()
        print "Unable to get hdf file %s", fname

def list_file(fname,path=None,display=True):
    h = get_file(fname,path=path)
    if h == None: return
    ll = []
    for a in h.walkNodes('/'):
        if display:
            print a
        else:
            ll.append(a)
    h.close()
    if display: return
    return ll

################################################################################
def _cleanup():
    tables.file.close_open_files()

################################################################################
def write_scandata(fname,data,setname=None,path=None,overwrite=True):
    """
    Write scan data
    """
    # get data as a dictionary of arrays.
    (d,image,med,xrf) = _scan_to_dict(data)
    h = get_file(fname,path)
    if h == None: return
    #try:
    _write_scan(h,d,setname,overwrite=overwrite)
    h.close()
    #except:
    #    _cleanup()
    #    print "Unable to write scandata"
    return

def _write_scan(h,data,setname,overwrite=True):
    """
    note create array stores any homogeneous data type,
    and remembers what it is so you get back correct
    type when you read it back in (except tuples->lists)
    """
    if not hasattr(h.root,'scans'):
        h.createGroup(h.root,'scans',"Scan Data")
    if setname == None:
        names = h.root.scans._v_children.keys()
        setname = calc_next_setname(names)
    if hasattr(h.root.scans,setname):
        print "Archive File %s: Setname %s already exists" % (h.filename,setname)
        if overwrite==False:
            print "Data is not overwritten"
            return
        else:
            print "Data being overwritten"
            h.removeNode(h.root.scans,setname,recursive=True)
            h.createGroup(h.root.scans,setname,"Scan Data")
    else:
        h.createGroup(h.root.scans,setname,"Scan Data")

    # add data
    grp0 = '/scans/' + setname
    h.createArray(grp0,'name',data['name'],data['name'])
    h.createArray(grp0,'dims',data['dims'],'dims')
    h.createArray(grp0,'primary_axis',data['primary_axis'],'primary_axis')
    h.createArray(grp0,'primary_det',data['primary_det'],'primary_det')
    # scalars
    h.createGroup(grp0,'scalar',"Scan Data Scalars")
    grp = grp0 + '/scalar'
    for j in range(len(data['scnames'])):
        name = data['scnames'][j]
        val = data['scvals'][j]
        h.createArray(grp,name,val,name)
    # positioners
    h.createGroup(grp0,'positioners',"Scan Data Positioners")
    grp = grp0 + '/positioners'
    for j in range(len(data['ponames'])):
        name = data['ponames'][j]
        val = data['povals'][j]
        h.createArray(grp,name,val,name)
    # state
    h.createGroup(grp0,'state',"Scan Data State Variables")
    grp = grp0 + '/state'
    for j in range(len(data['stnames'])):
        name = data['stnames'][j]
        val = data['stvals'][j]
        try:
            h.createArray(grp,name,val,name)
        except:
            print "Unable to write state variable: ", name, val
    
    # images (should add only if data has images..)
    if not hasattr(h.root,'images'):
        h.createGroup(h.root,'images',"Image Data")
    # xrf/med (should add only if data has med/xrf..)
    if not hasattr(h.root,'med'):
        h.createGroup(h.root,'med',"MED Data")
    if not hasattr(h.root,'xrf'):
        h.createGroup(h.root,'xrf',"XRF Data")

################################################################################
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
def _scan_to_dict(data):
    """
    Turn a scan data object into dictionary of arrays
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
        image = []
        for j in range(len(data.image.image)):
            # tedious
            #d['imparams']['aparam'][j] = d.image.aparam[j]
            #s[0].image._is_integrated    
            #s[0].image.bgrpar
            #s[0].image.im_max
            #s[0].image.peaks             
            #s[0].image.rois
            #s[0].image.rotangle          
            image.append(data.image.image[j])
    else:
        d['imparams'] = None
        image = None
    # need same for med and xrf
    med = None
    xrf = None
    return (d,image,med,xrf)

################################################################################
def read_scandata(file,setname,path=None):
    """
    Read data
    """
    h = get_file(file,path)
    if h==None: return
    data = {}
    try:
        grp = '/scans/' + setname
        #items = getattr(h,grp)
        items  = h.getNode('/scans',setname)
        names = items._v_children.keys()
    except:
        return data
    for n in names:
        node  = h.getNode(grp,n)
        data[n] = node.read()
    h.close()
    return data

################################################################################
################################################################################
if __name__ == '__main__':
    x = num.arange(100.)
    y = num.sin(x/num.pi)
    d = {'x':x,'y':y,'t':'a string','i':1,'p':(1,2),'Q':[1,2,3]}
    write_scandata('xx.hdf',d,setname='S1')
    zz = read_scandata('xx.hdf','S1')
    print 'zz', zz.keys()
    #
    list_file('xx.hdf')
    
