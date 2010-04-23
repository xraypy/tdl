"""
Saving and restoring data objects

Authors/Modifications:
----------------------
Tom Trainor (tptrainor@alaska.edu)

Notes:
------

Structure hdf files as follows:
root/scan_data/...
root/image_data/...
root/med_data/...
root/xrf_data/...
root/ctr_data/...


Examples:
--------
>>scandata.save(ctr,file='ctr.h5')
>>ctr = scandata.restore(file='ctr.h5',name='ctr')


Todo:
-----
- work in progress
- have nto done med or xrf yet
"""
#######################################################################

import types
import os
import copy
import numpy as num
import tables

import data as scandata
import image_data

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
        #print fname
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
def _cleanup():
    import sys
    xx, yy, zz = sys.exc_info()
    sys.excepthook(xx,yy,zz)
    #
    tables.file.close_open_files()

################################################################################
def write_scandata(fname,data,setname=None,path=None,overwrite=True):
    """
    Write scan data
    """
    h = get_file(fname,path)
    if h == None: return

    # Scan data
    try:
        d = _scan_data(data)
        _write_scan(h,d,setname,overwrite=overwrite)
    except:
        _cleanup()
        print "Unable to write scandata"
    
    # Image data
    if hasattr(data,'image'):
        try:
            im = _image_data(data.image)
            _write_image(h,im,setname,overwrite=overwrite)
        except:
            _cleanup()
            print "Unable to write imagedata"

    # Med data
    if hasattr(data,'med'):
        pass
    
    # Xrf data
    if hasattr(data,'xrf'):
        pass
    
    h.close()
    return

################################################################################
def _write_scan(h,data,setname,overwrite=True):
    """
    note create array stores any homogeneous data type,
    and remembers what it is so you get back correct
    type when you read it back in (except tuples->lists)
    """
    if not hasattr(h.root,'scan_data'):
        h.createGroup(h.root,'scan_data',"Scan Data")
    if setname == None:
        names = h.root.scans._v_children.keys()
        setname = calc_next_setname(names)
    if hasattr(h.root.scan_data,setname):
        print "Archive File %s: Setname %s already exists" % (h.filename,setname)
        if overwrite==False:
            print "Data is not overwritten"
            return
        else:
            print "Data being overwritten"
            h.removeNode(h.root.scan_data,setname,recursive=True)
            h.createGroup(h.root.scan_data,setname,"Scan Data")
    else:
        h.createGroup(h.root.scan_data,setname,"Scan Data")

    # add data
    grp0 = '/scan_data/' + setname
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
            print "Unable to write state variable: ", name

################################################################################
def _write_image(h,data,setname,overwrite=True):
    if not hasattr(h.root,'image_data'):
        h.createGroup(h.root,'image_data',"Image Data")

    if hasattr(h.root.image_data,setname):
        print "Image Archive File %s: Setname %s already exists" % (fname,setname)
        if overwrite==False:
            print "Data is not overwritten"
            return
        else:
            print "Data being overwritten"
            h.removeNode(h.root.image_data,setname,recursive=True)
            h.createGroup(h.root.image_data,setname,"Scan Data")
    else:
        h.createGroup(h.root.image_data,setname,"Scan Image Data")
    # images
    grp = '/image_data/' + setname
    h.createArray(grp,'images',data['images'],'Images')
    # image params
    for name,val in data['imparams'].items():
        h.createArray(grp,name,val,name)
    # bgr params
    for name,val in data['bparams'].items():
        h.createArray(grp,name,val,name)

################################################################################
def _write_med(h,data,setname,overwrite=True):
    if not hasattr(h.root,'med'):
        h.createGroup(h.root,'med',"MED Data")

################################################################################
def _write_xrf(h,data,setname,overwrite=True):
    if not hasattr(h.root,'xrf'):
        h.createGroup(h.root,'xrf',"XRF Data")

################################################################################
def _scan_data(data):
    """
    Turn a scan data object into dictionary of arrays
    """
    if not isinstance(data,scandata.ScanData):
        print "Warning data is not a ScanData instance"
    
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
    return d

################################################################################
def _image_data(imdata):
    """
    This parses image data
    """
    if not isinstance(imdata,image_data.ImageScan):
        print "Warning data is not a ImageScan instance"

    imparams = {}
    imparams['rois'] = []
    imparams['rotangle'] = []
    imparams['im_max'] = []
    imparams['I'] = []
    imparams['Ierr'] = []
    imparams['Ibgr'] = []
    imparams['I_c'] = []
    imparams['Ierr_c'] = []
    imparams['Ibgr_c'] = []
    imparams['I_r'] = []
    imparams['Ierr_r'] = []
    imparams['Ibgr_r'] = []
    
    bparams = {}
    bkey = imdata.bgrpar[0].keys()
    for key in bkey:
        bparams[key] = []

    images = []
    npts = len(imdata.image)

    for j in range(npts):
        images.append(imdata.image[j])
        imparams['rois'].append(imdata.rois[j])
        imparams['rotangle'].append(imdata.rotangle[j])
        imparams['im_max'].append(imdata.im_max[j])
        imparams['I'].append(imdata.peaks['I'][j])
        imparams['Ierr'].append(imdata.peaks['Ierr'][j])
        imparams['Ibgr'].append(imdata.peaks['Ibgr'][j])
        imparams['I_c'].append(imdata.peaks['I_c'][j])
        imparams['Ierr_c'].append(imdata.peaks['Ierr_c'][j])
        imparams['Ibgr_c'].append(imdata.peaks['Ibgr_c'][j])
        imparams['I_r'].append(imdata.peaks['I_r'][j])
        imparams['Ierr_r'].append(imdata.peaks['Ierr_r'][j])
        imparams['Ibgr_r'].append(imdata.peaks['Ibgr_r'][j])
        for key in bkey:
            bparams[key].append(imdata.bgrpar[j][key])

    data = {'images':images,'imparams':imparams,'bparams':bparams}
    return data

################################################################################
def _med_data(data):
    # need same for med and xrf
    med = None
    return med

################################################################################
def _xrf_data(data):
    xrf = None
    return xrf


################################################################################
def write_ctrdata(fname,data,setname=None,path=None,overwrite=True):
    pass

################################################################################
def read_scandata(file,setname,path=None):
    """
    Read data
    """
    h = get_file(file,path)
    if h==None: return
    data = {}
    try:
        grp = '/scan_data/' + setname
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
    sc = {'x':x,'t':'a string','i':1}
    po = {'y':y,'p':(1,2)}
    st = {'Q':[1,2,3]}
    d = scandata.ScanData(name='x',dims=[1,],scalers=sc,positioners=po,
                          primary_axis=['x'],primary_det=['y'],state=st)
    write_scandata('xx.hdf',d,setname='S1')
    zz = read_scandata('xx.hdf','S1')
    print 'zz', zz.keys()
    #
    list_file('xx.hdf')
    