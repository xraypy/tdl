"""
Saving and restoring data objects

Authors/Modifications:
----------------------
* Tom Trainor (tptrainor@alaska.edu)

Notes:
------
This module will allow writing scan data to
an hdf file, and read it back in.  This makes
the data portable ie readable outside python,
and provides a means to use disk files for
caching certain data types like images so they dont
have to be stored in memory.  

The hdf files are structured as follows:
root/scan_data/...
root/image_data/...
root/med_data/...  
root/xrf_data/...  
root/ctr_data/...  

Examples:
--------
>>write_ctrdata(ctr,file='ctr.h5')
>>ctr = read_ctrdata(file='ctr.h5',name='ctr')

Todo:
-----
* work in progress, have mostly completed writing ScanData, ImageData
  and CtrData to file.  
* have to do writing of med and xrf
* have to make read methods that will read the files and return the
  appropriate class instances.  
"""
################################################################################

import types
import os
import copy
import numpy as num

try:
    import tables
except:
    print "    ** Warning tables module (hdf interface) cannot be loaded in 'save_data'"
    
import scan_data
import image_data
import ctr_data

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
    """
    list file contents
    """
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
    """
    calculate the next set name
    """
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
    """
    make sure any open files are closed
    """
    import sys
    xx, yy, zz = sys.exc_info()
    sys.excepthook(xx,yy,zz)
    #
    tables.file.close_open_files()

################################################################################
def write_ctrdata(fname,ctr,setname='ctr',path=None,overwrite=True):
    """
    Write scan data
    """
    h = get_file(fname,path)
    if h == None: return

    # write all scans
    try:
        for j in range(len(ctr.scan)):
            d = _scan_data(ctr.scan[j])
            scanname = "%s_S%03d" % (setname,j)
            _write_scan(h,d,scanname,overwrite=overwrite)
            if hasattr(ctr.scan[j],'image'):
                im = _image_data(ctr.scan[j].image)
                _write_image(h,im,scanname,overwrite=overwrite)
    except:
        _cleanup()
        print "Unable to write scandata"
    # now write ctr data
    #try:
    d = _ctr_data(ctr)
    _write_ctr(h,d,setname,overwrite=overwrite)
    #except:
    #    _cleanup()
    #    print "Unable to write scandata"
    h.close()
    return d

################################################################################
def _write_ctr(h,data,setname='ctr',overwrite=True):
    """
    note create array stores any homogeneous data type,
    and remembers what it is so you get back correct
    type when you read it back in (except tuples->lists)
    """
    if not hasattr(h.root,'ctr_data'):
        h.createGroup(h.root,'ctr_data',"Ctr Data")
    if hasattr(h.root.ctr_data,setname):
        print "Archive File %s: Setname %s already exists" % (h.filename,setname)
        if overwrite==False:
            print "Data is not overwritten"
            return
        else:
            print "Data being overwritten"
            h.removeNode(h.root.ctr_data,setname,recursive=True)
            h.createGroup(h.root.ctr_data,setname,"Ctr Data")
    else:
        h.createGroup(h.root.ctr_data,setname,"Ctr Data")
    # add data
    """
    class CtrTable(table.IsDescription):
        I_lbl  = StringCol(16)
        In_lbl = StringCol(16)
        Ie_lbl = StringCol(16)
        Ib_lbl = StringCol(16)
        stype  = StringCol(16)
    """
    grp0 = '/ctr_data/' + setname
    if len(data['bad'])>0:
        h.createArray(grp0,'bad',data['bad'],'bad points')
    h.createArray(grp0,'scan_index',data['scan_index'],'scan_index')
    h.createArray(grp0,'scan_type',data['scan_type'],'scan_type')
    h.createArray(grp0,'Ilbl',data['I_lbl'],'I_lbl')
    h.createArray(grp0,'Inorm_lbl',data['In_lbl'],'Inorm_lbl')
    h.createArray(grp0,'Ierr_lbl',data['Ie_lbl'],'Ierr_lbl')
    h.createArray(grp0,'Ibgr_lbl',data['Ib_lbl'],'Ibgr_lbl')
    #
    h.createArray(grp0,'H',data['H'],'H')
    h.createArray(grp0,'K',data['K'],'K')
    h.createArray(grp0,'L',data['L'],'L')
    h.createArray(grp0,'I',data['I'],'I')
    h.createArray(grp0,'Inorm',data['In'],'Inorm')
    h.createArray(grp0,'Ierr',data['Ie'],'Ierr')
    h.createArray(grp0,'Ibgr',data['Ib'],'Ibgr')
    h.createArray(grp0,'ctot',data['c'],'ctot')
    h.createArray(grp0,'F',data['F'],'F')
    h.createArray(grp0,'Ferr',data['Fe'],'Ferr')
    h.createArray(grp0,'beam_slit_horz',data['beam_slit_horz'],'beam_slit_horz') 
    h.createArray(grp0,'beam_slit_vert',data['beam_slit_vert'],'beam_slit_vert')
    h.createArray(grp0,'det_slit_horz',data['det_slit_horz'],'det_slit_horz') 
    h.createArray(grp0,'det_slit_vert',data['det_slit_vert'],'det_slit_vert')
    h.createArray(grp0,'geom',data['geom'],'geom')
    h.createArray(grp0,'scale',data['scale'],'scale')
    h.createArray(grp0,'sdia',data['sdia'],'sdia')
    h.createArray(grp0,'spoly',data['spoly'],'spoly')
    if type(data['sangle']) == types.DictType:
        for name,val in data['sangle'].items():
            name = 'sangle_'+name
            h.createArray(grp0,name,val,name)
    else:
        h.createArray(grp0,'sangle',data['sangle'],'sangle')

################################################################################
def _ctr_data(ctr):
    """
    Turn a ctr data object into dictionary of arrays
    """
    if not isinstance(ctr,ctr_data.CtrData):
        print "Warning data is not a CtrData instance"
    npts = len(ctr.L)
    d = {}
    d['bad']        = ctr.bad
    d['scan_index'] = ctr.scan_index
    d['scan_type']  = ctr.scan_type
    d['I_lbl']      = ctr.labels['I']
    d['In_lbl']     = ctr.labels['Inorm']
    d['Ie_lbl']     = ctr.labels['Ierr']
    d['Ib_lbl']     = ctr.labels['Ibgr']
    d['H']  = ctr.H   
    d['K']  = ctr.K     
    d['L']  = ctr.L    
    d['I']  = ctr.I     
    d['In'] = ctr.Inorm 
    d['Ie'] = ctr.Ierr  
    d['Ib'] = ctr.Ibgr  
    d['c']  = ctr.ctot  
    d['F']  = ctr.F     
    d['Fe'] = ctr.Ferr
    #
    d['beam_slit_horz'] = []
    d['beam_slit_vert'] = []
    d['det_slit_horz'] = []
    d['det_slit_vert'] = []
    d['geom'] = []
    d['scale'] = []
    d['sdia'] = []
    d['spoly'] = []
    d['sangle'] = []
    for j in range(npts):
        d['geom'].append(ctr.corr_params[j].get('geom',''))
        d['scale'].append(ctr.corr_params[j].get('scale',1.0))
        if ctr.corr_params[j]['beam_slits'] != None:
            d['beam_slit_horz'].append(ctr.corr_params[j]['beam_slits'].get('horz',0.))
            d['beam_slit_vert'].append(ctr.corr_params[j]['beam_slits'].get('vert',0.))
        else:
            d['beam_slit_horz'].append(0.)
            d['beam_slit_vert'].append(0.)
        if ctr.corr_params[j]['det_slits'] != None:
            d['det_slit_horz'].append(ctr.corr_params[j]['det_slits'].get('horz',0.))
            d['det_slit_vert'].append(ctr.corr_params[j]['det_slits'].get('vert',0.))
        else:
            d['det_slit_horz'].append(0.)
            d['det_slit_vert'].append(0.)
        d['sdia'].append(ctr.corr_params[j]['sample'].get('dia',0.0))
        d['spoly'].append(ctr.corr_params[j]['sample'].get('polygon',0.0))
        d['sangle'].append(ctr.corr_params[j]['sample'].get('angles',0.0))
    if type(d['sangle'][0]) == types.DictType:
        tmp = {}
        names = d['sangle'][0].keys()
        for key in names:
            tmp[key] = []
        for sa in d['sangle']:
            for key in names:
                tmp[key].append(sa.get(key,0.))
        d['sangle'] = tmp
    return d

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
        names = h.root.scan_data._v_children.keys()
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
    """
    write image data
    """
    if not hasattr(h.root,'image_data'):
        h.createGroup(h.root,'image_data',"Image Data")
    if hasattr(h.root.image_data,setname):
        print "Image Archive File %s: Setname %s already exists" % (h.filename,setname)
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
    """
    write med
    """
    if not hasattr(h.root,'med'):
        h.createGroup(h.root,'med',"MED Data")

################################################################################
def _write_xrf(h,data,setname,overwrite=True):
    """
    write xrf
    """
    if not hasattr(h.root,'xrf'):
        h.createGroup(h.root,'xrf',"XRF Data")

################################################################################
def _scan_data(data):
    """
    Turn a scan data object into dictionary of arrays
    """
    if not isinstance(data,scan_data.ScanData):
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
    """
    med data
    """
    # need same for med and xrf
    med = None
    return med

################################################################################
def _xrf_data(data):
    """
    xrf data
    """
    xrf = None
    return xrf

################################################################################
def read_scandata(file,setname,path=None):
    """
    Read data
    """
    h = get_file(file,path)
    if h==None: return
    data = _read_scan(h,setname)
    h.close()
    return data

################################################################################
def _read_scan(h,setname):
    """
    read
    """
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
    return data


################################################################################
################################################################################
if __name__ == '__main__':
    """
    x = num.arange(100.)
    y = num.sin(x/num.pi)
    sc = {'x':x,'t':'a string','i':1}
    po = {'y':y,'p':(1,2)}
    st = {'Q':[1,2,3]}
    d = scandata.ScanData(name='x',dims=[1,],scalers=sc,positioners=po,
                          primary_axis=['x'],primary_det=['y'],state=st)
    write_scandata('xx.h5',d,setname='S1')
    zz = read_scandata('xx.h5','S1')
    print 'zz', zz.keys()
    #
    list_file('xx.h5')
    """
    pass
