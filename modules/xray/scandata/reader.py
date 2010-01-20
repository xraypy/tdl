#######################################################################
"""
Tom Trainor (tptrainor@alaska.edu)
Class to read in generic 'scan data'

Modifications:
--------------

"""
#######################################################################
"""
Todo
- add escan files, ascii_scan files ....

"""
########################################################################

import types
from   glob import glob
import os, copy
import numpy as num

from   specfile import SpecFile
from   data import ScanData
import image_data
import xrf_data
import med_data

########################################################################
def spec_scan(spec,sc_num):
    """
    Return scan data instance from a specfile / scan number
    (specfile can be a specfile instance or string file name)
    """
    # Define these to make sure we get things sorted
    # correclty. These may be instrument/geometry specific!!
    PSIC_POSITIONER_KEYS = ['phi','chi','eta','mu','nu','del']

    # get the spec scan data
    if type(spec) == types.StringType:
        spec = SpecFile(spec)
        if spec == None: return None
    d = spec.scan_dict(sc_num)

    # parse positioner and scaler vals
    # note if a positioner or scaler was list in the data
    # array then we append the array to the positioners/scalers
    # otherwise the positioner value will be what was given
    # in d['P'] which should be a single value
    scalers = {}
    positioners = copy.copy(d['P'])
    for key in d['data'].keys():
        if key in positioners.keys():
            positioners[key] = num.array(d['data'][key])
        elif key in PSIC_POSITIONER_KEYS:
            positioners[key] = num.array(d['data'][key])
        else:
            scalers[key] = num.array(d['data'][key])
    name  = d['file'] + ' Scan ' + str(int(sc_num))
    dims  = d['nrow']
    paxis = d['labels'][0]
    pdet  = d['labels'][-1]

    # Grab state info.  Note make 'A' the essential
    # gonio angles at the start of the scan
    # this is APS sector 13 specific, and may change...
    # We need a switch here based on the identification
    # of the beamline and instrument...
    A = {}
    try:
        A['delta'] = d['P'].get('TwoTheta')
        A['eta']   = d['P'].get('theta')
        A['chi']   = d['P'].get('chi')
        A['phi']   = d['P'].get('phi')
        A['nu']    = d['P'].get('Nu')
        A['mu']    = d['P'].get('Psi')
        A['keta']  = d['P'].get('Omega')
        A['kap']   = d['P'].get('Kappa')
        A['kphi']  = d['P'].get('Phi')
    except:
        pass
    state = {'G':d.get('G'),'Q':d.get('Q'),'A':A,
             'ATTEN':d.get('ATTEN'), 'ENERGY':d.get('ENERGY')}

    # create scan data object
    data = ScanData(name=name,dims=[dims],scalers=scalers,
                    positioners=positioners,primary_axis=[paxis],
                    primary_det=pdet,state=state)
    return data

########################################################################
def med_scan(med,name='med'):
    """
    Return scan data instance from a list of med's

    Note makes a dummy axis:
       x = num.arange(len(med))
    
    """
    if type(med) != types.ListType: med = [med]
    npts = len(med)
    x = num.arange(float(npts))
    data = ScanData(name=fname,
                    dims = [npts],
                    positioners={'x':x},
                    primary_axis = 'x',
                    primary_det = 'med',
                    med=med)
    return data

########################################################################
def xrf_scan(xrf,name='xrf',lines=None):
    """
    Return scan data instance from a list of xrf's

    Note makes a dummy axis:
       x = num.arange(len(xrf))
    """
    if type(xrf) != types.ListType: xrf = [xrf]
    npts = len(xrf)
    x = num.arange(float(npts))
    data = ScanData(name=fname,
                    dims = [npts],
                    positioners={'x':x},
                    primary_axis = 'x',
                    primary_det = 'xrf',
                    xrf=xrf,
                    xrf_lines=lines)
    return data

########################################################################
def image_scan(image,name='image',rois=None):
    """
    Return scan data instance from a list of images

    Note makes a dummy axis:
       x = num.arange(len(xrf))

    """
    data = ScanData(name=fname,
                    dims = [len(image)],
                    primary_axis = 'image',
                    primary_det  = 'image',
                    image        = image,
                    image_rois   = rois)
    return data

########################################################################
class Reader:
    """
    Use the reader to read data that should all be treated the same.
    Therefore we assume that if we read spectra, they all have the same
    set of params (bad dets', taus, etc..).  Same with images.

    This class knows about the following data/file types:
     * spec files
     * escan (to be added)
     * med file (CARS and other fmts)
     * xrf data
     * images 
     * other (ie ssrl files, old spec and super etc..)

    Arguments for initialization:
     * spec = spec files, string or list of strings
     * spec_path = string path for spec file locations
     * escan_path = string path for escan file locations
     * spectra_path = string path for med or xrf file locations
     * image_path = string path for image file locations

    Examples:
     >>r = scandata.Reader(spec='spec_file.spc')
     >>s1 = r.spec_scan(1,image=True)

    """
    ########################################################################
    def __init__(self,spec=None,spec_path=None,escan_path=None,
                 spectra_path=None,image_path=None):
        # Spec
        self.spec_path       = spec_path
        self.spec_files      = []

        # escan
        self.escan_path       = escan_path
        self.escan_files      = []

        # med/xrf parameters
        self.med            = False
        self.xrf            = False
        self.spectra_path   = spectra_path
        self.med_params     = {'bad_mca_idx':[],'total':True,'align':True,
                               'correct':True,'tau':None,'det_idx':0,
                               'emin':-1.0,'emax':-1.0,
                               'fmt':'CARS','nfmt':3}
        self.xrf_params     = {}
        self.xrf_lines      = None

        # image parameters
        self.image          = False
        self.image_path     = image_path
        self.image_params   = {'type':'TIFF','nfmt':3}
        self.image_rois     = None

        # load spec file(s) if passed
        if spec: self.read_spec(spec)

    ########################################################################
    def __repr__(self):
        lout = "Spec Files:"
        for s in self.spec_files:
            lout = "%s\n  %s, first=%i, last=%i" % (lout,s.fname,s.min_scan,s.max_scan)
        lout = "%s\nRead med   = %s, Read xrf = %s" % (lout,self.med,self.xrf)
        lout = "%s\nRead image = %s"   % (lout,self.image)
        lout = "%s\nPaths:"            % (lout)
        lout = "%s\n  Spec Path  = %s" % (lout,self.spec_path)
        lout = "%s\n  Med Path   = %s" % (lout,self.spectra_path)
        lout = "%s\n  Image Path = %s" % (lout,self.image_path)
        lout = "%s\n  Escan Path = %s" % (lout,self.escan_path)
        if self.xrf_lines:
            lout = "%s\n  Xrf Lines = %s" % (lout,str(self.xrf_lines))
        if self.image_rois:
            lout = "%s\n  Image Rois = %s" % (lout,str(self.image_rois))
        return lout

    ########################################################################
    def med_scan(self,fname,start=-1,end=-1,path=None):
        """
        Read (collection) of med files into a scan data object

        Note makes a dummy axis:
           x = num.arange(len(med))
        """
        if path: self.spectra_path = path
        med = self._read_spectra(fname,start=start,end=end,xrf=False)
        if med == None: return None
        data = med_scan(med,name=fname)
        return data
        
    ########################################################################
    def xrf_scan(self,fname,start=-1,end=-1,path=None):
        """
        Read (collection) of xrf files into a scan data object

        Note makes a dummy axis:
           x = num.arange(len(xrf))
        """
        if path: self.spectra_path = path
        xrf = self._read_spectra(fname,start=start,end=end,xrf=True)
        if xrf == None: return None
        data = xrf_scan(xrf,name=fname)
        return data

    ########################################################################
    def image_scan(self,fname,start=-1,end=-1,path=None):
        """
        Read (collection) of image files into a scan data object

        Note makes a dummy axis:
           x = num.arange(len(image))

        """
        if path: self.image_path = path
        image = self._read_image(fname,start=start,end=end)
        if image == None: return None
        data = image_scan(image,name=fname,rois=self.image_rois)
        return data

    ########################################################################
    def escan(self,fname):
        pass
    
    ########################################################################
    def ascii_scan(self,fname):
        pass

    ########################################################################
    def read_spec(self,spec,path=None):
        """
        add a spec file (or files)
        """
        if path != None: self.spec_path = path
        if type(spec) == types.StringType: spec = [spec]
        
        for s in spec:
            sfile = self._spec(file=s)
            sfile.read()

    ########################################################################
    def list_spec(self, show=True):
        """
        return a list of spec file scans
        """
        sc_list = []
        for s in self.spec_files:
            sc_list.extend(s.list_scans())
        if show:
            for line in sc_list:
                print '* ' + line
            return
        else:
            return sc_list
    
    ########################################################################
    def spec_scan(self,scan,file=None,med=None,xrf=None,image=None):
        """
        get data from a spec scan
        """
        if med!=None: self.med = med
        if xrf!=None: self.xrf = xrf
        if image!=None: self.image = image
        # File
        spec = self._spec(file=file)
        if not spec: return None
        if scan == None: return None
        data = spec_scan(spec,scan)
        if not data: return None
        # Spectra
        if self.med or self.xrf:
            spectra_pfx = self._spec_spectra_path(spec,scan)
            start = 0
            end   = data.dims[0] - 1
        if self.med:
            med = self._read_spectra(spectra_pfx,start=start,end=end,xrf=False)
            if med == None:
                print "Warning, med files not read"
            else:
                med = med_data.MedScan(med)
                data.med = med
            if len(med) !=  data.dims[0]:
                print "Warning, number of spectra dont match scan"
        if self.xrf:
            if len(data.med)>0:
                xrf = xrf_data.med2xrf(xrf_params=self.xrf_params,
                                       lines = self.xrf_lines,
                                       det_idx=self.med_params['det_idx'],
                                       emin=self.med_params['emin'],
                                       emax=self.med_params['emax'])
                xrf = xrf_data.XrfScan(xrf)
                data.xrf = xrf
            else:
                xrf = self._read_spectra(spectra_pfx,start=start,end=end,xrf=True)
                if xrf == None:
                    print "Warning, xrf files not read"
                else:
                    xrf = xrf_data.XrfScan(xrf)
                    data.xrf = xrf
                if len(xrf) !=  data.dims[0]:
                    print "Warning, number of spectra dont match scan"
        # Images
        if self.image:
            image_pfx = self._spec_image_path(spec,scan)
            start = 0
            end   = data.dims[0] - 1
            image = self._read_image(image_pfx,start=start,end=end)
            data.image = image
        return data

    ########################################################################
    def _spec(self,file=None):
        """
        get spec file instance
        """
        if (len(self.spec_files) == 0) and (file == None):
            return None
        if file == None:
            return self.spec_files[0]
        else:
            for s in self.spec_files:
                if s.fname == file:
                    return s
        # If cant find it add new
        if self.spec_path != None:
            self.spec_path = os.path.normpath(self.spec_path)
            file = os.path.join(self.spec_path,file)
        tmp = SpecFile(file)
        if tmp._ok==True:
            self.spec_files.insert(0,tmp)
            return tmp
        else:
            return None

    ########################################################################
    def _spec_spectra_path(self,spec,scan):
        """
        Compute path for spectra associated with spec scan
        Note this is IDC 13 specific!

        if self.spectra_path == None:
            fpfx = spec_path/xrf_files/spec_pfx/nnn/spec_file.spc_nnn
        else:
            fpfx = spec_file.spc_nnn

        With nnn = scan number

        Note the read funciton will add '.point_number' to all
        the file names...

        """
        sc_num = '%03d' % int(scan)
        fpfx = "%s_%s" % (spec.fname, sc_num)
        if self.spectra_path == None:
            spec_pfx = spec.fname.rsplit('.',1)[0]
            path = os.path.join(spec.path,'xrf_files',spec_pfx,sc_num)
            fpfx = os.path.join(path,fpfx)
        return fpfx
    
    ########################################################################
    def _spec_image_path(self,spec,scan):
        """
        Compute path for spectra associated with spec scan
        Note this is IDC 13 specific!

        if self.spectra_path == None:
            fpfx = spec_path/images/spec_pfx/Snnn/spec_file.spc_Snnn
        else:
            fpfx = image_path/spec_file.spc_Snnn

        With nnn = scan number

        Note the read function will add '_point_num.tif' to all
        the file names...
        """
        sc_num = '%03d' % int(scan)
        fpfx = "%s_S%s" % (spec.fname, sc_num)
        if self.image_path == None:
            spec_pfx = spec.fname.rsplit('.',1)[0]
            path = os.path.join(spec.path,'images',spec_pfx,'S'+sc_num)
            fpfx = os.path.join(path,fpfx)
        return fpfx

    ########################################################################
    def _read_spectra(self,fname,start=-1,end=-1,xrf=True,scan_obj=True):
        """
        read spectra files and return a list of med or xrf objects
        """
        bad        = self.med_params['bad_mca_idx']
        total      = self.med_params['total']
        align      = self.med_params['align']
        correct    = self.med_params['correct']
        tau        = self.med_params['tau']
        det_idx    = self.med_params['det_idx']
        emin       = self.med_params['emin']
        emax       = self.med_params['emax']
        fmt        = self.med_params['fmt']
        nfmt       = self.med_params['nfmt']

        if self.spectra_path != None:
            fname = os.path.join(self.spectra_path,fname)
        
        if start > -1:
            if end == -1:
                ret = self._spectra_range(fname)  
                if ret:
                    (start,end) = ret
                else:
                    print "No files found"
                    return None
            if xrf:
                spectra = xrf_data.read_files(fname,start=start,end=end,nfmt=nfmt,
                                              bad_mca_idx=bad,total=total,align=align,
                                              correct=correct,tau=tau,det_idx=det_idx,
                                              emin=emin,emax=emax,fmt=fmt,
                                              xrf_params=self.xrf_params,
                                              lines=self.xrf_lines)
            else:
                spectra = med_data.read_files(fname,start=start,end=end,nfmt=nfmt,
                                              bad_mca_idx=bad,total=total,align=align,
                                              correct=correct,tau=tau,det_idx=det_idx,
                                              emin=emin,emax=emax,fmt=fmt)
        else:
            if xrf:
                spectra = xrf_data.read_file(fname,bad_mca_idx=bad,total=total,align=align,
                                             correct=correct,tau=tau,det_idx=det_idx,
                                             emin=emin,emax=emax,fmt=fmt,
                                             xrf_params=self.xrf_params,
                                             lines=self.xrf_lines)
            else:
                spectra = med_data.read_file(fname,bad_mca_idx=bad,total=total,align=align,
                                             correct=correct,tau=tau,det_idx=det_idx,
                                             emin=emin,emax=emax,fmt=fmt)
        if spectra == None:
            return []
        if type(spectra) != types.ListType:
            spectra = [spectra]
        return spectra

    ########################################################################
    def _spectra_range(self,fname):
        fname = fname + ".*"
        files = glob(fname)
        n = len(files)
        if n>0:
            st = files[0].split('.')[-1]
            en = files[n-1].split('.')[-1]
            try:
                st = int(st)
                en = int(en)
                return (st,en)
            except:
                return None
        else:
            return None
    
    ########################################################################
    def _read_image(self,fname,start=-1,end=-1):
        """
        read image files
        """
        nfmt = self.image_params['nfmt']
        
        if self.image_path != None:
            fname = os.path.join(self.image_path,fname)
        
        if start > -1:
            if end == -1:
                ret=self._image_range(fname)
                if ret:
                    (start,end) = ret
                else:
                    print "No files found"
                    return None
            #print "%s, start=%i, end = %i, nfmt = %i" % (fname,start,end,nfmt)
            image = image_data.read_files(fname,start=start,end=end,nfmt=nfmt)
            image = image_data.read_files(fname,start=start,end=end,nfmt=nfmt)
        else:
            image = image_data.read_file(fname)

        return image

    ########################################################################
    def _image_range(self,fname):
        fname = fname + "_*.tif"
        files = glob(fname)
        n = len(files)
        if n>0:
            st = files[0].split('_')[-1]
            st = st.split('.')[0]
            en = files[n-1].split('_')[-1]
            en = en.split('.')[0]
            #try:
            st = int(st)
            en = int(en)
            return (st,en)
            #except:
            #    return None
        else:
            return None
    
############################################################################
############################################################################
if __name__ == '__main__':
    pass
