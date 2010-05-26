"""
Read 'scan data'

Authors/Modifications:
----------------------
* Tom Trainor (tptrainor@alaska.edu)


Todo:
-----
* add escan files, ascii_scan files ....
* use the 'geo' parameter to allow spec
  files from vairous beamlines / various
  gonio geometries to be used
"""
#######################################################################

import types
from   glob import glob
import os, copy
import numpy as num

from   specfile import SpecFile
from   scan_data import ScanData
from   detector import med
from   xrf import xrf_model
import image_data
import xrf_data
import med_data

########################################################################
def spec_scan(spec,sc_num,geo='PSIC_APS_S13'):
    """
    Return a ScanData instance from a specfile / scan number
    
    * spec can be a specfile instance or string file name
    * sc_num is the scan number
    * geo is a geometry label that helps the reader parse gonio
      angles depending on the particular beamline / geometry
      used for data collection.  

    Notes:
    ------
    The vector A is appended to the state info.  This
    is the set of essential gonio angles at the start of the scan
    These are defined based on the 'geo' label and are geometry/
    beamline specific

    Example:
    --------
    >>s = spec_scan('datafile.spc',12)
    >>plot(s['phi'],s['bicron'])
    
    """
    # Define to make sure we get things sorted correclty. 
    if geo=='PSIC_APS_S13':
        POSITIONER_KEYS = ['phi','chi','eta','mu','nu','del']
    else:
        POSITIONER_KEYS =[]
    
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
        elif key in POSITIONER_KEYS:
            positioners[key] = num.array(d['data'][key])
        else:
            scalers[key] = num.array(d['data'][key])
    name  = d['file'] + ' Scan ' + str(int(sc_num))
    dims  = d['nrow']
    paxis = d['labels'][0]
    pdet  = d['labels'][-1]

    # Grab state info.  
    try:
        A = []
        if geo=='PSIC_APS_S13':
            A.append(d['P'].get('TwoTheta'))
            A.append(d['P'].get('theta'))
            A.append(d['P'].get('chi'))
            A.append(d['P'].get('phi'))
            A.append(d['P'].get('Nu'))
            A.append(d['P'].get('Psi'))
            A.append(d['P'].get('Omega'))
            A.append(d['P'].get('Kappa'))
            A.append(d['P'].get('Phi'))
    except:
        A = []
    state = {'G':d.get('G'),'Q':d.get('Q'),'A':A,
             'ATTEN':d.get('ATTEN'), 'ENERGY':d.get('ENERGY')}

    # create scan data object
    data = ScanData(name=name,dims=[dims],scalers=scalers,
                    positioners=positioners,primary_axis=[paxis],
                    primary_det=pdet,state=state)
    return data

########################################################################
def spectra_scan(spectra,**kw):
    """
    Return ScanData instance holding a series of spectra.

    Parameters:
    -----------
    * spectra = string path name (maybe a prefix)
      or list of filenames

    Keywords and defaults:
    * med  = True -> return med's
    * xrf  = False -> return xrf's
    * sd   = True
      If sd == True then this fcn returns a scan data object
      Note it will create a dummy axis: x = num.arange(len(med))
      Otherwise either (XrfScan), (MedScan) or (XrfScan,MedScan)
      are returned depending on the med and xrf flags

    # files --> see the read function in med_data.py 
    * start   =  -1
    * end     =  -1
    * nfmt    =  3
    * fmt     = 'CARS'

    # For med data --> see med_data.py for more details
    * bad_mca_idx = []
    * total   = True
    * align   = True
    * correct = True
    * tau     = None

    # For xrf data --> see xrf_data.py for more details
    * det_idx = 0
    * emin    = -1.
    * emax    = -1.
    * xrf_params = {}
    * lines      = None
    
    """
    if kw.has_key('sd'): sd = kw.pop('sd')
    else: sd = True
    med = kw.get('med',True)
    xrf = kw.get('xrf',False)
    if med == False and xrf == False: return None

    # make sure we only read the files once!
    if med == True: kw['xrf'] = False
    if type(spectra) == types.StringType:
        name = os.path.split(spectra)[1].split('.')[0]
        spectra = _read_spectra(spectra,**kw)
    elif type(spectra) == types.ListType:
        name = os.path.split(spectra[0])[1].split('.')[0]
        kw['start']=-1
        spectra = _read_spectra(spectra,**kw)
    else:
        return None
    if spectra == None: return None
    #
    npts = len(spectra)
    med_spectra = None
    xrf_spectra = None
    if med == True and xrf == True:
        primary_det = 'xrf'
        med_spectra = spectra
        xrf_params  = kw.get('xrf_params',{})
        lines       = kw.get('lines')     
        det_idx     = kw.get('det_idx',0)
        emin        = kw.get('emin',-1.)
        emax        = kw.get('emax',-1.)
        xrf_spectra = xrf_data.med2xrf(med_spectra,xrf_params=xrf_params,
                                       lines=lines,det_idx=det_idx,
                                       emin=emin,emax=emax)
        xrf_spectra = xrf_data.XrfScan(xrf_spectra)
        med_spectra = med_data.MedScan(med_spectra)
        if sd == False: return (xrf_spectra,med_spectra)
    elif xrf == True and med == False:
        primary_det = 'xrf'
        xrf_spectra = xrf_data.XrfScan(spectra)
        if sd == False: return (xrf_spectra)
    elif med == True and xrf == False:
        primary_det = 'med'
        med_spectra = med_data.MedScan(spectra)
        if sd == False: return (med_spectra)
    x = num.arange(float(npts))
    data = ScanData(name=name,
                    dims = [npts],
                    positioners={'x':x},
                    primary_axis = 'x',
                    primary_det = primary_det,
                    med=med_spectra,
                    xrf=xrf_spectra)
    return data

########################################################################
def _read_spectra(fname,**kw):
    """
    read spectra files and return a list of med or xrf objects
    """
    start   = kw.get('start',-1)
    end     = kw.get('end',-1)
    nfmt    = kw.get('nfmt',3)
    fmt     = kw.get('fmt','CARS')
    #
    bad     = kw.get('bad_mca_idx',[])
    total   = kw.get('total',True)
    align   = kw.get('align',True)
    correct = kw.get('correct',True)
    tau     = kw.get('tau',None)
    #
    xrf     = kw.get('xrf',False)
    det_idx = kw.get('det_idx',0)
    emin    = kw.get('emin',-1.)
    emax    = kw.get('emax',-1.)
    xrf_params = kw.get('xrf_params',{})
    lines      = kw.get('lines',None)
        
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
                                          emin=emin,emax=emax,xrf_params=xrf_params,
                                          lines=lines,fmt=fmt)
        else:
            spectra = med_data.read_files(fname,start=start,end=end,nfmt=nfmt,
                                          bad_mca_idx=bad,total=total,align=align,
                                          correct=correct,tau=tau,fmt=fmt)
    else:
        if xrf:
            spectra = xrf_data.read_file(fname,bad_mca_idx=bad,total=total,align=align,
                                         correct=correct,tau=tau,det_idx=det_idx,
                                         emin=emin,emax=emax,xrf_params=xrf_params,
                                         lines=lines,fmt=fmt)
        else:
            spectra = med_data.read_file(fname,bad_mca_idx=bad,total=total,align=align,
                                         correct=correct,tau=tau,fmt=fmt)
    if spectra == None: return None
    if type(spectra) != types.ListType:
        spectra = [spectra]
    return spectra

########################################################################
def _spectra_range(self,fname):
    """
    Finds the range of numbered spectrum files.
    assume the fmt of files is: fname.nnn
    """
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
def image_scan(image,**kw):
    """
    Return scan data instance holding a series of images.

    Parameters:
    -----------
    * image = string path name (maybe prefix)
      or list of filenames

    Keywords and defaults:
    * sd   = True
      If sd == True then is returns a scan data object
      Note it will create a dummy axis: x = num.arange(len(image))
      Otherwise a ImageScan instance is returned

    # files
    * start = -1
    * end   = -1
    * nfmt  = 3
    * fmt   = 'tif'
    
    # for image data
    * rois  = None
    * archive is a dictionary with archive info
      archive['file'] = file name for image archive
      archive['path'] = path for image archive
      archive['setname'] = set name for image archive
      archive['descr'] = description of data for image archive

    Examples:
    >>s = image_scan('Image_file_pfx',start=0,end=5)
        
    """
    sd      = kw.pop('sd',True)
    rois    = kw.pop('rois',None)
    archive = kw.pop('archive',None)
    
    if type(image) == types.StringType:
        name  = os.path.split(image)[1].split('.')[0]
        image = _read_image(image,**kw)
    elif type(image) == types.ListType:
        name  = os.path.split(image[0])[1].split('.')[0]
        kw['start']=-1
        image = _read_image(image,**kw)
    else:
        return None
    if image == None: return None
    npts = len(image)
    image = image_data.ImageScan(image,rois=rois,archive=archive)
    if sd == False: return image
    x = num.arange(float(npts))
    data = ScanData(name=name,
                    dims = [npts],
                    positioners={'x':x},
                    primary_axis = 'x',
                    primary_det  = 'image',
                    image        = image)
    return data

########################################################################
def _read_image(fname,**kw):
    """
    read image files and return as a list of images
    """
    start = kw.get('start',-1)
    end   = kw.get('end',-1)
    nfmt  = kw.get('nfmt',3)
    fmt   = kw.get('fmt','tif')
    #
    if start > -1:
        if end == -1:
            ret=self._image_range(fname)
            if ret:
                (start,end) = ret
            else:
                print "No files found"
                return None
        image = image_data.read_files(fname,start=start,end=end,nfmt=nfmt)
    else:
        image = image_data.read(fname)
    if image == None: return None
    if type(image) != types.ListType:
        image = [image]
    return image

########################################################################
def _image_range(fname):
    """
    Finds the range of numbered spectrum files.
    assume the fmt of files is: fname.nnn
    """
    fname = fname + "_*.tif"
    files = glob(fname)
    n = len(files)
    if n>0:
        st = files[0].split('_')[-1]
        st = st.split('.')[0]
        en = files[n-1].split('_')[-1]
        en = en.split('.')[0]
        st = int(st)
        en = int(en)
        return (st,en)
    else:
        return None

########################################################################
class Reader:
    """
    Use the reader to read data files and return ScanData objects.
    
    The purpose of the reader is to hold parameters so that
    the data is processed in the same manner for each file/scan
    that is read.that should all be treated the same.
    
    For example we assume that if we read spectra, they all have the same
    set of params (bad dets', taus, etc..).  Same with images.

    This class knows about the following data/file types:
    * spec files
    * escan (to be added)
    * med file (CARS and other fmts)
    * xrf data
    * images 
    * other (ie ssrl files, old spec and super etc.. to be added)

    Attributes:
    -----------
    # Spec
    * spec_path is the path to locate spec files
    * spec_files is a list of spec files
    * spec_params is a dicitonary: {'image': False,'xrf':False,'med':False}

    # escan
    * escan_path is the path to locate escan files
    * escan_files is a list of escan files

    # med/xrf parameters
    * spectra_path is the path to locate spectra files
    * spectra_params is a dictionary:
      {'bad_mca_idx':[],'total':True,'align':True,
       'correct':True,'tau':None,'det_idx':0,
       'emin':-1.0,'emax':-1.0,'xrf_params':{},
       'lines':None,'fmt':'CARS','nfmt':3}

    # image parameters
    * image_path path to locate image files
    * image_params is a dictionary:
      {'rois':None,'fmt':'tif','nfmt':3,'archive':None}

    Examples:
    ---------
    >>r = scandata.Reader(spec='spec_file.spc')
    >>s1 = r.spec_scan(1,image=True)
    """
    ########################################################################
    def __init__(self,spec=None,spec_path=None,escan_path=None,
                 spectra_path=None,image_path=None):
        """
        Parameters:
        -----------
        * spec = spec files, string or list of strings
        * spec_path = string path for spec file locations
        * escan_path = string path for escan file locations
        * spectra_path = string path for med or xrf file locations
        * image_path = string path for image file locations
        """
        # Spec
        self.spec_path       = spec_path
        self.spec_files      = []
        self.spec_params     = {'image': False,'xrf':False,'med':False}

        # escan
        self.escan_path       = escan_path
        self.escan_files      = []

        # med/xrf parameters
        self.spectra_path   = spectra_path
        self.spectra_params = {'bad_mca_idx':[],'total':True,'align':True,
                               'correct':True,'tau':None,'det_idx':0,
                               'emin':-1.0,'emax':-1.0,'xrf_params':{},
                               'lines':None,'fmt':'CARS','nfmt':3}

        # image parameters
        self.image_path     = image_path
        self.image_params   = {'rois':None,'fmt':'tif','nfmt':3,
                               'archive':None}

        # load spec file(s) if passed
        if spec: self.read_spec(spec)

    ########################################################################
    def __repr__(self):
        """display"""
        lout = "Spec Files:"
        for s in self.spec_files:
            lout = "%s\n  %s, first=%i, last=%i" % (lout,s.fname,s.min_scan,s.max_scan)
        lout = "%s\nRead spec med   = %s, Read spec xrf = %s" % \
               (lout,self.spec_params['med'],self.spec_params['xrf'])
        lout = "%s\nRead spec image = %s" % (lout,self.spec_params['image'])
        lout = "%s\nPaths:"            % (lout)
        lout = "%s\n  Spec Path  = %s" % (lout,self.spec_path)
        lout = "%s\n  Med Path   = %s" % (lout,self.spectra_path)
        lout = "%s\n  Image Path = %s" % (lout,self.image_path)
        lout = "%s\n  Escan Path = %s" % (lout,self.escan_path)
        if self.spectra_params['lines'] != None:
            lout = "%s\n  Xrf Lines = %s" % (lout,str(self.spectra_params['lines']))
        if self.image_params['rois'] != None:
            lout = "%s\n  Image Rois = %s" % (lout,str(self.image_params['rois']))
        if self.image_params['archive'] != None:
            lout = "%s\n  Image archive = %s" % (lout,str(self.image_params['archive']))
        return lout

    ########################################################################
    def med_scan(self,fname,start=-1,end=-1,path=None):
        """
        Read (a collection of) med files into a scan data object

        Parameters:
        -----------
        * fname is the file name or prefix (or list of file names)
        * start and end are the initial and final numbers of
          file suffix if they are formatted numerically
        * path updates the med path setting

        Note:
        -----
        The data has a dummy axis: x = num.arange(len(med))
        """
        if path: self.spectra_path = path
        if self.spectra_path != None:
            fname = os.path.join(self.spectra_path,fname)

        data = spectra_scan(fname,med=True,xrf=False,sd=True,
                            start=start,end=end,
                            nfmt=self.spectra_params['nfmt'],
                            fmt=self.spectra_params['fmt'],
                            bad_mca_idx=self.spectra_params['bad_mca_idx'],
                            total=self.spectra_params['total'],
                            align=self.spectra_params['align'],
                            correct=self.spectra_params['correct'],
                            tau=self.spectra_params['tau'])
        return data
        
    ########################################################################
    def xrf_scan(self,fname,start=-1,end=-1,path=None,med=True):
        """
        Read (a collection of) xrf files into a scan data object

        Parameters:
        -----------
        * fname is the file name or prefix (or list of file names)
        * start and end are the initial and final numbers of
          file suffix if they are formatted numerically
        * path updates the med path setting
        * med is a flag that if True results in the ScanData object
          holding both med and xrf instances

        Notes:
        -----
        This data object holds both xrf and med data by default.
        The data has dummy axis: x = num.arange(len(xrf))
        """
        if path: self.spectra_path = path
        if self.spectra_path != None:
            fname = os.path.join(self.spectra_path,fname)
        data = spectra_scan(fname,med=med,xrf=True,sd=True,
                            start=start,end=end,
                            nfmt=self.spectra_params['nfmt'],
                            fmt=self.spectra_params['fmt'],
                            bad_mca_idx=self.spectra_params['bad_mca_idx'],
                            total=self.spectra_params['total'],
                            align=self.spectra_params['align'],
                            correct=self.spectra_params['correct'],
                            tau=self.spectra_params['tau'],
                            det_idx=self.spectra_params['det_idx'],
                            emin=self.spectra_params['emin'],
                            emax=self.spectra_params['emax'],
                            xrf_params=self.spectra_params['xrf_params'],
                            lines=self.spectra_params['lines'])
        return data

    ########################################################################
    def image_scan(self,fname,start=-1,end=-1,path=None):
        """
        Read (a collection of) image files into a scan data object

        Parameters:
        -----------
        * fname is the file name or prefix (or list of file names)
        * start and end are the initial and final numbers of
          file suffix if they are formatted numerically
        * path updates the image path setting

        Note:
        -----
        The data has a dummy axis: x = num.arange(len(med))
        """
        if path: self.image_path = path
        if self.image_path != None:
            fname = os.path.join(self.image_path,fname)
        data = image_scan(fname,sd=True,
                          start=start,end=end,
                          nfmt=self.image_params['nfmt'],
                          fmt=self.image_params['fmt'],
                          rois=self.image_params['rois'],
                          archive=self.image_params['archive'])
        return data

    ########################################################################
    def escan(self,fname):
        """
        not yet implemented
        """
        pass
    
    ########################################################################
    def ascii_scan(self,fname):
        """
        not yet implemented
        """
        pass

    ########################################################################
    def read_spec(self,spec,path=None):
        """
        Add a spec file (or files) to the reader

        The spec files are read and cached by the
        reader instance.  See list_spec and spec_scan
        for more information 

        Paramters:
        ----------
        * spec is a string file name (or list of file names)
        * path updates the spec path setting
        """
        if path != None: self.spec_path = path
        if type(spec) == types.StringType: spec = [spec]
        
        for s in spec:
            sfile = self._spec(file=s)
            sfile.read()

    ########################################################################
    def list_spec(self, show=True):
        """
        Return a list of spec file scans
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
        Get data from a spec scan

        Parameters:
        -----------
        * scan  = spec scan number
        * file  = spec file name (default = first file read in)
        * med   = read med files, True/False (default False)
        * xrf   = read xrf files, True/False (default False)
        * image = read image files, True/False (default False)

        Notes:
        ------
        When spec files are read in (see read_spec) they get appended
        to a list. To access them using this method pass in the file
        name and scan number.  If the file has been updated since last read,
        it will be re-read.  If you have only read in a single file, or you
        want a scan from the first file read in then you dont need to
        supply the file argument.  If you would like to resort the list
        of spec files you can use something like the following to move the
        third file to the top of the list:
        >>reader.spec_files.insert(0,reader.spec_files.pop(2))

        Also note that if you are reading med/xrf or image files and an
        explicit path has not been set, then one will be computed using
        an assumed path relative to the spec files.  See the
        _spec_spectr_path and _spec_image_path methods
        """
        if med!=None: self.spec_params['med'] = med
        if xrf!=None: self.spec_params['xrf'] = xrf
        if image!=None: self.spec_params['image'] = image
        med = self.spec_params['med']
        xrf = self.spec_params['xrf']
        img = self.spec_params['image']

        # File
        spec = self._spec(file=file)
        if not spec: return None
        if scan == None: return None
        data = spec_scan(spec,scan)
        if not data: return None
        # Spectra
        if med or xrf:
            spectra_pfx = self._spec_spectra_path(spec,scan)
            start = 0
            end   = data.dims[0] - 1
            spectra = spectra_scan(spectra_pfx,sd=False,
                                   med=med,xrf=xrf,
                                   start=start,end=end,
                                   nfmt=self.spectra_params['nfmt'],
                                   fmt=self.spectra_params['fmt'],
                                   bad_mca_idx=self.spectra_params['bad_mca_idx'],
                                   total=self.spectra_params['total'],
                                   align=self.spectra_params['align'],
                                   correct=self.spectra_params['correct'],
                                   tau=self.spectra_params['tau'],
                                   det_idx=self.spectra_params['det_idx'],
                                   emin=self.spectra_params['emin'],
                                   emax=self.spectra_params['emax'],
                                   xrf_params=self.spectra_params['xrf_params'],
                                   lines=self.spectra_params['lines'])
            if spectra == None:
                print "Warning, med files not read"
            elif med == True and xrf == True:
                (data.xrf,data.med) = spectra
                if len(data.xrf.xrf) != data.dims[0] or \
                   len(data.med.med) != data.dims[0]:
                    print "Warning, number of spectra dont match scan"
            elif med == True and xrf == False:
                data.med = spectra
                if len(data.med.med) != data.dims[0]:
                    print "Warning, number of spectra dont match scan"
            elif med == False and xrf == True:
                data.xrf = spectra
                if len(data.xrf.xrf) != data.dims[0]:
                    print "Warning, number of spectra dont match scan"
        # Images
        if img:
            image_pfx = self._spec_image_path(spec,scan)
            start = 0
            end   = data.dims[0] - 1
            if self.image_params['archive'] != None:
                #if self.image_params['archive'].has_key('file') == False:
                imfile = spec.fname + '_images.h5'
                self.image_params['archive']['file'] = imfile
                self.image_params['archive']['setname'] = 'S%03d' % int(scan)
            image = image_scan(image_pfx,sd=False,
                               start=start,end=end,
                               nfmt=self.image_params['nfmt'],
                               fmt=self.image_params['fmt'],
                               rois=self.image_params['rois'],
                               archive=self.image_params['archive'])
            if image == None:
                print "Warning, image files not read"
            else:
                data.image = image
                if len(data.image.image) != data.dims[0]:
                    print "Warning, number of images dont match scan"
        
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
        
        Notes:
        ------
        * This is IDC 13 specific!
            if self.spectra_path == None:
              fpfx = spec_path/xrf_files/spec_pfx/nnn/spec_file.spc_nnn
            else:
              fpfx = self.spectra_path/spec_file.spc_nnn
          With nnn = scan number

        * The read function will add '.point_number' to all
          the file names...
        """
        sc_num = '%03d' % int(scan)
        fpfx = "%s_%s" % (spec.fname, sc_num)
        if self.spectra_path == None:
            spec_pfx = spec.fname.rsplit('.',1)[0]
            path = os.path.join(spec.path,'xrf_files',spec_pfx,sc_num)
            fpfx = os.path.join(path,fpfx)
        else:
            fpfx = os.path.join(self.spectra_path,fpfx)
        return fpfx
    
    ########################################################################
    def _spec_image_path(self,spec,scan):
        """
        Compute path for spectra associated with spec scan
        
        Notes:
        ------
        * This is IDC 13 specific!
            if self.spectra_path == None:
              fpfx = spec_path/images/spec_pfx/Snnn/spec_file.spc_Snnn
            else:
              fpfx = image_path/spec_file.spc_Snnn
          With nnn = scan number

        * The read function will add '_point_num.tif' to all
          the file names...
        """
        sc_num = '%03d' % int(scan)
        fpfx = "%s_S%s" % (spec.fname, sc_num)
        if self.image_path == None:
            spec_pfx = spec.fname.rsplit('.',1)[0]
            path = os.path.join(spec.path,'images',spec_pfx,'S'+sc_num)
            fpfx = os.path.join(path,fpfx)
        return fpfx
    
############################################################################
############################################################################
if __name__ == '__main__':
    pass
