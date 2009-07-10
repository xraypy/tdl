#######################################################################
"""
Tom Trainor (fftpt@uaf.edu)
Class to read in generic 'scan data'

Modifications:
--------------

"""
#######################################################################
"""
Todo
- add escan files
- other files...

"""
########################################################################

import types
from   glob import glob
import os, copy
import numpy as Num

from   specfile import SpecFile
from   data import ScanData
import image_data
import xrf_ops

########################################################################
def spec_scan(spec,sc_num):
    """
    return scan data from a specfile instance
    """
    d = spec.scan_dict(sc_num)
    scalers = {}
    positioners = d['P']
    for key in d['data'].keys():
        if key in positioners.keys():
            positioners[key] = Num.array(d['data'][key])
        else:
            scalers[key] = Num.array(d['data'][key])
    name  = d['file'] + ' Scan ' + str(int(sc_num))
    dims  = d['nrow']
    paxis = d['labels'][0]
    pdet  = d['labels'][-1]

    state = {'G':d['G'],'Q':d['Q'], 'ATTEN':d['ATTEN']}

    sd = ScanData(name=name,dims=[dims],scalers = scalers,
                  positioners=positioners,primary_axis=[paxis],
                  primary_det=pdet,state=state)
    return sd


########################################################################
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
     * spec = spec files. string or list of strings
     * spec_path = string path for spec file locations
     * escan_path = string path for escan file locations
     * med_path = string path for escan file locations
     * image_path = string path for escan file locations

    Examples:
     >>r = scandata.Reader(spec='spec_file.spc')
     >>s1 = r.spec_scan(1,img=True)

    """
    ########################################################################
    def __init__(self,spec=None,spec_path=None,escan_path=None,
                 med_path=None,image_path=None):
        # Spec
        self.spec_path       = spec_path
        self.spec_files      = []

        # escan
        self.escan_path       = escan_path
        self.escan_files      = []

        # med/xrf parameters
        self.med_path       = med_path
        self.med            = False
        self.med_params     = {'bad_mca_idx':[],'total':True,'align':True,
                               'correct':True,'tau':None,'det_idx':0,
                               'emin':-1.0,'emax':-1.0,
                               'fmt':'CARS','nfmt':3}
        self.xrf            = False
        self.xrf_params     = {}
        self.xrf_lines      = None

        # image parameters
        self.image_path     = image_path
        self.img            = False
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
        lout = "%s\nRead image = %s"   % (lout,self.img)
        lout = "%s\nPaths:"            % (lout)
        lout = "%s\n  Spec Path  = %s" % (lout,self.spec_path)
        lout = "%s\n  Med Path   = %s" % (lout,self.med_path)
        lout = "%s\n  Image Path = %s" % (lout,self.image_path)
        lout = "%s\n  Escan Path = %s" % (lout,self.escan_path)
        if self.xrf_lines:
            lout = "%s\n  Xrf Lines = %s" % (lout,str(self.xrf_lines))
        if self.image_rois:
            lout = "%s\n  Image Rois = %s" % (lout,str(self.image_rois))
            
        return lout

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
    def spec_scan(self,scan,file=None,med=None,xrf=None,img=None):
        """
        get data from a spec scan
        """
        # keywords
        if med!=None: self.med = med
        if xrf!=None: self.xrf = xrf
        if img!=None: self.img = img

        # file
        spec = self._spec(file=file)
        if not spec: return None
        if scan == None: return None

        data = spec_scan(spec,scan)
        if not data: return None

        # set more parameters from reader
        data.xrf_lines  = self.xrf_lines
        data.image_rois = self.image_rois
        
        # Spectra
        if self.med or self.xrf:
            fmt_scan_num = '%03d' % int(scan)
            med_pfx = spec.fname
            med_pfx = "%s_%s" % (med_pfx, fmt_scan_num)
            if self.med_path == None:
                spec_pfx = spec.fname.rsplit('.',1)[0]
                path = os.path.join(spec.path,spec_pfx)
                path = os.path.join(path,fmt_scan_num)
                med_pfx = os.path.join(path,med_pfx)
            # get range for scan... ie first and last idx
            start = 0
            end   = data.dims[0] - 1

        if self.med:
            med = self._read_spectra(med_pfx,start=start,end=end,xrf=False)
            if med == None:
                print "Warning, med files not read"
            else:
                data.med = med
            if len(med) !=  data.dims[0]:
                print "Warning, number of spectra dont match scan"

        if self.xrf:
            if len(data.med)>0:
                data.med2xrf(xrf_params=self.xrf_params,
                             det_idx=self.med_params['det_idx'],
                             emin=self.med_params['emin'],
                             emax=self.med_params['emax'])
            else:
                xrf = self._read_spectra(med_pfx,start=start,end=end,xrf=True)
                if xrf == None:
                    print "Warning, xrf files not read"
                else:
                    data.xrf = xrf
                if len(xrf) !=  data.dims[0]:
                    print "Warning, number of spectra dont match scan"

        # Images
        if self.img:
            fmt_scan_num = '%03d' % int(scan)
            #image_pfx = spec.fname.rsplit('.',1)[0]
            #image_pfx = "%s_%s" % (image_pfx, fmt_scan_num)
            image_pfx = "%s_S%s" % (spec.fname, fmt_scan_num)
            if self.image_path == None:
                spec_pfx = spec.fname.rsplit('.',1)[0]
                path = os.path.join(spec.path,'images',spec_pfx)
                fmt_scan_num2 = 'S'+fmt_scan_num
                path = os.path.join(path,fmt_scan_num2)
                #path = os.path.join(path,fmt_scan_num)
                tmp_xxx = "S%s" % fmt_scan_num
                path = os.path.join(path,tmp_xxx)
                image_pfx = os.path.join(path,image_pfx)
            # get range for scan... ie first and last idx
            start = 0
            end   = data.dims[0] - 1
            images = self._read_image(image_pfx,start=start,end=end)
            data.image = images
            
        return data

    ########################################################################
    def _spec(self,file=None):

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
    def read_med(self,fname,start=-1,end=-1,path=None):
        """
        read med files
        """
        if path: self.med_path = path
            
        med = self._read_spectra(fname,start=start,end=end,xrf=False)
        if med == None: return None
        data = ScanData(name=fname,
                        dims = [len(med)],
                        primary_axis = 'med',
                        primary_det = 'med',
                        med=med,
                        xrf_lines = self.xrf_lines,
                        image_rois = self.image_rois)
        return data
        
    ########################################################################
    def read_xrf(self,fname,start=-1,end=-1,path=None):
        """
        read xrf files
        """
        if path: self.med_path = path
            
        xrf = self._read_spectra(fname,start=start,end=end,xrf=True)
        if xrf == None: return None
        data = ScanData(name=fname,
                        dims = [len(xrf)],
                        primary_axis = 'xrf',
                        primary_det = 'xrf',
                        xrf=xrf,
                        xrf_lines = self.xrf_lines,
                        image_rois = self.image_rois)
        return data

    ########################################################################
    def _read_spectra(self,fname,start=-1,end=-1,xrf=True):
        """
        read spectra files
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
        if xrf:
            xrf_params = self.xrf_params
            xrf_lines = self.xrf_lines
        else:
            xrf_params = None
            xrf_lines = None

        if self.med_path != None:
            fname = os.path.join(self.med_path,fname)
            
        if start > -1:
            if end == -1:
                ret = self._spectra_range(fname)  
                if ret:
                    (start,end) = ret
                else:
                    print "No files found"
                    return None
            spectra = xrf_ops.read_files(fname,start=start,end=end,nfmt=nfmt,
                                         bad_mca_idx=bad,total=total,align=align,
                                         correct=correct,tau=tau,det_idx=det_idx,
                                         emin=emin,emax=emax,fmt=fmt,
                                         xrf_params=xrf_params,lines=xrf_lines)
        else:
            spectra = xrf_ops.read_file(fname,bad_mca_idx=bad,total=total,align=align,
                                        correct=correct,tau=tau,det_idx=det_idx,
                                        emin=emin,emax=emax,fmt=fmt,
                                        xrf_params=xrf_params,lines=xrf_lines)
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
    def read_image(self,fname,start=-1,end=-1,path=None):
        """
        read image files
        """
        if path: self.image_path = path
        image = self._read_image(fname,start=start,end=end)
        if image == None: return None
        data = ScanData(name=fname,
                        dims = [len(image)],
                        primary_axis = 'image',
                        primary_det  = 'image',
                        image        = image,
                        image_rois   = self.image_rois)
        return data

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
    
    ########################################################################
    def read_escan(self,scan):
        pass
    
    ########################################################################
    def read_column(self,scan):
        pass

########################################################################
