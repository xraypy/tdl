# ScanData
# T. Trainor (fftpt@uaf.edu)
# Class to store and operate on generic 'scan data'
# Note spec.get_scan should ret one of these....
#################

import types
import numpy as Num

#################

"""
ScanData Class holds:
name = scan name
scan_dims -> would be [npts1] for a 1-d scan
             [npts1,npts2] for a 2-d scan
             where npts1 is the outer loop
             and npts1 is for the inner loop
scalers = {'io':[],'i1':[]}  -> each array has dims = scan_dims
positioners  = {'E':[]}      -> these may be single values, or same size as scalers
primary_axis = ['E']         -> if multi-dim, list of outer, inner
primary_det  = 'i1'          -> use for default plotting
spectra = []                 -> array has dims = scan_dims or is [] for No mca
images = []                  -> array has dims = scan_dims or is [] for No images
state = {}                   -> arbitraty dictionary of additional state information

Need to add methods for extracting, merging, appending
"""

class ScanData:
    def __init__(self,name='',scan_dims=[],scalers={},positioners={},
                 primary_axis=[],primary_det=None,spectra=[], images=[],state={}):
        self.name         = name
        self.scan_dims    = scan_dims
        self.scalers      = scalers
        self.positioners  = positioners
        self.primary_axis = primary_axis
        self.primary_det  = primary_det
        self.spectra      = spectra
        self.images       = images
        self.state        = state

    """
    def __setitem__(self,arg,val):
        # note this currenlty craps for sub arrays.
        # ie x[0][0].  In such a case arg is a tuple
        # first is the name, the rest are the subscripts
        #print arg, type(arg)
        #print val, type(val)
        if type(arg) == types.StringType:
            if arg.startswith("_"):
                raise TypeError, "Private attribute"
        #if hasattr(self,arg):
        #    setattr(self,arg,val)
        setattr(self,arg,val)
        return
    
    def __getitem__(self,arg):
        #print arg
        if type(arg) == types.StringType:
            if arg.startswith("_"):
                raise TypeError, "Private attribute"
        return getattr(self,arg)
    """
    
    def __getitem__(self,arg):
        #print arg
        
        # order matters,
        # see if match with spectra or images then 
        # scalers, positioners, state
        if type(arg) == types.StringType:
            if arg == 'spectra':
                return self.spectra
            if arg == 'images':
                return self.images
            for key in self.scalers.keys():
                if key == arg:
                    return self.scalers[key]
            for key in self.positioners.keys():
                if key == arg:
                    return self.positioners[key]
            for key in self.state.keys():
                if key == arg:
                    return self.state[key]

        return None

    def __repr__(self):
        lout = "* Scan Name = %s\n" % self.name
        lout = lout + "* Scan Dimensions = %s\n" % str(self.scan_dims)

        lout = lout + "* Scalers:\n"
        ct = 0
        for sc in self.scalers.keys():
            lout = lout + "%12s" % sc
            ct = ct + 1
            if ct > 7:
                lout = lout + '\n'
                ct = 0

        lout = lout + "\n* Positioners:\n"
        ct = 0
        for p in self.positioners.keys():
            lout = lout + "%12s" % p
            ct = ct + 1
            if ct > 7:
                lout = lout + '\n'
                ct = 0

        if len(self.state) > 0:
            lout = lout + "\n* Additional State Variables:\n"  
            ct = 0
            for p in self.state.keys():
                lout = lout + "%12s" % p
                ct = ct + 1
                if ct > 7:
                    lout = lout + '\n'
                    ct = 0

        lout = lout + "\n* Primary scan axis = %s\n"  % str(self.primary_axis)
        lout = lout + "* Primary detector = %s\n"  % self.primary_det

        if self.spectra != []:
            lout = lout + "* Scan includes %i spectra\n" % len(self.spectra)
        else:
            lout = lout + "* Scan does not include spectra\n"

        if self.images != []:
            lout = lout + "* Scan includes %i images\n" % len(self.images)
        else:
            lout = lout + "* Scan does not include images\n"
        
        return lout
    
    def get_scaler(self,label=None):
        if label == None:
            label = self.primary_det[0]
        return self.scalers.get(label)
        
    def get_positioner(self,label=None):
        if label == None:
            label = self.primary_axis[0]
        return self.positioners.get(label)

    def get_spectrum(self,pnt=0):
        pnt = int(pnt)
        if self.spectra == []:
            return None
        if pnt not in range(self.scan_dims[0]):
            return None
        return self.spectra[pnt]

    def get_spectra_OCR(self):
        ndet = self.spectra[0].ndet
        npnt = len(self.spectra)
        ocr = Num.zeros((npnt,ndet))
        for j in range(npnt):
            cts = self.spectra[j].get_count_totals()
            for k in range(ndet):
                ocr[j][k] = cts[k]['OCR']
        return Num.transpose(ocr)
        #return ocr

##############################################################################
