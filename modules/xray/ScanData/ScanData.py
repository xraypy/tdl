# ScanData
# T. Trainor (fftpt@uaf.edu)
# Class to store and operate on generic 'scan data'
# Note spec.get_scan should ret one of these....
#################

import XRF
import CarsMcaFile
import os, os.path

#################

"""
ScanData Class holds:
name = scan name
scan_dims -> would be [npts1] for a 1-d scan
            [npts1,npts2] for a 2-d scan
             where npts1 is the outer loop
             and npts1 is for the inner loop
mcas = []                    -> array has dims = scan_dims or is [] for No mca
scalars = {'io':[],'i1':[]}  -> each array has dims = scan_dims
positioners  = {'E':[]}      -> these may be single values, or same size as scalars
primary_axis = ['E']         -> if multi-dim, list of outer, inner
primary_det  = 'i1'          -> use for default plotting
"""

class ScanData:
    def __init__(self,name='',scan_dims=[],spectra=[],scalars={},
                 positioners={},primary_axis=[],primary_det=None):
        self.name         = name
        self.scan_dims    = scan_dims
        self.spectra      = spectra
        self.scalars      = scalars
        self.positioners  = positioners
        self.primary_axis = primary_axis
        self.primary_det  = primary_det

    def __repr__(self):
        lout = "Scan Name = %s\n" % self.name
        lout = lout + "Scan Dimensions = %s\n" % str(self.scan_dims)
        lout = lout + "Scalars:\n"
        for sc in self.scalars.keys():
            lout = lout + "    %s" % sc
        lout = lout + "\nPositioners:\n"
        for p in self.positioners.keys():
            lout = lout + "    %s" % p
        if self.spectra != []:
            lout = lout + "\nScan includes spectra data:\n"
        else:
            lout = lout + "\nScan does not include mca data:\n"
        lout = lout + "Primary scan axis = %s\n"  % str(self.primary_axis)
        lout = lout + "Primary detector = %s\n"  % self.primary_det
        return lout

    def get_scalar(self,label=None):
        if label == None:
            label = self.primary_det[0]
        return self.scalars.get(label)
        
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
        ocr = []
        for s in self.spectra:
            cts = s.get_count_totals()
            det_ocr = []
            for j in range(ndet):
                det_ocr.append(cts[j]['OCR'])
            ocr.append(det_ocr)
        return ocr

##############################################################################
def read_xrf_scan(first='',bad_mca_idx=[],
                  total=True,align=True,tau=None):
    fname = first
    spectra = []
    while os.path.exists(fname):
        xrf =  XRF.read_xrf_file(file=fname,bad_mca_idx=bad_mca_idx,
                                    total=total,align=align,tau=tau)
        spectra.append(xrf) 
        fname = CarsMcaFile.increment_filename(fname)

    npts = len(spectra)
    sd = ScanData(name=first,scan_dims=[npts],spectra=spectra)
    return sd


def get_spectrum(sd,pnt,**kw):
    return sd.get_spectrum(pnt)

def get_spectra_OCR(sd,**kw):
    return sd.get_spectra_OCR()

##############################################################################
# tell tdl to create these groups:
_groups_ = [('scan',True)]

_func_ = {"scan.xrf_scan":(read_xrf_scan,None),
          "scan.get_spectrum":(get_spectrum,None),
          "scan.get_ocr":(get_spectra_OCR,None)}

#
#_var_  = {"x":[1,2,3],
#          "test.dat":[1.3,6.5]}
#
# code to run on initialization (no args, but will get a 'tdl reference')
#_init_ = libtest_init
