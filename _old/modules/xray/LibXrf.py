# T. Trainor 
# Wrapper functions for XRF
# T. Trainor, 8-11-2006
#################################################################################################

##########################################################################
from Num import Num, num_version

import os
import sys
import types
import copy
from Util import datalen

import XRF


title = "CARS XRF Library "

##############################################################################
def read_xrf(file='',bad_mca_idx=[],total=True,align=True,tau=None,
             det_idx=0,emin=-1.0,emax=-1.0,fmt='CARS',xrf_params={},
             tdl=None,**kws):
    """
    Read detector files
    >>m = xrf.read(file="file_name",bad_mca_idx=[],total=True,align=True,tau=None)

    Returns an xrf object.  This function is appropriate for reading single and
    multi-element detectors.  (We always assume that the detector may be a
    multi-element detector)

    Inputs:
        file:
            File name
            
        bad_mca_idx:
            A list of bad detectors to be used.  An empty
            list (default) means use all the detectors.  Note detector indexing
            starts at zero!

        total:
            Set this keyword to toggle the total flag in xrf.  Processing will work
            on the sum of detectors.
            
        align:
            Set this keyword to return spectra which have been shifted and
            and stretched to match the energy calibration parameters of the
            first detector.

        tau:
            None or a list of deadtime factors for each detector (dimesion should
            be the same as the total number of detectors)
                     ocr = icr * exp(-icr*tau)
                     cor = (icr/ocr)*(rt/lt)
                     max_icr = 1/tau
                     max_ocr = max_icr*exp(-1)
    """
    xrf = XRF.read_xrf(file=file,bad_mca_idx=bad_mca_idx,total=total,align=align,
                       tau=tau,det_idx=det_idx,emin=emin,emax=emax,fmt=fmt,
                       xrf_params=xrf_params)
    return xrf

##############################################################################
def read_xrf_files(sd=None,prefix='xrf',start=0,end=100,nfmt=3,
                   bad_mca_idx=[],total=True,align=True,tau=None,
                   det_idx=0,emin=-1.0,emax=-1.0,fmt='CARS',xrf_params={}):
    """
    read multiple xrf spectra
    """
    spectra = read_xrf_files(prefix=prefix,start=start,end=end,nfmt=nfmt,
                   bad_mca_idx=bad_mca_idx,total=total,align=align,tau=tau,
                   det_idx=det_idx,emin=emin,emax=emax,fmt=fmt,
                   xrf_params=xrf_params)

    if sd == None:
        npts = len(spectra)
        sd = ScanData.ScanData(name=prefix,scan_dims=[npts],spectra=spectra)
        return sd
    else:
        sd.spectra = spectra
        return sd

##################################################################################

_groups_ = [('xrf',True)]

# tdl functions
_func_ = {"xrf.read":(read_xrf),
          "xrf.read_scan":read_xrf_files}

#_scripts_ = ['xrf_scripts.tdl']

