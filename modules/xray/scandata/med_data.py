##############################################################################
"""
T. Trainor (fftpt@uaf.edu)
Functions for operating on / working with
med's in scandata

Modifications:
--------------

"""
##############################################################################

import types
import numpy as num
from matplotlib import pyplot

from detector import deadtime
from detector import medfile_cars
from detector import medfile_emsa
import xrf_data

##############################################################################
def read(file,bad_mca_idx=[],total=True,align=True,
         correct=True,tau=None,fmt='CARS'):
    """
    Read detector files
    >>m = med.read(file="file_name",bad_mca_idx=[],
                   total=True,align=True,tau=None)
    """
    if fmt == 'CARS':
        rd = medfile_cars.read_med
    elif fmt == 'EMSA':
        rd = medfile_emsa.read_med
    else:
        return
    if type(file) == types.StringType:
        med = rd(file=file,bad_mca_idx=bad_mca_idx,
                 total=total,align=align,correct=correct,tau=tau)
    elif type(file) == types.ListType:
        med = []
        for f in file:
            tmp = rd(file=f,bad_mca_idx=bad_mca_idx,
                     total=total,align=align,correct=correct,tau=tau)
            med.append(tmp)
    else:
        return None
    return med

##############################################################################
def read_files(prefix,start=0,end=100,nfmt=3,bad_mca_idx=[],
               total=True,align=True,correct=True,tau=None,
               fmt='CARS'):
    """
    Read multiple files
    returns a list of med objects
    """
    if fmt == 'CARS':
        rd = medfile_cars.read_med_files
    elif fmt == 'EMSA':
        rd = medfile_emsa.read_med_files
    else:
        return
    med = rd(prefix,start=start,end=end,nfmt=nfmt,
             bad_mca_idx=bad_mca_idx,total=total,
             align=align,correct=correct,tau=tau)
    return med

########################################################################
def med_plot(data,scan_pnt=0,hold=False,ylog=True):
    """
    plot med spectra
    """
    tiny = 1.e-9
    en = data.med[scan_pnt].get_energy()
    ch = data.med[scan_pnt].get_data()
    bad = data.med[scan_pnt].bad_mca_idx
    nchan = len(ch)
    if not hold:
        pyplot.clf()
    if nchan == 1:
        pyplot.plot(en[0],ch[0])
    else:
        for j in range(nchan):
            if j not in bad:
                label = str(j)
                pyplot.plot(en[j],ch[j]+tiny,label=label)
        pyplot.legend()
    if ylog:
        pyplot.semilogy()
        pyplot.ylim(ymin=1.)

########################################################################
class MedScan:
    """
    Class to hold a collection of meds associated with a scan
    ie one med per scan point

    Keywords:
    bad_mca_idx: A list of bad mca's, data will be zeros.  An empty
               list (default) means all the detectors are ok.
               Note detector indexing starts at zero!

    total: Set this keyword to return the sum of the spectra from all
           of the Mcas as a 1-D numeric array.
        
    align: Set this keyword to return spectra which have been shifted and
           and stretched to match the energy calibration parameters of the
           first (good) detector.  This permits doing arithmetic on a
           "channel-by-channel" basis. This keyword can be used alone
           or together with the TOTAL keyword, in which case the data
           are aligned before summing.
           
    correct:
        True means to apply deadtime correction, false ignores it

    tau:  mca deadtime tau values
           None --> recompute correction factor
           []   --> Turn off correction, ie set taus to -1
           single value (or single valued list) --> assign to all mcas
           list (or array) --> assign to individual mcas
    
    """
    ########################################################################
    def __init__(self,med=[],**kws):
        if type(med) != types.ListType: med = [med]
        self.med = med
        self.init_params(params=kws)

    ########################################################################
    def init_params(self,params={}):
        """
        (Re)initialize all the med params to those given
        """
        if len(params) > 0:
            for m in self.med:
                m.init_params(params)

    ################################################################
    def update_tau(self,tau):
        """
        update med tau factors
        """
        npnt = len(self.med)
        for j in range(npnt):
            self.med[j].update_correction(tau)

    ################################################################
    def get_tau(self):
        """
        return tau values
        """
        return self.med[0].get_tau()
        tau = num.zeros(ndet)
        for j in range(ndet):
            tau[j] = self.med[0].mca[j].tau
        return tau

    ################################################################
    def ocr(self):
        """
        return outgoing count rate (ocr) values from meds
        """
        ndet = self.med[0].n_detectors
        npnt = len(self.med)
        ocr  = num.zeros((npnt,ndet))
        for j in range(npnt):
            for k in range(ndet):
                tot = self.med[j].mca[k].total_counts
                lt  = self.med[j].mca[k].live_time
                ocr[j][k] = float(tot)/float(lt)
        return num.transpose(ocr)

    ################################################################
    def med2xrf(self,xrf_params={},lines=None,det_idx=0,emin=-1.,emax=-1.):
        """
        convert med objects to xrf objects.
        returns an XrfScan object
        """
        xrf = xrf_data.med2xrf(self.med,xrf_params=xrf_params,
                               lines=lines,det_idx=det_idx,
                               emin=emin,emax=emax)
        #return xrf
        return xrf_data.XrfScan(xrf)

########################################################################
########################################################################
if __name__ == '__main__':
    pass
