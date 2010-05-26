"""
Functions for operating on / working with
med's in scandata

Authors/Modifications:
----------------------
* T. Trainor (tptrainor@alaska.edu)

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

    Parameters:
    -----------
    * file is the med filename
    * bad_mca_idx is a list of bad detector elements (ie to ignore)
      indexing starts at zero.
    * total is a flag to indicate if data should be summed when returned
    * correct is a flag to indicate if deadtime corrections should be
      applied when data is returned
    * align is a flag to indicate that data should be aligned based
      on individual element calibrations
    * tau are detector deadtime parameters (one number would apply to all
      or a list with a tau value for each detector element)
    * fmt is a string indicating file format = 'CARS' or 'EMSA'

    Returns:
    -------
    * med object

    Examples:
    ---------
    >>m = med.read(file="file_name",bad_mca_idx=[],
                   total=True,align=True,tau=None)
    >>d = m.get_data()
    >>e = m.get_energy()
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
    Read multiple files with numerical suffix in file names

    See read_file    
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

    data is a ScanData or MedScan object    
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

    Attributes:
    -----------
    * med is a list of med objects

    """
    ########################################################################
    def __init__(self,med=[],**kws):
        """
        Parameters:
        * med is a list of med objects

        The remaining (optional) keyword parameters are applied to
        each med object in the list
        
        * bad_mca_idx: A list of bad mca's, data will be zeros.  An empty
          list (default) means all the detectors are ok.
          Note detector indexing starts at zero!

        * total: Set this keyword to return the sum of the spectra from all
          of the Mcas as a 1-D numeric array.
            
        * align: Set this keyword to return spectra which have been shifted and
          and stretched to match the energy calibration parameters of the
          first (good) detector.  This permits doing arithmetic on a
          "channel-by-channel" basis. This keyword can be used alone
          or together with the TOTAL keyword, in which case the data
          are aligned before summing.
               
        * correct: True means to apply deadtime correction, false ignores it

        * tau:  mca deadtime tau values
          None --> recompute correction factor
          []   --> Turn off correction, ie set taus to -1
          single value (or single valued list) --> assign to all mcas
          list (or array) --> assign to individual mcas
        """
        if type(med) != types.ListType: med = [med]
        self.med = med
        self.init_params(params=kws)

    ########################################################################
    def init_params(self,params={}):
        """
        (Re)initialize all the med params

        Parameters:
        -----------
        params is a dictionary of paramaters that may include the following:

        * bad_mca_idx: A list of bad mca's, data will be zeros.  An empty
          list (default) means all the detectors are ok.
          Note detector indexing starts at zero!

        * total: Set this keyword to return the sum of the spectra from all
          of the Mcas as a 1-D numeric array.
            
        * align: Set this keyword to return spectra which have been shifted and
          and stretched to match the energy calibration parameters of the
          first (good) detector.  This permits doing arithmetic on a
          "channel-by-channel" basis. This keyword can be used alone
          or together with the TOTAL keyword, in which case the data
          are aligned before summing.
               
        * correct: True means to apply deadtime correction, false ignores it

        * tau:  mca deadtime tau values
          None --> recompute correction factor
          []   --> Turn off correction, ie set taus to -1
          single value (or single valued list) --> assign to all mcas
          list (or array) --> assign to individual mcas

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
        #tau = num.zeros(ndet)
        #for j in range(ndet):
        #    tau[j] = self.med[0].mca[j].tau
        #return tau

    ################################################################
    def ocr(self):
        """
        return outgoing count rate (ocr) values from meds

        Outputs:
        --------
        * the returned value is an array of ocr values.  the first index
          is the detector index. the second index is the point

        Example:
        -------
        >>ocr = medscan.ocr()
        >>ocr_detector_1 = ocr[1]
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

        Parameters:
        -----------
        * xrf_params - see XrfData object
        * lines is a list of xrf lines
        * det_idx is the detector idx to be used to generate the data
          if the detector is not summed (for an med) you can use a sepcified
          detector to make the xrf object
        * emin and emax are the energy range in keV
        
        Returns:
        -------
        * an XrfScan object
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
