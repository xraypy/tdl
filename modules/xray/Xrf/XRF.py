#######################################################################
"""
XRF class for data handling and fitting
T. Trainor, fftpt@uaf.edu

--------------
 Modifications
--------------

"""
#############################################################################
"""
Todo

- Improve above doc and below doc/comments
- Test!

"""
########################################################################
"""
Inputs for Med objects:
    file:
        File name or list of file names
        
    bad_mca_idx:
        A list of bad detectors.  An empty list (default) means use all
        the detectors.  Note detector indexing starts at zero!

    total:
        Set this keyword to toggle the total flag in xrf.  Processing will work
        on the sum of detectors.
        
    align:
        Set this keyword to return spectra which have been shifted and
        and stretched to match the energy calibration parameters of the
        first detector.

    tau:
        List of deadtime factors for each detector.  Dimesion should
        be the same as the total number of detectors, or empty list for no
        deatdtime correction
                 ocr = icr * exp(-icr*tau)
                 cor = (icr/ocr)*(rt/lt)
                 max_icr = 1/tau
                 max_ocr = max_icr*exp(-1)

    det_idx:
        If you chose to work on a single detector (and total is not True)

    emin/emax:
        Set emin/emax of data range to work on (keV)

    fmt:
        'CARS' or 'EMSA' files


Inputs for Xrf objects:

    xrf_params:
        Format is that returned by Xrf.get_params:

    xrf_params = {'fit':fit_par,'bgr':bgr_par,'pk':peak_par}

    fit_par = {'energy_offset':XrfSpectrum.energy_offset,
               'energy_slope': XrfSpectrum.energy_slope,
               'energy_flag':  XrfSpectrum.energy_flag,
               'fwhm_offset':  XrfSpectrum.fwhm_offset,
               'fwhm_slope':   XrfSpectrum.fwhm_slope,
               'fwhm_flag':    XrfSpectrum.fwhm_flag,
               'chi_exp':      XrfSpectrum.chi_exp,
               'max_eval':     XrfSpectrum.max_eval,
               'max_iter':     XrfSpectrum.max_iter,
               'tolerance':    XrfSpectrum.tolerance}

    bgr_par = {'bottom_width':      Background.bottom_width,
               'bottom_width_flag': Background.bottom_width_flag,
               'top_width':         Background.top_width,
               'top_width_flag':    Background.top_width_flag,
               'exponent':          Background.exponent,
               'tangent':           Background.tangent,
               'compress':          Background.compress}

    peak_par = [peak_par[0],peak_par[1],etc...]
    peak_par[i] = {'label':       XrfPeak.label,  
                   'ignore':      XrfPeak.ignore,
                   'energy':      XrfPeak.energy,
                   'energy_flag': XrfPeak.energy_flag,
                   'fwhm':        XrfPeak.fwhm,
                   'fwhm_flag':   XrfPeak.fwhm_flag,
                   'ampl':        XrfPeak.ampl,
                   'ampl_factor': XrfPeak.ampl_factor,
                   'max_sigma':   XrfPeak.max_sigma,
                   'area':        XrfPeak.area}
        
"""
##############################################################################

import numpy as Num
import types
import time
import XrfPeaks
import XrfBgr
from   XrayData import xrf_lookup
from   Mca import McaCalib as calib
from   Mca import MedFile_CARS
from   Mca import MedFile_EMSA


##############################################################################
def read_xrf_files(prefix='xrf',start=0,end=100,nfmt=3,
                   bad_mca_idx=[],total=True,align=True,tau=None,
                   det_idx=0,emin=-1.0,emax=-1.0,fmt='CARS',xrf_params={}):
    """
    returns a list of xrf objects
    """
    if fmt == 'CARS':
        rd = MedFile_CARS.read_med_files
    elif fmt == 'EMSA':
        rd = MedFile_EMSA.read_med_files
    else:
        return
    med = rd(prefix,start=start,end=end,nfmt=nfmt,
             bad_mca_idx=bad_mca_idx,total=total,
             align=align,correct=correct,tau=tau)
    # generate xrf from med 
    xrf = Med2Xrf(med,xrf_params=xrf_params,
                  det_idx=det_idx,emin=emin,emax=emax)
    return xrf

##############################################################################
def read_xrf(file='',bad_mca_idx=[],total=True,align=True,tau=None,
             det_idx=0,emin=-1.0,emax=-1.0,fmt='CARS',xrf_params={}):
    """
    Read detector files
    >>m = xrf.read(file="file_name",bad_mca_idx=[],total=True,align=True,tau=None)

    """
    if fmt == 'CARS':
        rd = MedFile_CARS.read_med
    elif fmt == 'EMSA':
        rd = MedFile_EMSA.read_med
    else:
        return

    if type(file) == types.StringType:
        med = rd(file=file, bad_mca_idx=bad_mca_idx,
                 total=total, align=align, tau=tau)
    elif type(file) == types.ListType:
        med = []
        for f in file:
            tmp = rd(file=f, bad_mca_idx=bad_mca_idx,
                     total=total, align=align, tau=tau)
            med.append(tmp)
    else:
        return None

    # generate xrf from med 
    xrf = Med2Xrf(med,xrf_params=xrf_params,
                  det_idx=det_idx,emin=emin,emax=emax)
    return xrf

##############################################################################
def Med2Xrf(med,xrf_params={},det_idx=0,emin=-1.,emax=-1.):
    """
    Given an med object (or list of med objects) return a
    new xrf object (or list of xrf objects).
    
    Note xrf_params should have the same format as returned
    by Xrf.get_params
    """
    if type(med) == types.ListType:
        xrf = []
        for m in med:
            tmp = _med2xrf(med=med,xrf_params=xrf_params,
                           det_idx=det_idx,emin=emin,emax=emax)
            xrf.append(tmp)
    else:
        xrf = _med2xrf(med=med,xrf_params=xrf_params,
                       det_idx=det_idx,emin=emin,emax=emax)
    return xrf

def _med2xrf(med,xrf_params={},det_idx=0,emin=-1.,emax=-1.):
    # XRF just fits one data/energy array    
    if med.total == True:
        det_idx=0
    if (det_idx < 0) or (det_idx > med.n_detectors -1):
        det_idx=0
    
    # med will always ret data, energy and cal as arrays
    data = med.get_data()[det_idx]
    en   = med.get_energy()[det_idx]
    cal  = med.get_calib_params()[det_idx]
    
    # if energy range specified, truncate
    if ((emin + emax) > 0.) and (emin != emax):
        idx  = calib.energy_idx(en,emin=emin,emax=emax)
        data = data[idx]
        en   = en[idx]
    
    # Make sure we have energy calib params.
    # use those passed in if given, otherwise
    # use those from med.  
    if not xrf_params.has_key('fit'):
        xrf_params['fit'] = {}

    keys = xrf_params['fit'].keys()
    if 'energy_offset' not in keys:
        xrf_params['fit']['energy_offset'] = cal['offset']
    if 'energy_slope' not in keys:
        xrf_params['fit']['energy_slope'] = cal['slope']

    # convert energy back to chans
    chans = calib.energy_to_channel(en,offset=xrf_params['fit']['energy_offset'],
                                       slope=xrf_params['fit']['energy_slope'])
    
    # Return an XrfSpectrum object
    xrf =  XRF(data=data,chans=chans,params=xrf_params)

    return xrf

#################################################################################
def xrf_auto_fit(xrf_list,xrf_params={},use_prev_fit=False,fit_init=-1):
    """
    given a list of xrf objects fit them all
    """
    if fit_init >=0:
        xrf[fit_int].fit()
        params = xrf[fit_init].get_params()
        init = True
    elif len(xrf_params) > 0:
        params = xrf_params
        init = True
    else:
        init = False
        
    for j in range(len(xrf)):
        if (j > 0) and (use_prev_fit == True):
            params = xrf[j-1].get_params()
            xrf[j].init(params)
        else:
            if init:
                xrf[j].init(params)
            xrf[j].fit()
    return

#################################################################################
def xrf_fit_scan_files(prefix='xrf',start=0,end=100,nfmt=3,
                       bad_mca_idx=[],total=True,align=True,tau=None,
                       det_idx=0,emin=-1.0,emax=-1.0,fmt='CARS',xrf_params={},
                       use_prev_fit=False,fit_init=-1):
    """
    given an med object and set of xrf params - create xrf object and fit it

    work on this - logic etc... make sure does what needed...

    """
    xrf = read_xrf_files(prefix=prefix,start=start,end=end,nfmt=nfmt,
                         bad_mca_idx=bad_mca_idx,total=total,align=align,tau=tau,
                         det_idx=det_idx,emin=emin,emax=emax,fmt=fmt,
                         xrf_params=xrf_params)
    
    xrf_auto_fit(xrf,xrf_params=xrf_params,use_prev_fit=use_prev_fit,fit_init=fit_init)
    return xrf

#################################################################################
#################################################################################
class XRF(XrfPeaks.XrfSpectrum):
    """
    XRF parameters
    parameters = {'fit':fit_par,'bgr':bgr_par,'pk':peak_par}

    fit_par = {'energy_offset':XrfSpectrum.energy_offset,
               'energy_slope': XrfSpectrum.energy_slope,
               'energy_flag':  XrfSpectrum.energy_flag,
               'fwhm_offset':  XrfSpectrum.fwhm_offset,
               'fwhm_slope':   XrfSpectrum.fwhm_slope,
               'fwhm_flag':    XrfSpectrum.fwhm_flag,
               'chi_exp':      XrfSpectrum.chi_exp,
               'max_eval':     XrfSpectrum.max_eval,
               'max_iter':     XrfSpectrum.max_iter,
               'tolerance':    XrfSpectrum.tolerance}

    bgr_par = {'bottom_width':      Background.bottom_width,
               'bottom_width_flag': Background.bottom_width_flag,
               'top_width':         Background.top_width,
               'top_width_flag':    Background.top_width_flag,
               'exponent':          Background.exponent,
               'tangent':           Background.tangent,
               'compress':          Background.compress}

    peak_par = [peak_par[0],peak_par[1],etc...]
    peak_par[i] = {'label':       XrfPeak.label,  
                   'ignore':      XrfPeak.ignore,
                   'energy':      XrfPeak.energy,
                   'energy_flag': XrfPeak.energy_flag,
                   'fwhm':        XrfPeak.fwhm,
                   'fwhm_flag':   XrfPeak.fwhm_flag,
                   'ampl':        XrfPeak.ampl,
                   'ampl_factor': XrfPeak.ampl_factor,
                   'max_sigma':   XrfPeak.max_sigma,
                   'area':        XrfPeak.area}
    """
    #########################################################################
    def init_bgr(self,params={}):
        """
        Reset bgr parameters. 
            
        Inputs:
            params = {bottom_width=4.,top_width=0.,
                      exponent=2, tangent=0, compress=4}
        see XrfBgr.py for more details            
        """
        self.bgr = XrfBgr.Background(**params)

        return
    
    #########################################################################
    def set_bgr(self,params={}):
        """
        Set background parameters
        See Bgr module 
        """
        if self.bgr:
            self.bgr.init(params=params)
        return
    
    #########################################################################
    def init_lines(self,lines=[],guess=True,**kws):
        """
        Reset all peaks based on line identifiers.
        All peak parameters are defaults.
        
        lines:
            A list of xrf lines to be fit. The elements of the list
            should be strings such as ['Mn ka', 'Fe ka', 'Fe kb'] and/or
            a list of numbers eg. ['Mn ka', 'Fe ka', 7.12]
            See xrf_lookup.py for more details.

        guess:
            If guess is true the peak amplitude and fwhm will be estimated
            
        """
        self.peaks = []
        for line in lines:
            self.init_peak(line,guess=guess,**kws)

        return

    #######################################################################
    def init_peak(self,line,guess=True,**kws):
        """
        Create a new peak given an xrf line

        Inputs:
        
        line:
            Can be either a string representation of the line e.g. 'Mn ka',
            'Fe ka', 'Fe kb' (see xrf_lookup.py) or an energy (as float).
            
        guess:
            If guess is true the peak amplitude and fwhm will be estimated

        kws:
            Passed to peak

        """
        # creat peak_params dict
        peak_params = {}
        peak_params.update(kws)

        # get energy from line or
        # convert line to energy
        if type(line)==types.StringType:
            line = line.strip()
            en = xrf_lookup.lookup_xrf_line(line)
        else:
            try:
                en = float(line)
            except:
                en = None
        if en == None:
            print "Error setting energy for line"
            return
        else:
            peak_params['energy'] = en


        # create label if not specified
        if peak_params.has_key('label') == False:
            label = str(line)
            label.strip()
            peak_params['label']=label
            
        self._init_peak(peak_params=peak_params,guess=guess)

        return
    
    #######################################################################
    def _init_peak(self,peak_params={},guess=True):
        """
        Create a new peak given an xrf line
        """
        try:
            label = peak_params['label']
        except:
            print "Peak params is required to have a string label"
            return
        
        # see if its a duplicate 
        dup = self._peak_idx(label)
        if dup > -1:
            self.peaks.pop(dup)

        # create a new XrfPeak
        self.peaks.append( XrfPeaks.XrfPeak(**peak_params) )
        idx = len(self.peaks) - 1
        self._initSinglePeak(idx=idx,guess=guess)

        # return pk_idx of new peak_param
        return (len(self.peaks) - 1)

    #########################################################################
    def _peak_idx(self,label):
        """ See if a peak label exists and ret idx"""
        if label == None:
            return -1

        for j in range(len(self.peaks)):
            if self.peaks[j].label == label:
                return j
        return -1        

    #########################################################################
    def set_peak(self,label,params={}):
        """
        Inputs:
        
        label:
            String label for the peak.

        peak_params:
            Dictionary of peak parameters
            
        """
        idx = self._peak_idx(label)
        if idx < 0: return

        self._set_peak(idx,params=params)

    #######################################################################
    def _set_peak(self,idx=0,params={}):
        """
        Set peak parameters for the specified peak
        pk_idx:
            Integer index of the peak (see _peak_idx).  
        """
        if not idx in range(len(self.peaks)): return
        self.peaks[idx].init(params=params)
        
        return

    #########################################################################
    def get_peaks(self):
        peak_results = []
        for peak in self.peaks:
            peak_results.append(peak.get_params())
        return peak_results
    
    
##############################################################################
##############################################################################
if __name__ == "__main__":
    import pylab
    xrf = read_xrf(file='test.xrf',bad_mca_idx=[0,2,13],emin=4.,emax=9.)
    pylab.plot(xrf.get_energy(),xrf.get_data(),'-k')
    #pylab.show()
    #
    xrf.init_lines(['Fe ka',7.12])
    xrf.init_bgr()
    #
    xrf.calc()
    pylab.plot(xrf.get_energy(),xrf.predicted,'-r')
    #
    xrf.fit()
    en = xrf.get_energy()
    pylab.plot(en,xrf.predicted,'-b')
    #
    #for peak in xrf.peaks:
    #    pylab.plot(en,1+peak.calc(en),'.')
    #
    cnts = xrf.calc_peaks()
    for j in range(len(xrf.peaks)):
        pylab.plot(en,cnts[j],'.')
    #
    print xrf
    pylab.semilogy()
    pylab.show()
    
