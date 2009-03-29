#######################################################################
"""
T. Trainor (fftpt@uaf.edu)
Xrf class for data handling and fitting

Modifications:
--------------

"""
########################################################################
"""
Todo

- Improve documentation
- Test!

"""
########################################################################

import types
import numpy as Num

from   detector import mca_calib as calib
import xrf_peaks
import xrf_bgr
from   xtab import xrf_lookup

#################################################################################
class Xrf(xrf_peaks.XrfSpectrum):
    """
    Xrf subclasses XrfSpectrum and includes the following data:
    
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
        see xrf_bgr.py for more details            
        """
        self.bgr = xrf_bgr.Background(**params)

        return
    
    #########################################################################
    def set_bgr(self,params={}):
        """
        Set background parameters
        See Bgr module 
        """
        if self.bgr:
            self.bgr.init(params=params)
        else:
            self.init_bgr(params=params)
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
        # create peak_params dict
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
        self.peaks.append( xrf_peaks.XrfPeak(**peak_params) )
        idx = len(self.peaks) - 1
        self._initSinglePeak(idx=idx,guess=guess)

        # return pk_idx of new peak_param
        return (len(self.peaks) - 1)

    #########################################################################
    def _peak_idx(self,label):
        """
        See if a peak label exists and ret idx
        """
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
        """
        get peaks
        """
        peak_results = []
        for peak in self.peaks:
            peak_results.append(peak.get_params())
        return peak_results
