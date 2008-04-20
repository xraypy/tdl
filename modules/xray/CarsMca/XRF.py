########################################################################
# XRF class for data handling and fitting
# T. Trainor
########################################################################

import numpy as Num
import types
import time
import copy
import Med
import Mca
import fitPeaks
import fitBgr
import xrf_lookup
import CarsMcaFile
import EmsaFile

"""
Conventions:
 - If total == true then data array etc has only a single entry --> the sum.
  otherwise the data array has an entry for all mca's (bad det's are just zeros)

- therefore det indexing always starts at zero (which referes to either the first
  mca, or to the sum depending on the total flag)

- self.ndet refers to the length of the data array, which is either 1 if just sum, or
  equal to the number of mca's

"""

##############################################################################
def read_xrf_file(file=None,bad_mca_idx=[],total=True,align=True,tau=[],fmt='CARS'):
    """
    Read detector files
    >>m = xrf.read("file_name",bad_mca_idx=[],total=True,align=True,tau=None)

    Returns an xrf object.  This function is appropriate for reading single and
    multi-element detectors.  (We always assume that the detector may be a
    multi-element detector)

    Inputs:
        file:
            File name
            
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

    """

    if file:
        if fmt == 'CARS':
            med = CarsMcaFile.read_med(file=file)
        elif fmt == 'EMSA':
            med = EmsaFile.read_med(file=file)
    else:
        med = Med.Med()
    xrf = XRF(med=med,bad_mca_idx=bad_mca_idx,total=total,align=align,tau=tau)
    return xrf


#####################################################
class XRF:
    def __init__(self,med=None,bad_mca_idx=[],total=True,align=True,tau=[]):

        self.med         = med

        # params, set these in init_data ...        
        self.ndet        = 0
        self.total       = True
        self.align       = True
        self.correct     = True
        self.bad_mca_idx = []
        self.tau         = []
        self.emin        = 0.0
        self.emax        = 0.0

        # lists for corrected/processed data
        #self.sum         = []
        self.data        = []
        self.bgr         = []
        self.predicted   = []
        
        # lists for fitting parameters
        self.bgr_params  = []
        self.fit_params  = []
        self.peak_params = []

        # update the data arrays
        self.init_data(bad_mca_idx=bad_mca_idx,total=total,align=align,
                       correct=self.correct,tau=tau,init_params=True)

        return

    #########################################################################
    # needs work -- should provide summary of fitting parameters, nicer fmt
    def __repr__(self):
        if self.med == None: return "No Data"
        lout = 'XRF Data Object:\n'
        lout = lout + '  Name = %s\n' % self.med.name
        lout = lout + '  Num detectors = %i\n' % self.med.n_detectors
        lout = lout + '  Bad Detectors = %s\n' % str(self.bad_mca_idx)
        lout = lout + '  Detector Taus = %s\n' % str(self.tau)
        lout = lout + '  Align = %s\n' % str(self.align)
        lout = lout + '  Total = %s\n' % str(self.total)
        lout = lout + '  Correct = %s\n' % str(self.correct)
        ##
        for pk in self.peak_params:
            for p in pk:
                lout = lout + '  Peak=%s, Energy=%s,FWHM=%s,Amp=%s,Area=%s,Ignore=%s\n'  % \
                       (p.label, str(p.energy),str(p.fwhm),str(p.ampl),str(p.area),str(p.ignore))
        return lout

    def show_rois(self):
        if self.med == None: return "No Data"
        lout = self.med.__repr__()
        return lout

    #########################################################################
    def init_data(self,bad_mca_idx=None,total=None,align=None,correct=None,
                  tau=None,init_params=True):
        """
        Initialize or re-initialize the data.
        Note during re-initialization only passed parameters will modified, ie
        the orginal values set during file read will only be changed if set explicitly.

        Inputs:

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

        correct:
            Set this keyword to apply deadtime corrections to data

        tau:
            List of deadtime factors for each detector. Dimesion should
            be the same as the total number of detectors, or [] for no tau values.
            If tau is None, previously passed tau is used...  
              If a detector tau > 0 this will be used in the correction factor calculation 
              If a detector tau = 0 then we assume ocr = icr in the correction factor calculation, ie only lt correction
              If a detector tau < 0 (or None):
                if input_counts > 0 this will be used for icr in the factor calculation
                if input_counts <= 0 we assume ocr = icr in the correction factor calculation, ie only lt correction
            Note:
             ocr = icr * exp(-icr*tau)
             cor = (icr/ocr)*(rt/lt)
             max_icr = 1/tau
             max_ocr = max_icr*exp(-1)

        init_parameters:
            Set to True to reset all fitting parameters

        """

        if self.med == None: return

        # update parameters
        if bad_mca_idx is not None:
            self.bad_mca_idx = bad_mca_idx
        if total is not None:
            self.total = total
        if align is not None:
            self.align = align
        if correct is not None:
            self.correct = correct

        # update deadtime correction factors
        if tau is not None:
            self.tau = tau
            self.med.update_correction(tau=self.tau)
        else:
            self.med.update_correction(tau=None)

        # recalc rois
        self.med.update_rois(correct=correct)
 
        # clear arrays
        self.data = []

        #get data 
        self.data = self.med.get_data(bad_mca_idx=self.bad_mca_idx,
                                      total=self.total,
                                      align=self.align,
                                      correct=self.correct)

        self.ndet = len(self.data)
        
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # reinit the bgr and peak parameters
        # note this blows away everything!
        # we should be more clever about this
        # and only blow stuff away if ndet changes...
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if init_params == True:
            self.init_params()

    #########################################################################
    def init_params(self,lines = [], bottom_width=4.,top_width=0.,
                    exponent=2, tangent=0, compress=4):
        """
        Reset all peak and bgr parameters to defaults. Note xrf lines and background
        parameters are copied to all detectors (if total = False).  See the following
        XRF methods for more details/control over fitting parameters
            init_peak
            set_peak
            set_bgr
            
        Inputs:

        lines:
            A list of xrf lines to be fit. The elements of the list should be strings
            such as ['Mn ka', 'Fe ka', 'Fe kb'].  See xrf_lookup.py for more details.

         bottom_width:
            Specifies the width of the polynomials which are concave downward.
            The bottom_width is the full width in energy units at which the
            magnitude of the polynomial is 100 counts. The default is 4.

        top_width:
            Specifies the width of the polynomials which are concave upward.
            The top_width is the full width in energy units at which the
            magnitude of the polynomial is 100 counts. The default is 0, which
            means that concave upward polynomials are not used.

         exponent:
            Specifies the power of polynomial which is used. The power must be
            an integer. The default is 2, i.e. parabolas. Higher exponents,
            for example EXPONENT=4, results in polynomials with flatter tops
            and steeper sides, which can better fit spectra with steeply
            sloping backgrounds.

         tangent:
            Specifies that the polynomials are to be tangent to the slope of the
            spectrum. The default is vertical polynomials. This option works
            best on steeply sloping spectra. It has trouble in spectra with
            big peaks because the polynomials are very tilted up inside the
            peaks.

         compress:
            Compression factor to apply before fitting the background.
            Default=4, which means, for example, that a 2048 channel spectrum
            will be rebinned to 512 channels before fitting.
            The compression is done on a temporary copy of the input spectrum,
            so the input spectrum itself is unchanged.
            The algorithm works best if the spectrum is compressed before it
            is fitted. There are two reasons for this. First, the background
            is constrained to never be larger than the data itself. If the
            spectrum has negative noise spikes they will cause the fit to be
            too low. Compression will smooth out such noise spikes.
            Second, the algorithm requires about 3*N^2 operations, so the time
            required grows rapidly with the size of the input spectrum. On a
            200 MHz Pentium it takes about 3 seconds to fit a 2048 channel
            spectrum with COMPRESS=1 (no compression), but only 0.2 seconds
            with COMPRESS=4 (the default).

            see fitBgr.py for more details            
        """

        self.fit_params  = []
        self.bgr_params  = []
        self.peak_params = []
        self.bgr         = []
        self.predicted   = []
        for i in range(self.ndet):
            self.fit_params.append([])
            self.bgr_params.append([])
            self.peak_params.append([])
            self.bgr.append([])
            self.predicted.append([])

        for i in range(self.ndet):
            self.fit_params[i] = fitPeaks.McaFit()
            calib = self.get_calibration(i)
            self.bgr_params[i] = fitBgr.McaBackground(slope=calib.slope,
                                                      exponent=exponent,
                                                      top_width=top_width,
                                                      bottom_width=bottom_width,
                                                      tangent=tangent,
                                                      compress=compress)
            for line in lines:
                self.init_peak(line,det_idx=i)

        return

    #########################################################################
    def get_params(self):
        """ Returns the fit_params, bgr_params and peak_params""" 
        return (self.fit_params, self.bgr_params, self.peak_params)

    #def set_params(): etc...

    #########################################################################
    def get_data(self):
        """
        Returns the data array.
        The data is always a list of dim [ndetectors, nchannels]
        Note use the method init_data to change how the data is processed
        """
        return self.data

    #########################################################################
    def get_energy(self):
        """
        Returns the energy array
        The energy is always a list of dim [ndetectors, nchannels]
        Note use the method init_data to change how the energy is handled
        """
        en = self.ndet*[[]]
        if self.total:
            calib    = self.get_calibration(0)
            channels = Num.arange(len(self.data[0]))
            en[0]    = calib.channel_to_energy(channels)
        elif self.align:
            calib    = self.get_calibration(0)
            channels = Num.arange(len(self.data[0]))
            ref_en   = calib.channel_to_energy(channels)
            for i in range(self.ndet):
                en[i] = ref_en
        else:
            for i in range(self.ndet):
                calib    = self.get_calibration(i)
                channels = Num.arange(len(self.data[i]))
                en[i]    = calib.channel_to_energy(channels)
        return en

    #########################################################################
    def get_calibration(self,det_idx):
        "get an mca's calibration"
        # note if you did a total without aligning
        # there is no way to return a "correct" calibration
        # so it self.total == True will always return the
        # calibration of the first good mca
        idx = self._get_calibration_idx(det_idx)
        return self.med.mcas[idx].calibration

    #########################################################################
    def _get_calibration_idx(self,det_idx):
        "get an mca's calibration index"
        # note if you did a total without aligning
        # there is no way to return a "correct" calibration
        # so if self.total == True will always return the
        # calibration of the first good mca
        if self.total or self.align:
            idx = self.med.get_align_idx(self.bad_mca_idx)
            return idx
        else:
            return det_idx 

    #######################################################################
    def init_peak(self,line,label=None,det_idx=[]):
        """
        Create new peak params given an xrf line

        Inputs:
        
        line:
            Can be either a string representation of the line e.g. 'Mn ka',
            'Fe ka', 'Fe kb' (see xrf_lookup.py) or an energy (as float).
            Note if line == None it will clear all the peaks
        
        label:
            String label for the peak.  If label == None it will use a
            string representation of line
        
        det_idx:
            The index of the detector to add the peak to.  Should be an integer
            in the range 0 to ndetectors, or a list of integers.  If this is None
            or an empty list (default) the peak_param is created for all detectors
            
        """
        if det_idx == None or len(det_idx) == 0:
            det_idx = range(self.ndet)
        if type(det_idx) == types.ListType:
            for det in det_idx:
                self._init_peak(line,label=label, det_idx=det)
        else:
            try:
                det = int(det_idx)
                self._init_peak(line,label=label, det_idx=det)
            except:
                print "Error setting line for detector: ", det_idx
        return
    
    #######################################################################
    def _init_peak(self,line,label=None,det_idx=0):
        """
        Create a new peak_param for the specified detector given an xrf line
        and label
        """
        if not det_idx in range(self.ndet): return

        # if line == None, blow away all peaks for this detector
        if line == None:
            self.peak_params[det_idx] = []
            return

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
            return 

        # get label
        if label == None:
            label = str(line)
            label.strip()
        
        # see if its a duplicate 
        dup = self._peak_idx(line, det_idx=det_idx)
        if dup > -1:
            self.peak_params[det_idx].pop(dup)

        # create a new McaPeak structure
        tmp = fitPeaks.McaPeak()
        tmp.label          = label
        tmp.initial_energy = en
        tmp.energy         = en
        self.peak_params[det_idx].append(tmp)

        # return pk_idx of new peak_param
        return (len(self.peak_params[det_idx]) - 1)

    #########################################################################
    def _peak_idx(self,label,det_idx=0):
        """ See if a peak label exists and ret idx"""
        if not det_idx in range(self.ndet):
            return -1
        if label == None:
            return -1

        for j in range(len(self.peak_params[det_idx])):
            if self.peak_params[det_idx][j].label == label:
                return j
        return -1        

    #########################################################################
    def set_peak(self,label,energy=None,ampl=None,fwhm=None,energy_flag=None,
                 fwhm_flag=None,ampl_factor=None,ignore=False,det_idx=[]):
        """
        Inputs:
        
        label:
            String label for the peak. Note if the peak does not already
            exist it this function will create it.  

        (Note: use None for the below peak paramters to keep current value unchanged)
        
        energy:
            Peak energy in keV

        ampl:
            Peak amplitude

        fwhm:
            Peak full width at half max

        energy_flag:
            Flag for fitting energy.  0 = fix the peak energy, 1 = optimize the peak
            energy.  Note peak energy can be fixed, while optimizing the spectra energy
            calibration.  See the energy_flag in fit method.  

        fwhm_flag:
            Flag for fitting peak fwhm. 0 = Fix fwhm to global curve. 1 = optimize
            the peak fwhm, 2 = fix the fwhm to input value.

        ampl_factor:
            Fixed amplitude ratio to previous peak
            0.  = Optimize amplitude of this peak
            >0. = Fix amplitude to this value relative
                  to amplitude of previous free peak
            -1.0 = Fix amplitude at 0.0

        ignore:
            Flag to ignore peak (True/False), default = False

        det_idx:
            The index of the detector to add the peak parameters to.  Should be an integer
            in the range 0 to ndetectors, or a list of integers.  If this is None
            or an empty list (default) the parameters are applied to all detectors

        """
        if det_idx == None or len(det_idx) == 0:
            det_idx = range(self.ndet)

        if type(det_idx) == types.ListType:
            for det in det_idx:
                pk_idx = _peak_idx(label,det_idx=det)
                if pk_idx < 0:
                    pk_idx = self._init_peak(line,label=label,det_idx=det)
                self._set_peak(pk_idx=pk_idx,det_idx=det_idx,energy=energy,ampl=ampl, 
                               fwhm=fwhm,energy_flag=energy_flag,fwhm_flag=fwhm_flag,
                               ampl_factor=ampl_factor,ignore=ignore)
        else:
            try:
                det = int(det_idx)
                pk_idx = _peak_idx(label,det_idx=det)
                if pk_idx < 0:
                    pk_idx = self._init_peak(line,label=label,det_idx=det)
                self._set_peak(pk_idx=pk_idx,det_idx=det_idx,energy=energy,ampl=ampl, 
                               fwhm=fwhm,energy_flag=energy_flag,fwhm_flag=fwhm_flag,
                               ampl_factor=ampl_factor,ignore=ignore)
            except:
                print "Error setting parameters for detector: ", det_idx


    #######################################################################
    def _set_peak(self,pk_idx=0,det_idx=0,energy=None,ampl=None,fwhm=None,
                  energy_flag=None,fwhm_flag=None,ampl_factor=None,ignore=False):
        """
        Set peak parameters for the specified detector and peak
        
        pk_idx:
            Integer index of the peak (see _peak_idx).  If its a sting, assume
            its the peak label and try to find the corresponding peak index

        det_idx:
            Integer index of the detector
        """
        if not det_idx in range(self.ndet): return

        if type(pk_idx) == types.StringType:
            pk_idx = self._peak_idx(pk_idx.strip(),det_idx=det_idx)
        
        if not pk_idx in range(len(self.peak_params[det_idx])): return

        if energy is not None:
            self.peak_params[det_idx][pk_idx].initial_energy = energy
            self.peak_params[det_idx][pk_idx].energy = energy

        if ampl is not None:
            self.peak_params[det_idx][pk_idx].initial_ampl = ampl
            self.peak_params[det_idx][pk_idx].ampl = ampl
            
        if fwhm is not None:
            self.peak_params[det_idx][pk_idx].initial_fwhm = fwhm
            self.peak_params[det_idx][pk_idx].fwhm = fwhm

        if energy_flag is not None:
            self.peak_params[det_idx][pk_idx].energy_flag = energy_flag
            
        if fwhm_flag is not None:
            self.peak_params[det_idx][pk_idx].fwhm_flag = fwhm_flag

        if ampl_factor is not None:
            self.peak_params[det_idx][pk_idx].ampl_factor = ampl_factor

        self.peak_params[det_idx][pk_idx].ignore = ignore

        return

    #########################################################################
    def set_bgr(self,slope=None,exponent=None,top_width=None,bottom_width=None,
                tangent=None,compress=None,det_idx=[]):
        """
        Set background parameters for a detector
        
        Inputs:

        (Note: use None for the below background paramters to keep current value unchanged)

        slope:
            Slope for channel to energy calculation
        
        bottom_width:
            Specifies the width of the polynomials which are concave downward.
            The bottom_width is the full width in energy units at which the
            magnitude of the polynomial is 100 counts. The default is 4.

        top_width:
            Specifies the width of the polynomials which are concave upward.
            The top_width is the full width in energy units at which the
            magnitude of the polynomial is 100 counts. The default is 0, which
            means that concave upward polynomials are not used.

         exponent:
            Specifies the power of polynomial which is used. The power must be
            an integer. The default is 2, i.e. parabolas. Higher exponents,
            for example EXPONENT=4, results in polynomials with flatter tops
            and steeper sides, which can better fit spectra with steeply
            sloping backgrounds.

         tangent:
            Specifies that the polynomials are to be tangent to the slope of the
            spectrum. The default is vertical polynomials. This option works
            best on steeply sloping spectra. It has trouble in spectra with
            big peaks because the polynomials are very tilted up inside the
            peaks.

         compress:
            Compression factor to apply before fitting the background.
            Default=4, which means, for example, that a 2048 channel spectrum
            will be rebinned to 512 channels before fitting.
            The compression is done on a temporary copy of the input spectrum,
            so the input spectrum itself is unchanged.
            The algorithm works best if the spectrum is compressed before it
            is fitted. There are two reasons for this. First, the background
            is constrained to never be larger than the data itself. If the
            spectrum has negative noise spikes they will cause the fit to be
            too low. Compression will smooth out such noise spikes.
            Second, the algorithm requires about 3*N^2 operations, so the time
            required grows rapidly with the size of the input spectrum. On a
            200 MHz Pentium it takes about 3 seconds to fit a 2048 channel
            spectrum with COMPRESS=1 (no compression), but only 0.2 seconds
            with COMPRESS=4 (the default).

            see fitBgr.py for more details

        det_idx:
            The index of the detector to add the peak parameters to.  Should be an integer
            in the range 0 to ndetectors, or a list of integers.  If this is None
            or an empty list (default) the parameters are applied to all detectors
  
        """
        if det_idx == None or len(det_idx) == 0:
            det_idx = range(self.ndet)
        if type(det_idx) == types.ListType:
            for det in det_idx:
                self._set_bgr(det_idx=det,slope=slope,exponent=exponent,top_width=top_width,
                              bottom_width=bottom_width, tangent=tangent,compress=compress)
        else:
            try:
                det = int(det_idx)
                self._set_bgr(det_idx=det,slope=slope,exponent=exponent,top_width=top_width,
                              bottom_width=bottom_width, tangent=tangent,compress=compress)
            except:
                print "Error setting background for detector: ", det_idx
        return

    #########################################################################
    def _set_bgr(self,det_idx=0,slope=None,exponent=None,top_width=None, 
                 bottom_width=None,tangent=None,compress=None,):
        """
        Set bgr parameters for a detector
        """
        if not det_idx in range(self.ndet): return
        if slope        is not None: self.bgr_params[det_idx].slope        = slope
        if exponent     is not None: self.bgr_params[det_idx].exponent     = exponent
        if top_width    is not None: self.bgr_params[det_idx].top_width    = top_width
        if bottom_width is not None: self.bgr_params[det_idx].bottom_width = bottom_width
        if tangent      is not None: self.bgr_params[det_idx].tangent      = tangent
        if compress     is not None: self.bgr_params[det_idx].compress     = compress

        return

    #########################################################################
    def fit(self,fwhm_flag=1,energy_flag=1,chi_exp=0.0,fit_bgr=True):
        """
        Fit the data.
        Inputs:

        fwhm_flag:
            0 = Fix global FWHM coefficients
            1 = Optimize global FWHM coefficients (Default)

        energy_flag:
            0 = Fix energy calibration coefficients
            1 = Optimize energy calibration coefficients (Default)

        chi_exp:
            Exponent of chi (Default = 0.0). The fit function assumes that:
                sigma[i] = y_obs[i] ** chi_exponent
            e.g. that the standard deviation in each channel is equal to the counts
            in the channel to some power. For photon counting spectra where Poisson
            statistics apply chi_exponent=0.5. Setting chi_exponent=0. will set all
            of the sigma[i] values to 1., and the fit
            would then be minimizing the sum of the squares of the residuals. This
            should tend to result in a better fit for the large peaks in a spectrum
            and a poorer fit for the smaller peaks. Setting chi_exponent=1.0 will
            result in a minimization of the sum of the squares of the relative error
            in each channel. This should tend to weight the fit more strongly toward
            the small peaks.

        fit_bgr: (True/False)
            Flag to indicate if the background should be fit (and removed)
            before fitting peaks (default = True)
        """
        for i in range(self.ndet):
            #print i, fit_bgr
            if fit_bgr: self._fit_bgr(i)
            self._fit_peaks(i,fwhm_flag=fwhm_flag,energy_flag=energy_flag,chi_exp=chi_exp)
        return 

    #########################################################################
    def fit_bgr(self):
        """
        Fit the background
        """
        for i in range(self.ndet):
            self._fit_bgr(i)
        return

    #########################################################################
    def _fit_bgr(self,det_idx):

        if not det_idx in range(self.ndet): return
        if det_idx in self.bad_mca_idx:
            self.bgr[det_idx] = Num.zeros(len(self.data[det_idx]))
        else:
            calib = self.get_calibration(det_idx)
            self.bgr_params[det_idx].slope = calib.slope
            self.bgr[det_idx] = fitBgr.fit_background(self.data[det_idx],
                                                      self.bgr_params[det_idx])

        return

    #########################################################################
    def fit_peaks(self,fwhm_flag=1,energy_flag=1,chi_exp=0.0):

        for i in range(self.ndet):
            self._fit_peaks(i,fwhm_flag=fwhm_flag,energy_flag=energy_flag,chi_exp=chi_exp)
        return 

    def _fit_peaks(self,det_idx,fwhm_flag=1,energy_flag=1,chi_exp=0.0):

        if not det_idx in range(self.ndet): return
        if self.total == False:
            if det_idx in self.bad_mca_idx:
                print "skip det ", det_idx
                self.predicted[det_idx] = Num.zeros(len(self.data[det_idx]))
                return
        observed    = self.data[det_idx]
        background  = self.bgr[det_idx]
        
        # in case background was not subtracted
        if len(background) == 0:
            background = Num.zeros(len(observed))
        calib_idx    = self._get_calibration_idx(det_idx)
        peaks        = self.peak_params[det_idx]
        fit          = self.fit_params[det_idx]
        #print "peak params for ", det_idx
        #print "length of peaks", len(peaks)

        if len(peaks) == 0:
            #print "det ", det_idx, " has no peak parameters"
            self.predicted[det_idx] = Num.zeros(len(self.data[det_idx]))
            return

        # Copy parameters to fit
        fit.npeaks                = len(peaks)
        fit.initial_energy_offset = self.med.mcas[calib_idx].calibration.offset
        fit.initial_energy_slope  = self.med.mcas[calib_idx].calibration.slope
        fit.nchans                = len(observed)
        fit.last_chan             = fit.nchans-1
        fit.fwhm_flag             = fwhm_flag
        fit.energy_flag           = energy_flag
        fit.chi_exp               = chi_exp
        
        t0 = time.time()
        [fit, peaks, fit_counts] = fitPeaks.fitPeaks(fit, peaks, observed - background)
        t1 = time.time()
        
        self.predicted[det_idx]   = fit_counts + background
        self.fit_params[det_idx]  = fit
        self.peak_params[det_idx] = peaks
        self.med.mcas[calib_idx].calibration.offset = fit.energy_offset
        self.med.mcas[calib_idx].calibration.slope  = fit.energy_slope

        return
    
    #########################################################################
    def calc_peaks(self,calc_bgr=True):

        for i in range(self.ndet):
            if calc_bgr: self._fit_bgr(i)
            self._calc_peaks(i)
        return 

    #########################################################################
    def _calc_peaks(self,det_idx):

        if not det_idx in range(self.ndet): return

        if self.total == False:
            if det_idx in self.bad_mca_idx:
                self.predicted[det_idx] = Num.zeros(len(self.data[det_idx]))
                return
        
        # Bgr
        background   = self.bgr[det_idx]
        if len(background) == 0:
            background = Num.zeros(len(self.data[det_idx]))
        
        calibration  = self.get_calibration(det_idx)
        peaks        = self.peak_params[det_idx]
        fit          = self.fit_params[det_idx]

        # Copy parameters to fit
        fit.npeaks        = len(peaks)
        fit.energy_offset = calibration.offset
        fit.energy_slope  = calibration.slope
        fit.nchans        = len(self.data[det_idx])
        fit.last_chan     = fit.nchans-1
        
        fit_counts = fitPeaks.predict_gaussian_spectrum(fit, peaks)
        self.predicted[det_idx] = fit_counts + background
        
        return

    #########################################################################
    def calc_pk(self,det_idx,pk_idx,add_bgr = False):
        "calc and ret single peak"
# work on this
        if not det_idx in range(self.ndet): return

        if type(pk_idx) == types.StringType:
            pk_idx = self._peak_idx(pk_idx)
            
        if not pk_idx in range(len(self.peak_params[det_idx])): return

        if self.total == False:
            if det_idx in self.bad_mca_idx:
                self.predicted[det_idx] = Num.zeros(len(self.data[det_idx]))
                return

        # Bgr
        if add_bgr: self._fit_bgr(det_idx)
        background   = self.bgr[det_idx]
        if len(background) == 0:
            background = Num.zeros(len(self.data[det_idx]))

        calibration  = self.get_calibration(det_idx)
        peaks        = [self.peak_params[det_idx][pk_idx]]
        fit          = self.fit_params[det_idx]

        # Copy parameters to fit
        fit.npeaks        = 1
        fit.energy_offset = calibration.offset
        fit.energy_slope  = calibration.slope
        fit.nchans        = len(self.data[det_idx])
        fit.last_chan     = fit.nchans-1
        
        fit_counts = fitPeaks.predict_gaussian_spectrum(fit, peaks)
        
        if add_bgr:
            return fit_counts + background
        else:
            return fit_counts

    #########################################################################
    def get_count_totals(self):
        lst  = []
        for mca in self.med.mcas:
            rt = mca.elapsed.real_time
            lt = mca.elapsed.live_time
            OCR = mca.elapsed.total_counts/mca.elapsed.live_time
            ICR = mca.elapsed.input_counts/mca.elapsed.live_time
            ICR_CALC = mca.elapsed.icr_calc
            COR = mca.elapsed.cor_factor
            lst.append({'rt':lt,'lt':lt,'OCR':OCR,'ICR':ICR,'ICR_CALC':ICR_CALC,'COR':COR})
        return lst
            
    #########################################################################
    def data_dict(self):
        dict = {'energy':[],'counts':[],'background':[],'predicted':[],
                'peak_areas':[],'peaks':[]}
        energy    = self.get_energy()
        counts    = self.get_data()
        bgr       = self.bgr
        predicted = self.predicted

        for mca_idx in range(self.ndet):
            dict['energy'].append(energy[mca_idx])
            dict['counts'].append(counts[mca_idx])
            dict['background'].append(bgr[mca_idx])
            dict['predicted'].append(predicted[mca_idx])

            # peak areas
            peak_areas = {}
            for peak in self.peak_params[mca_idx]:
                peak_areas[peak.label] = peak.area
            dict['peak_areas'].append(peak_areas)

            # get each peak function
            peaks = {}
            for peak_params in self.peak_params[mca_idx]:
                calc = fitPeaks.predict_gaussian_spectrum(self.fit_params[mca_idx],
                                                          [peak_params])
                peaks[peak_params.label] = calc
            dict['peaks'].append(peaks)

        return dict

    #########################################################################
    def get_peaks(self):

        results = []
        for mca_idx in range(self.ndet):
            peak_results = {}
            for peak in self.peak_params[mca_idx]:
                peak_results[peak.label] = {'energy':0.0,'ampl':0.0,'fwhm':0.0,'area':0.0}
                peak_results[peak.label]['energy'] = peak.energy
                peak_results[peak.label]['ampl'] = peak.ampl
                peak_results[peak.label]['fwhm'] = peak.fwhm
                peak_results[peak.label]['area'] = peak.area
            results.append(peak_results)
        return results

    #########################################################################
    def get_rois(self,background_width=1,correct=True):

        # if total we could calc rois for the total...
        if self.total:
            ret = {}
            if self.total:
                idx = self.med.get_align_idx(self.bad_mca_idx)
                rois = copy.copy(self.med.mcas[idx].rois)
                for roi in rois:
                    roi.update_counts(self.data[0], background_width=background_width)
                    ret[roi.label] = (roi.total, roi.net)
                return [ret]
        else:
            return self.med.get_roi_counts_lbl(background_width=background_width,
                                               correct=correct)

    #########################################################################
    def set_roi(self,label,lrn=[],mcas=[],units='keV'):

        if mcas == []: mcas = range(self.med.n_detectors)

        # delete rois
        for i in mcas:
            d = int(i)
            if isinstance(label, types.IntType):
                med.mcas[d].delete_roi(label)
            elif type(label) == types.StringType:
                idx =  med.mcas[d].find_roi_label(label=label)
                if idx > 0:
                    med.mcas[d].delete_roi(idx)
            else:
                print "roi needs to be string label or integer index"
                return

        if lrn == []:
            return

        # add rois
        if len(lrn) == 2:
            left = lrn[0]
            right = lrn[1]
            nbgr  = 1
        elif len(lrn) == 3:
            left = lrn[0]
            right = lrn[1]
            nbgr  = int(lrn[2])
        else:
            print "lrn needs either 2 or 3 vals"
            return

        roi = Mca.McaROI(units=units,left=left,right=right,bgd_width=nbgr,label=str(label))
        for i in mcas:
            d = int(i)
            med.mcas[d].add_roi(roi)

