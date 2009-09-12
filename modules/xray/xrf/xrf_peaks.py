#######################################################################
"""
Mark Rivers, GSECARS
Fit a spectrum to a set of Gaussian peaks.  

Modifications:
--------------

 * Mark Rivers, October 21, 1998.
    This is the latest re-write of a routine which has a long history, begun
    at X-26 at the NSLS.  The original version was written in a program called
    SPCALC, and was then ported to IDL.
    These early versions used IMSL for the least-squares routine.  The port to
    CURVEFIT, so that no external software package is required, was done in
    1998.

 * Mark Rivers, Nov. 12, 1998.
    Significant re-write to use MPFITFUN in place of CURVEFIT

 * Mark Rivers, Feb. 1, 2001.
    Changed amplitude ratio calculation so that the AREA of the two peaks
    has the specified ratio, rather than the AMPLITUDE.  This is done by
    adjusting the constrained ratio by the relative peak widths.

 * Mark Rivers, Sept 18, 2002.
    Converted from IDL to Python.

 * Mark Rivers, Sept 25., 2002 
    Previously several fields in the peaks could be clobbered if a
    peak was outside the energy range of the spectrum.  Added new .ignore
    field to McaPeak to work around this problem and use that field here.

 * See http://cars9.uchicago.edu/software/python/index.html

 * TPT: update for tdl
    use numpy and new mpfit
 
"""
########################################################################
"""
Todo:
 - check fitting: bgr handled correctly for all options?
 - check _fit_peaks: is this correct way/place to apply weights
"""
#######################################################################
"""
Procedure:
    In general a Gaussian peak has 3 adjustable parameters: position 
    (or energy), sigma (or FWHM), and amplitude (or area).  For many
    applications, however, not all of these parameters should be
    adjustable during the fit.  For example, in XRF analysis the energy of
    the peaks is known, and should not be optimized.  However, the overall
    energy calibration coefficients for the entire spectrum, which relate
    channel number to energy, might well be optimized during the fit.
    Similarly, the FWHM of XRF peaks are not independent, but rather
    typically follow a predictable detector response function:
         FWHM = A + B*sqrt(energy)
    Finally, even the amplitude of an XRF peak might not be a free
    parameter, since, for example one might want to constrain the K-beta 
    peak to be a fixed fraction of the K-alpha.  Such constraints allow 
    one to fit overlapping K-alpha/K-beta peaks with much better accuracy.

    This procedure is designed to be very flexible in terms of which
    parameters are fixed and which ones are optimized.  The constraints are
    communicated via the Fit and Peaks structures.

    The energy of each channel is assumed to obey the relation: 
         energy = energy_offset + (channel * energy_slope)

    These parameters control the fit for peaks whose energy is fixed, 
    rather than being a fit parameter.
    If Fit.energy_flag is 1 then these energy calibration coefficients
    will be optimized during the fitting process. If it is 0 then these
    energy calibration coefficients are assumed to be correct and are not 
    optimized.  Not optimizing the energy calibration coefficients can 
    both speed up the fitting process and lead to more stable results when 
    fitting small peaks.  This function does a sanity check and will not
    optimize these energy calibration coefficients unless at least 2 peaks
    have their .energy_flag field set to 0, so that they use these global
    calibration coefficients.

    The FWHM of the peaks is assumed to obey the relation:
         fwhm = fwhm_offset + (fwhm_slope * sqrt(energy))
    These parameters control the fit for peaks whose FWHM is neither fixed 
    nor a fit parameter.
    If Fit.fwhm_flag is 1 then these coefficients will be optimized during
    the fitting process. If it is 0 then the specified coefficients are 
    assumed to be correct and are not optimized. Not optimizing the FWHM
    coeffcients can both speed up the fitting process and lead to more 
    stable results when fitting very small peaks. This function does a 
    sanity check and will not optimize these FWHM calibration coefficients 
    unless at least 2 peaks have their .fwhm_flag field set to 0, so that 
    they use these global calibration coefficients.

    This function also optimizes the following parameters:
         - The amplitudes of all peaks whose .ampl_factor field is 0
         - The energies of all peaks whose .energy_flag field is 1
         - The FWHM of all peaks whose .fwhm_flag field is 1

    The parameter which is the minimized during the fitting process is 
    chi**2, defined as:
                                                 2
        2            y_obs[i]    -     y_pred[i]
    chi  = sum (  ---------------------------- )
              i              sigma[i]

    where y_obs[i] is the observed counts in channel i, y_pred is the
    predicted counts in channel i, and sigma[i] is the standard deviation
    of y_obs[i].

    This function assumes that:

    sigma[i] = y_obs[i] ** chi_exponent

    e.g. that the standard deviation in each channel is equal to the counts
    in the channel to some power. For photon counting spectra where Poisson
    statistics apply chi_exponent=0.5. Setting chi_exponent=0.0 (default)
    will set all of the sigma[i] values to 1., and the fit
    would then be minimizing the sum of the squares of the residuals. This
    should tend to result in a better fit for the large peaks in a spectrum
    and a poorer fit for the smaller peaks. Setting chi_exponent=1.0 will
    result in a minimization of the sum of the squares of the relative error
    in each channel. This should tend to weight the fit more strongly toward
    the small peaks.

    If .ampl_factor for a peak is 0., then the amplitude of the peak is a 
    fit parameter. If the amplitude_factor is non-zero then the amplitude 
    of this peak is not a fit parameter, but rather is constrained to
    be equal to the amplitude of the last previous peak in the array which 
    had an amplitude factor of zero, times the amplitude_factor. This can 
    be used, for instance, fit K-alpha and K-beta x-ray lines when the 
    alpha/beta ratio is known, and one wants to add this known constraint 
    to the fitting process.
    For example:
         peaks = replicate({mca_peak}, 3)
         # Fe Ka is the "reference" peak
         peaks[0].initial_energy=6.40 & peaks[0].ampl_factor=0.0 
         # Si-Ka escape peak is 3% of Fe Ka at 4.66 keV
         peaks[1].initial_energy=4.66 & peaks[1].ampl_factor=0.03
         # Fe-Kb is 23% of Fe Ka
         peaks[2].initial_energy=7.06 & peaks[2].ampl_factor=0.23
    In this example the amplitude of the Fe-Ka peak will be fitted, but the
    amplitudes of the escape peak and the Fe-Kb peak are constrained to
    be fixed fractions of the Fe-Ka peak.  The reference peak is always the
    closest preceding peak in the array for which ampl_factor is 0.


"""
########################################################################

import numpy as num
import copy

#import mpfit
from   utils.mpfit import nmpfit as mpfit
from   detector import mca_calib as calib
import xrf_bgr

SIGMA_TO_FWHM = 2.35482

########################################################################
class XrfPeak:
    """
    Class defining a single xrf peak.
    The following inputs may be set by kw argument.
        self.label          = ""      # Peak label
        self.ignore         = False   # Don't fit peak
        self.energy         = 0.      # Peak energy
        self.energy_flag    = 1       # Flag for fitting energy
                                      #   0 = Optimize energy
                                      #   1 = Fix energy 
        self.fwhm           = 0.      # Peak FWHM
        self.fwhm_flag      = 1       # Flag for fitting FWHM
                                      #   0 = Optimize FWHM
                                      #   1 = Fix FWHM to global curve
                                      #   2 = Fix FWHM to input value
        self.ampl           = 0.      # Peak amplitude
        self.ampl_factor    = 0.      # Fixed amplitude ratio to previous peak
                                      #   0.  = Optimize amplitude of this peak
                                      #   >0. = Fix amplitude to this value relative
                                      #         to amplitude of previous free peak
                                      #  -1.0 = Fix amplitude at 0.0
        self.max_sigma      = 8.      # Max sigma for peak contribution
    """
    ###########################################################################
    def __repr__(self):
        lout = 'Xrf Peak:  Label = %s, Ignore = %s\n' % (self.label, str(self.ignore))
        lout = lout + '   energy    = %10.3f, flag = %i\n' % (self.energy, self.energy_flag)
        lout = lout + '   fwhm      = %10.3f, flag = %i\n' % (self.fwhm, self.fwhm_flag)
        lout = lout + '   amplitude = %10.3f, amp factor = %i\n' % (self.ampl, self.ampl_factor)
        lout = lout + '   area      = %10.3f\n' % (self.area)
        return lout

    ###########################################################################
    def __init__(self,**kws):
        self.label          = ""      # Peak label
        self.ignore         = False   # Don't fit peak
        self.energy         = 0.      # Peak energy
        self.energy_flag    = 1       # Flag for fitting energy
                                      #   0 = Optimize energy
                                      #   1 = Fix energy 
        self.fwhm           = 0.      # Peak FWHM
        self.fwhm_flag      = 1       # Flag for fitting FWHM
                                      #   0 = Optimize FWHM
                                      #   1 = Fix FWHM to global curve
                                      #   2 = Fix FWHM to input value
        self.ampl           = 0.      # Peak amplitude
        self.ampl_factor    = 0.      # Fixed amplitude ratio to previous peak
                                      #   0.  = Optimize amplitude of this peak
                                      #   >0. = Fix amplitude to this value relative
                                      #         to amplitude of previous free peak
                                      #  -1.0 = Fix amplitude at 0.0
        self.max_sigma      = 8.      # Max sigma for peak contribution
        self.area           = 0.      # Area of peak

        self.init(params=kws) 
    
    ###########################################################################
    def init(self,params=None):
        if params:
            keys = params.keys()
            if 'label' in keys:       self.label       = str(params['label'])
            if 'ignore' in keys:      self.ignore      = bool(params['ignore'])
            if 'energy' in keys:      self.energy      = float(params['energy'])
            if 'energy_flag' in keys: self.energy_flag = int(params['energy_flag'])
            if 'fwhm' in keys:        self.fwhm        = float(params['fwhm'])
            if 'fwhm_flag' in keys:   self.fwhm_flag   = int(params['fwhm_flag'])
            if 'ampl' in keys:        self.ampl        = float(params['ampl'])
            if 'ampl_factor' in keys: self.ampl_factor = float(params['ampl_factor'])
            if 'max_sigma' in keys:   self.max_sigma   = float(params['max_sigma'])
            if 'area' in keys:        self.area        = float(params['area'])

    ###########################################################################
    def get_params(self):
        """
        Return a dictionary of peak parameters
        """
        params = {'label':self.label,  
                  'ignore':self.ignore,
                  'energy':self.energy,
                  'energy_flag':self.energy_flag,
                  'fwhm':self.fwhm,
                  'fwhm_flag':self.fwhm_flag,
                  'ampl':self.ampl,
                  'ampl_factor':self.ampl_factor,
                  'max_sigma':self.max_sigma,
                  'area':self.area}
        return params

    ###########################################################################
    def calc(self, energy, compute_area=True):
        """
        Calculate the gaussian line shape.
        Inputs:  energy
        Output:  counts
            Returns a num array containing the predicted counts.
        """

        sigma = self.fwhm/SIGMA_TO_FWHM
        counts = self.ampl * num.exp(-((energy - self.energy)**2 / (2. * sigma**2)))

        if compute_area:
            #self.area = counts.sum()
            self.area = num.trapz(counts,energy)

        return( counts )
    
    ###########################################################################
    def _calc_range(self, energy, compute_area=True):
        """
        Calculate the gaussian line shape only within energy range that it
        makes a significant contribution.
        Inputs: energy
        Output: (counts, (idx_min,idx_max))
            Returns a num array containing the predicted counts.
            The array is limited in length to only the range that the peak makes
            a significant contribution.  The (idx_min, idx_max) tuple are the indicies
            of the energy array corresponding to the peak limits
        """
        (idx_min, idx_max) = self._en_range(energy)
        sigma = self.fwhm/SIGMA_TO_FWHM
        counts = self.ampl * num.exp(-((energy[idx_min:idx_max] - self.energy)**2 / (2. * sigma**2)))

        if compute_area:
            #self.area = counts.sum()
            self.area = num.trapz(counts,energy[idx_min:idx_max])

        return( counts, (idx_min,idx_max) )

    ######################################################################
    def _en_range(self, energy):
        """
        Calculate the range over which the peak makes a significant contribution.
        """

        sigma = self.max_sigma*(self.fwhm/SIGMA_TO_FWHM)
        nchan = len(energy)

        # find the index closest to the peak energy
        del_e    = num.abs(energy - self.energy)
        idx_cen  = num.where( del_e == min(del_e) )
        idx_cen  = idx_cen[0][0]
        
        # assume energy in ascending order and linear
        # slope is energy/channels
        slope = (energy[nchan-1] - energy[0]) / nchan
        del_sig_chan = abs(sigma / slope)
        
        idx_min = idx_cen -  int(del_sig_chan/2.)
        if idx_min > nchan -1:
            idx_min = nchan-1
        elif idx_min < 0:
            idx_min = 0
        
        idx_max = idx_cen +  int(del_sig_chan/2.)
        if idx_max > nchan -1:
            idx_max = nchan-1
        elif idx_max < 0:
            idx_max = 0

        return( (idx_min,idx_max) )

#######################################################################################
class XrfSpectrum:
    """
    Class defining an xrf spectrum.
    Fields:
        self.data                  = []        # Data to be fit
        self.channels              = []        # Corresponding channel numbers
        self.peaks                 = []        # List of XrfPeak instances
        self.bgr                   = None      # Background model, see xrf_bgr module
                                               #   if none, we'll assume here that the data
                                               #   has been background subtracted
        ### Parameters 
        self.energy_offset         =  0.       # Energy calibration offset (keV)
        self.energy_slope          =  1.       # Energy calibration slope 
        self.energy_flag           =  0        # Energy flag
                                               #   0 = Optimize energy calibration coefficients
                                               #   1 = Fix energy calibration coefficients
        self.fwhm_offset           =  0.15     # FWHM model offset (keV)
        self.fwhm_slope            =  0.       # FWHM model slope
        self.fwhm_flag             =  0        # Fwhm flag
                                               #   0 = Optimize FWHM coefficients
                                               #   1 = Fix FWHM coefficients
        self.chi_exp               =  0.       # Exponent of data for weights
                                               #   w = 1/y**chi_exp
                                               #   0.0 = ones
                                               #   0.5 = sqrt (poisson stats)
        self.max_eval              =  0        # Maximum number of function evaluations
                                               # Default sets this to 0 which
                                               # does not limit the number of function 
                                               # evaluations
        self.max_iter              =  20       # Maximum number of iterations
        self.tolerance             =  1.e-4    # Convergence tolerance. The fitting
                                               # process will stop when the value of 
                                               # chi**2 changes by a relative amount
                                               # less than tolerance on two successive 
                                               # iterations. 
    #################################################################################
    init and __init__ functions can take parameters dictionary with the following
    (same as returned by get_params)
    
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

    ###########################################################################
    def __repr__(self):
        lout = 'Xrf Spectrum Fit Statistics:\n'
        lout = lout + '   nchan   = %i\n' % self.nchan 
        lout = lout + '   npeaks  = %i\n' % self.npeaks 
        lout = lout + '   nparams = %i\n' % self.nparams
        lout = lout + '   neval   = %i\n' % self.n_eval
        lout = lout + '   n_iter  = %i\n' % self.n_iter
        lout = lout + '   chisqr  = %f\n' % self.chisqr
        lout = lout + '   status  = %i\n' % self.status
        lout = lout + '   err_str = %s\n' % self.err_string
        lout = lout + 'Xrf Spectrum Fit Parameters:\n'
        lout = lout + '   energy_offset = %f\n' % self.energy_offset
        lout = lout + '   energy_slope  = %f\n' % self.energy_slope
        lout = lout + '   energy_flag   = %i\n' % self.energy_flag
        lout = lout + '   fwhm_offset   = %f\n' % self.fwhm_offset
        lout = lout + '   fwhm_slope    = %f\n' % self.fwhm_slope
        lout = lout + '   fwhm_flag     = %i\n' % self.fwhm_flag
        lout = lout + '   chi_exp       = %f\n' % self.chi_exp
        lout = lout + '   max_eval      = %i\n' % self.max_eval
        lout = lout + '   max_iter      = %i\n' % self.max_iter
        lout = lout + '   tolerance     = %f\n' % self.tolerance
        if self.bgr:
            lout = lout + self.bgr.__repr__()
        for peak in self.peaks:
            lout = lout + peak.__repr__()
        return lout
    
    ################################################################################
    def __init__(self,data=None,chans=None,params={},guess=True):
        # data
        self.data                  = []        # Data to be fit
        self.channels              = []        # Corresponding channel numbers

        #Peaks
        self.peaks                 = []        # list of peaks instances

        # Background
        self.bgr                   = None      # background model

        # Parameters
        self.energy_offset         =  0.       # Energy calibration offset
        self.energy_slope          =  1.       # Energy calibration slope
        self.energy_flag           =  0        # Energy flag
                                               #   0 = Optimize energy calibration coefficients
                                               #   1 = Fix energy calibration coefficients
        self.fwhm_offset           =  0.15     # FWHM offset
        self.fwhm_slope            =  0.       # FWHM slope
        self.fwhm_flag             =  0        # Fwhm flag
                                               #   0 = Optimize FWHM coefficients
                                               #   1 = Fix FWHM coefficients
        self.chi_exp               =  0.       # Exponent of data for weights
                                               #   w = 1/y**chi_exp
                                               #   0.0 = ones
                                               #   0.5 = sqrt (poisson stats)
        self.max_eval              =  0        # Maximum number of function evaluations
        self.max_iter              =  20       # Maximum number of iterations
        self.tolerance             =  1.e-4    # Convergence tolerance
        
        # Determined by fit
        self.predicted             = []        # Theory curve
        self.weights               = []        # Fit weights
        self.parinfo               = []        # Fit parameter info (transient)
        self.nchan                 =  0        # number of channels to fit
        self.npeaks                =  0        # number of peaks to fit
        self.nparams               =  0        # number of fit parameters
        self.n_eval                =  0        # Actual number of function evalutions
        self.n_iter                =  0        # Actual number of iterations
        self.chisqr                =  0.       # Chi-squared on output
        self.status                =  0        # Output status code
        self.err_string            =  ''       # Output error string

        # init
        self.init(data=data,chans=chans,params=params,guess=guess)

    ##########################################################################################
    def init(self,data=None,chans=None,params={},guess=False,bgr_model='Default',calc=False):
        """
        Set/Reset data and parameters.

        Should be safe to pass data/chans = None and not change if already set

        Note params here are assumed to be the same as output by get_params 

        Note default guess = False ==>  do not
        guess at peak parameters.  Therefore, pass guess=True
        to have routine guess at peak params
        """
        ###### set data
        if data != None:
            self.data     = num.asarray(data,dtype=float)
        if chans != None:
            self.channels = num.asarray(chans,dtype=num.int)

        # Check data/channels numbers
        # if data is set, but not channels, assume that channels
        # corresponds to array indicies of data
        if (len(self.channels) == 0) and len(self.data > 0):
            self.channels = num.asarray(range(len(self.data)),dtype=num.int)

        # make sure arrays match
        self.nchan = len(self.channels)
        if len(self.data) != self.nchan:
            if len(self.data) > 0:
                print len(self.data), len(self.channels)
                raise "Data/Channels array mismatch error in peak fit"
        if self.nchan < 1:
            print "Warning no channels"

        ###### get parameters (see self.get_params)
        fit_par  = params.get('fit')
        bgr_par  = params.get('bgr')
        peak_par = params.get('pk')

        ###### set/reset fit parameters
        if fit_par != None:
            keys = fit_par.keys()
            if 'energy_offset' in keys: self.energy_offset = float(fit_par['energy_offset'])
            if 'energy_slope' in keys: self.energy_slope   = float(fit_par['energy_slope'])
            if 'energy_flag'  in keys: self.energy_flag    = int(fit_par['energy_flag'])
            if 'fwhm_offset' in keys: self.fwhm_offset     = float(fit_par['fwhm_offset'])
            if 'fwhm_slope' in keys: self.fwhm_slope       = float(fit_par['fwhm_slope'])
            if 'fwhm_flag' in keys: self.fwhm_flag         = int(fit_par['fwhm_flag'])
            if 'chi_exp' in keys: self.chi_exp             = float(fit_par['chi_exp'])
            if 'max_eval' in keys: self.max_eval           = int(fit_par['max_eval'])
            if 'max_iter' in keys: self.max_iter           = int(fit_par['max_iter'])
            if 'tolerance' in keys: self.tolerance         = float(fit_par['tolerance'])
            
        ####### set peaks if passed
        # note if passed this assumes we are
        # setting all of them!
        if peak_par != None:
            self.peaks = []
            # each entry should be the result
            # of peak.get_param
            for par in peak_par:
                tmp = XrfPeak()
                tmp.init(params=par)
                self.peaks.append(tmp)

        ####### set bgr if passed
        if bgr_par != None:
            if bgr_model == 'Default':
                self.bgr = xrf_bgr.Background()
                self.bgr.init(params=bgr_par)
                
        ####### Now some more init of internal stuff
        self._initFitParams()
        self._initPeaks(guess=guess)

        ####### Calc the model if flag is set
        if calc: self.calc()

    ################################################################################
    def _initFitParams(self,):
        """ Reset the internal data/parameters """
        self.predicted = num.zeros(self.nchan)
        self.weights   = num.ones(self.nchan)
        self.npeaks    = len(self.peaks)
        self.parinfo   = []
        self.nparams   = 0
        self.n_eval    = 0
        self.n_iter    = 0
        self.chisqr    = 0.
        self.n_iter    = 0
        self.chisqr    = 0.
        self.status    = 0
        self.err_string = ''

    ################################################################################
    def _initPeaks(self,guess=True):
        """ Init peaks"""
        #for peak in self.peaks:
        for j in range(len(self.peaks)):
            self._initSinglePeak(idx=j,guess=guess)

    def _initSinglePeak(self,idx=0,guess=True):
        """ Init a single peak"""
        if self.data == None: guess = False
        
        peak = self.peaks[idx]
        
        # Don't fit peaks outside the energy range of the data
        chan = calib.energy_to_channel(peak.energy,
                                       offset=self.energy_offset,
                                       slope=self.energy_slope) 
        if ((chan < self.channels[0]) or (chan > self.channels[self.nchan-1])):
            peak.ignore = True

        # Peak fwhm.
        # Note fwhm flag:
        #   0 = Optimize FWHM
        #   1 = Fix FWHM to global curve
        #   2 = Fix FWHM to input value
        if (peak.fwhm_flag == 1):
            peak.fwhm = (self.fwhm_offset + self.fwhm_slope*num.sqrt(peak.energy))
        elif (peak.fwhm_flag == 0) and (guess==True):
            # use the same approx for guess fwhm?
            peak.fwhm = (self.fwhm_offset + self.fwhm_slope*num.sqrt(peak.energy))

        # Peak ampl
        # Note peak.ampl_factor
        #   0.  = Optimize amplitude of this peak
        #   >0. = Fix amplitude to this value relative
        #         to amplitude of previous free peak
        #  -1.0 = Fix amplitude at 0.0
        if (peak.ignore == True) or (peak.ampl_factor < 0.):
            peak.ampl = 0.
        elif (peak.ampl_factor == 0.) and (guess == True):
            chan = int( (peak.energy - self.energy_offset) / self.energy_slope )
            chan = min( max(chan, 0), (self.nchan-1) )
            peak.ampl = max(self.data[chan], 0.)
            last_opt_peak = peak
        elif (peak.ampl_factor > 0.):
            peak.ampl = last_opt_peak.ampl * peak.ampl_factor
            # Don't correct for FWHM here, this is just initial value???
            # peak.ampl = peak.ampl * (last_opt_peak.fwhm / max(peak.fwhm, .001))

    ################################################################################
    def get_params(self,):
        """
        Return fit parameters in a dictionary
        Should we also get peak and bgr params???
        """
        # fit parameters 
        fit_par = {'energy_offset':self.energy_offset,
                   'energy_slope':self.energy_slope,
                   'energy_flag':self.energy_flag,
                   'fwhm_offset':self.fwhm_offset,
                   'fwhm_slope':self.fwhm_slope,
                   'fwhm_flag':self.fwhm_flag,
                   'chi_exp': self.chi_exp,
                   'max_eval': self.max_eval,
                   'max_iter': self.max_iter,
                   'tolerance': self.tolerance}

        # bgr params
        if self.bgr:
            bgr_par = self.bgr.get_params()
        else:
            bgr_par = {}

        # peak params
        peak_par = []
        for pk in self.peaks:
            peak_par.append(pk.get_params())

        return ({'fit':fit_par,'bgr':bgr_par,'pk':peak_par}) 

    #########################################################################
    def get_data(self):
        """
        Returns the data array.
        """
        return self.data

    #########################################################################
    def get_energy(self):
        """
        Returns the energy array
        """
        #calc energy with updated params
        en = calib.channel_to_energy(self.channels, offset=self.energy_offset,
                                     slope=self.energy_slope)
        return en

    ####################################################################################
    def calc(self, compute_areas=True):
        """
        Predicts a Gaussian spectrum 
        """
        self.predicted = num.zeros(self.nchan, dtype=num.float)
        energy = self.get_energy()
        
        for peak in self.peaks:
            (counts, (idx_min, idx_max) )   = peak._calc_range(energy, compute_area=compute_areas)
            self.predicted[idx_min:idx_max] = self.predicted[idx_min:idx_max] + counts

        if self.bgr:
            self.bgr.calc(self.data,slope=self.energy_slope)
            self.predicted = self.predicted + self.bgr.bgr

    ####################################################################################
    def calc_peaks(self,):
        """
        Return array with predicted values for each peak 
        """
        npks   = len(self.peaks)
        cnts   = num.zeros((npks,self.nchan), dtype=num.float)
        energy = self.get_energy()

        for j in range(npks):
            (tmp, (idx_min, idx_max) ) = self.peaks[j]._calc_range(energy)
            cnts[j,idx_min:idx_max] = tmp

        # could use this also for calulating the total predicted signal
        # self.predicted = cnt.sum(axis=1) + self.bgr.calc(self.data,slope=self.energy_slope)
        
        return cnts

    ####################################################################################
    def fit(self,guess=True,opt_bgr=True, quiet=1):
        """
        Fit the data
        """
        if self.bgr:
            self.bgr.calc(self.data,slope=self.energy_slope)
            if opt_bgr == False:
                data = copy.copy(self.data)
                self.data = data - self.bgr.bgr
                bgr_model = self.bgr
                self.bgr  = None
        else:
            bgr_model = None
            opt_bgr   = False

        # Prep and call lsq
        self._preFit(guess=guess)
        functkw = {'fit':self}
        m = mpfit.mpfit(_fit_peaks, parinfo=self.parinfo, functkw=functkw, 
                        quiet=quiet, xtol=self.tolerance, maxiter=self.max_iter)

        # Make sure final results are updated
        self._update(m.params)
        self.calc(compute_areas=True)
        if opt_bgr == False and bgr_model != None:
            self.bgr  = bgr_model
            self.data = data
            self.predicted = self.predicted + self.bgr.bgr

        # some of the results
        if (m.status <= 0): print m.errmsg
        self.n_iter = m.niter
        self.n_eval = m.nfev
        self.chisqr = m.fnorm
        self.status = m.status
        self.err_string = m.errmsg

    ############################################################################################
    def _preFit(self,guess=True):
        """
        Constructs fit weights and param info
        """
        self._initFitParams()
        self._initPeaks(guess=guess)
        
        # Compute sigma of observations to computed weighted residuals
        # Default is ones set in _init_fit_
        if  (self.chi_exp > 0.0):
            self.weights = num.asarray(self.data,dtype=float)

            # get rid of zeros
            idx = num.where(self.weights < 1.)
            if len(idx) > 0:  self.weights[idx] = num.ones(len(idx))

            # Treat special cases of self.chi_exp=.5, 1.
            if (self.chi_exp == 0.5):
                self.weights = 1./num.sqrt(self.weights)
            elif (self.chi_exp == 1.0):
                self.weights = 1./self.weights
            else:
                self.weights = 1./((self.weights)**self.chi_exp)

        # Total number of fit parameters
        self.nparams = self.npeaks*3 + 4

        # Create the Parameter info structure for initial guesses and constraints
        # and other parameter info to be used by mpfit
        self.parinfo = []
        for i in range(self.nparams):
            self.parinfo.append({'value':0., 'fixed':0, 'limited':[0,0],
                                 'limits':[0., 0.], 'step':0.})

        # Energy calibration offset
        np = 0
        self.parinfo[np]['value']   = self.energy_offset
        self.parinfo[np]['parname'] = 'Energy offset'
        if (self.energy_flag == 1):
            self.parinfo[np]['fixed']=1

        # Energy calibration slope
        np = np+1
        self.parinfo[np]['value']   = self.energy_slope
        self.parinfo[np]['parname'] = 'Energy slope'
        if (self.energy_flag == 1):
            self.parinfo[np]['fixed']=1

        # Global FWHM offset
        np = np+1
        self.parinfo[np]['value']   = self.fwhm_offset
        self.parinfo[np]['parname'] = 'FWHM offset'
        if (self.fwhm_flag == 1):
            self.parinfo[np]['fixed']=1

        # Global FWHM slope
        np = np+1
        self.parinfo[np]['value']   = self.fwhm_slope
        self.parinfo[np]['parname'] = 'FWHM slope'
        if (self.fwhm_flag == 1):
            self.parinfo[np]['fixed']=1

        # Peaks
        for peak in self.peaks:
            # Peak energy
            # Note energy flag:
            # 0 = Optimize energy
            # 1 = Fix energy 
            np = np+1
            self.parinfo[np]['value']   = peak.energy
            self.parinfo[np]['parname'] = peak.label + ' energy'
            if (peak.energy_flag == 1):
                self.parinfo[np]['fixed']=1

            # Peak fwhm
            # Note fwhm flag:
            #   0 = Optimize FWHM
            #   1 = Fix FWHM to global curve
            #   2 = Fix FWHM to input value
            np = np+1
            self.parinfo[np]['value']   = peak.fwhm
            self.parinfo[np]['parname'] = peak.label + ' FWHM'
            if (peak.fwhm_flag != 0):
                self.parinfo[np]['fixed']=1
            else: 
                # Limit the FWHM to .1 to 10 times initial guess
                self.parinfo[np]['limited'] =[1,1]
                self.parinfo[np]['limits']  =[peak.fwhm/10., peak.fwhm*10.]

            # Peak ampl
            # Note peak.ampl_factor
            #   0.  = Optimize amplitude of this peak
            #   >0. = Fix amplitude to this value relative
            #         to amplitude of previous free peak
            #  -1.0 = Fix amplitude at 0.0
            np = np+1
            self.parinfo[np]['value'] = peak.ampl
            self.parinfo[np]['parname'] = peak.label + ' amplitude'
            if (peak.ignore == True) or (peak.ampl_factor < 0.):
                self.parinfo[np]['fixed'] = 1
            elif (peak.ampl_factor == 0.):
                # Limit the amplitude to non-negative values
                self.parinfo[np]['limited'] =[1,0]
                self.parinfo[np]['limits']  =[0.,0.]
            elif (peak.ampl_factor > 0.):
                self.parinfo[np]['fixed']=1

        # get bgr fit parameter info
        if self.bgr:
            nbgr = len(self.bgr.parinfo)
            self.nparams = self.nparams + nbgr
            for j in range(nbgr):
                np = np+1
                self.parinfo.append(self.bgr.parinfo[j])

        # number of fit parameters and degrees of freedom
        # self.nparams = np
        if len(self.parinfo) != self.nparams :
            print "Why are these not the same?"
            print len(self.parinfo), self.nparams

    ################################################################################
    def _update(self, parameters):
        """
        Update model parameters from lsq parameter vector
        Note order of params same as above
        """
        # Energy calibration offset
        np = 0
        self.energy_offset = parameters[np]

        # Energy calibration slope
        np = np + 1
        self.energy_slope = parameters[np]

        # Global FWHM offset
        np = np + 1
        self.fwhm_offset = parameters[np]

        # Global FWHM slope
        np = np + 1
        self.fwhm_slope = parameters[np]

        # Peaks
        for peak in self.peaks:
            # Peak energy
            np = np + 1
            peak.energy = parameters[np]

            # Peak fwhm
            # Note fwhm flag:
            #   0 = Optimize FWHM
            #   1 = Fix FWHM to global curve
            #   2 = Fix FWHM to input value
            np = np + 1
            if (peak.fwhm_flag == 1):
                peak.fwhm = (self.fwhm_offset + 
                             self.fwhm_slope*num.sqrt(peak.energy))
            else:
                peak.fwhm = parameters[np]

            # Peak amplitude
            # Note peak.ampl_factor
            #   0.  = Optimize amplitude of this peak
            #   >0. = Fix amplitude to this value relative
            #         to amplitude of previous free peak
            #  -1.0 = Fix amplitude at 0.0
            np = np + 1
            if (peak.ignore == True) or (peak.ampl_factor < 0.):
                peak.ampl = 0.
            elif (peak.ampl_factor == 0.):
                peak.ampl = parameters[np]
                last_opt_peak = peak
            elif (peak.ampl_factor > 0.):
                peak.ampl = (last_opt_peak.ampl * peak.ampl_factor)
                peak.ampl = peak.ampl * (last_opt_peak.fwhm / max(peak.fwhm, .001))

        # pass the remaining params
        # to the bgr model
        if self.bgr:
            np = np+1
            self.bgr._update(parameters[np:])

#########################################################################
def _fit_peaks(parameters, fjac = None, fit = None):
    """ Private function """
    fit._update(parameters)
    fit.calc(compute_areas=False)
    status = 0
    res = (fit.predicted - fit.data) * fit.weights
    return (status, res )

########################################################################
########################################################################
########################################################################
def test_peak():
    import pylab
    # make some dat
    import _test_dat as test_dat
    chans = num.arange(2048)
    offset = 1.0
    slope = .01
    en = offset + slope*chans
    data = test_dat.data1(en)
    pylab.plot(en,data,'ko')
    #
    p1 = XrfPeak(label='1',energy=4.,ampl=550,fwhm=.5)
    print p1.get_params()
    p2 = XrfPeak(label='1',energy=7.,ampl=850,fwhm=.6)
    print p2.get_params()
    pylab.plot(en,p1.calc(en),'r')
    pylab.plot(en,p2.calc(en),'r')
    #
    (y,(mi,ma)) = p1._calc_range(en)
    print mi,ma
    pylab.plot(en[mi:ma],y,'g')
    #
    pylab.show()

########################################################################
def test_fit():
    import pylab
    import _test_dat as test_dat
    from   xrf_bgr import Background
    #############################
    # Data
    #############################
    chans = num.arange(2048)
    offset = 1.0
    slope = .01
    en = offset + slope*chans
    data = test_dat.data1(en)
    pylab.subplot(211)
    pylab.plot(en,data,'ko')

    #############################
    # fit
    #############################
    p1    = XrfPeak(label='p1',energy=4.)
    p2    = XrfPeak(label='p2',energy=7.)
    bgr   = Background(bottom_width=4,compress=4)
    xspec = XrfSpectrum(data=data,chans=chans,peaks=[p1,p2],bgr=bgr,
                        energy_offset=1.,energy_slope=0.01,chi_exp=0.5,
                        guess=True)

    # lets see how good the initial guess is
    pylab.plot(en,xspec.peaks[0].calc(en),'r')
    pylab.plot(en,xspec.peaks[1].calc(en),'r')
    
    # Now do the fit
    xspec.fit(opt_bgr=False,quiet=0)
    print xspec
    pylab.subplot(212)
    en = xspec.energy_offset + xspec.energy_slope*chans
    pylab.plot(en,data,'ko')
    pylab.plot(en,xspec.predicted,'g')
    pylab.plot(en,xspec.bgr.bgr,'k-')

    #############################
    pylab.show()
    
########################################################################
if __name__ == "__main__":
    #test_peak()
    test_fit()
