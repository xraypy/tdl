
from Num import Num, num_version
import types
import time
import Med
import Mca
import fitPeaks
import fitBgr
import xrf_lookup

#####################################################
class XRF:
    def __init__(self,med=None,detectors=[],total=True,align=True,tau=None):

        self.med         = med
        self.array_len   = 0
        self.detectors   = map(int,detectors)
        self.total       = total
        self.align       = align

        # lists for corrected/processed data
        self.data        = []
        self.calibration = []
        self.bgr         = []
        self.predicted   = []
        
        # lists for fitting parameters
        self.bgr_params  = []
        self.fit_params  = []
        self.peak_params = []

        # rois for corrected data
        # self.rois        = []

        # update the data arrays
        self.init_data(tau=tau)

        return

    #########################################################################
    # needs work -- should provide summary of fitting parameters etc..
    def __repr__(self):
        if self.med == None: return "No Data"
        lout = 'XRF Data Object:\n'
        lout = lout + '  Name = %s\n' % self.med.name
        lout = lout + '  Num detectors = %i\n' % self.med.n_detectors
        lout = lout + '  Detectors used = %s\n' % str(self.detectors)
        lout = lout + '  Align = %s\n' % str(self.align)
        lout = lout + '  Total = %s\n' % str(self.total)
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

    #def __copy__(self):
    #    new = XRF()
    #    return new

    #def __deepcopy__(self,visit):
    #    new = XRF()
    #    return new

    #########################################################################
    def init_data(self,detectors=None,total=None,align=None,correct=True,tau=None):

        if self.med == None: return

        # update stuff
        if total is not None:
            self.total = total
        if align is not None:
            self.align = align

        # note tweak this so if any of the detectors are negative
        # that means skip
        if detectors is not None:
            self.detectors = []
            if len(detectors) > 0 and len(detectors) < self.med.n_detectors:
                self.detectors = map(int,detectors)
        
        # make sure we have a detector list
        if self.detectors == []:
            self.detectors = range(self.med.n_detectors) 
        #for j in range(len(self.bad)):
        #    if bad[j] in detectors: detectors.remove(bad[j])

        # update deadtime correction factors
        self.med.update_correction(tau=tau)
 
        # clear arrays
        self.data = []
        self.calibration = []

        #get data 
        self.data = self.med.get_data(detectors=self.detectors,
                                      total=self.total,
                                      align=self.align,
                                      correct=correct)
        #self.array_len = len(self.detectors)
        self.array_len = len(self.data)

        # get calibration from med
        if self.total:
            #self.array_len = 1
            idx = int(self.detectors[0])
            self.calibration.append(self.med.mcas[idx].calibration)
        elif self.align:
            idx = int(self.detectors[0])
            for d in self.detectors:
                idx = int(d)
                self.calibration.append(self.med.mcas[idx].calibration)
        else:
            for d in self.detectors:
                idx = int(d)
                self.calibration.append(self.med.mcas[idx].calibration)

        # recalc rois if needed
        # ....

        # reinit the bgr and peak parameters
        # note this blows away everything!  
        self.init_params()

        return


    #########################################################################
    def detector_idx(self,d):
        """given a detector number, find which index it corresponds
        to in the detector list"""
        idx = Num.where(self.detectors==d)
        try:
            idx = idx[0]
        except:
            idx = -1
        return idx
    
    #########################################################################
    def init_params(self,lines = [], bottom_width=4.,top_width=0.,
                    exponent=2, tangent=0, compress=4):

        # reset all peak and bgr parameters to defaults
        
        self.fit_params  = self.array_len * [[]]
        self.bgr_params  = self.array_len * [[]]
        self.peak_params = self.array_len * [[]]

        self.bgr         = self.array_len * [[]]
        self.predicted   = self.array_len * [[]]

        for i in range(self.array_len):
            self.fit_params[i] = fitPeaks.McaFit()
            self.bgr_params[i] = fitBgr.McaBackground(slope=self.calibration[i].slope,
                                                      exponent=exponent,
                                                      top_width=top_width,
                                                      bottom_width=bottom_width,
                                                      tangent=tangent,
                                                      compress=compress)
            for line in lines:
                #self.set_peak(idx=i,line=line)
                self.init_peak_line(idx=i,line=line)
                
        return

    #########################################################################
    def get_params(self):
        return (self.fit_params, self.bgr_params, self.peak_params)

    #######################################################################
    def init_peak_line(self,idx=0,line=None,lookup=True):
        """
        create a new peak_param given an xrf line
        """
        
        if not idx in range(self.array_len): return

        # if line = None, blow away all peaks for this detector
        #if line == None:
        #    self.peak_params[idx] = []
        #    return

        if lookup==True and type(line)==types.StringType:
            label = line
            en = xrf_lookup.lookup_xrf_line(line)
        else:
            #label = 'peak:' + str(line)
            #label = str(line)
            #en = line
            return -1

        # see if its a duplicate 
        dup = self.peak_idx(idx=idx,label=label)
        if dup > -1:
            self.peak_params[idx].pop(dup)

        # create a new McaPeak structure
        tmp = fitPeaks.McaPeak()
        tmp.label = label
        tmp.initial_energy = en
        tmp.energy         = en
        self.peak_params[idx].append(tmp)

        # return pk_idx of new peak_param
        return (len(self.peak_params[idx]) - 1)

    #######################################################################
    def init_peak_en(self,idx=0,label=None,energy=0.0):
        """
        create a new peak_param given a name and energy
        """
        
        if not idx in range(self.array_len): return

        # if line = None, blow away all peaks for this detector
        #if line == None:
        #    self.peak_params[idx] = []
        #    return

        if label == None:
            label = str(energy)

        # see if its a duplicate 
        dup = self.peak_idx(idx=idx,label=label)
        if dup > -1:
            self.peak_params[idx].pop(dup)

        # create a new McaPeak structure
        tmp = fitPeaks.McaPeak()
        tmp.label = label
        tmp.initial_energy = energy
        tmp.energy         = energy
        self.peak_params[idx].append(tmp)

        # return pk_idx of new peak_param
        return (len(self.peak_params[idx]) - 1)


    #########################################################################
    def peak_idx(self,idx=0,label=None):
        """ See if a peak label exists and ret idx"""
        if not idx in range(self.array_len):
            return -1
        if label == None:
            return -1

        for j in range(len(self.peak_params[idx])):
            if self.peak_params[idx][j].label == label:
                return j
        return -1        

    #########################################################################
    def set_bgr(self,idx=0,slope=None,exponent=None,top_width=None,
                bottom_width=None, tangent=None, compress=None):
        
        # note idx refers to corrected/processed arrays not to the detector number

        if not idx in range(self.array_len): return
        if slope        is not None: self.bgr_params[idx].slope        = slope
        if exponent     is not None: self.bgr_params[idx].exponent     = exponent
        if top_width    is not None: self.bgr_params[idx].top_width    = top_width
        if bottom_width is not None: self.bgr_params[idx].bottom_width = bottom_width
        if tangent      is not None: self.bgr_params[idx].tangent      = tangent
        if compress     is not None: self.bgr_params[idx].compress     = compress

        return

    #########################################################################
    def set_peak(self,idx=0,pk_idx=0,energy=None,ampl=None,fwhm=None,
                 energy_flag=None,fwhm_flag=None,ampl_factor=None,ignore=False):

        # note idx refers to corrected/processed arrays not to the detector number
        # add peak fitting parameters...

        if not idx in range(self.array_len): return
        if not pk_idx in range(len(self.peak_params[idx])): return

        if energy is not None:
            self.peak_params[idx][pk_idx].initial_energy = energy
            self.peak_params[idx][pk_idx].energy = energy

        if ampl is not None:
            self.peak_params[idx][pk_idx].initial_ampl = ampl
            self.peak_params[idx][pk_idx].ampl = ampl
            
        if fwhm is not None:
            self.peak_params[idx][pk_idx].initial_fwhm = fwhm
            self.peak_params[idx][pk_idx].fwhm = fwhm

        if energy_flag is not None:
            self.peak_params[idx][pk_idx].energy_flag = energy_flag
            
        if fwhm_flag is not None:
            self.peak_params[idx][pk_idx].fwhm_flag = fwhm_flag

        if ampl_factor is not None:
            self.peak_params[idx][pk_idx].ampl_factor = ampl_factor

        self.peak_params[idx][pk_idx].ignore = ignore

        return

    #########################################################################
    def get_data(self):
        return self.data

    #########################################################################
    def get_energy(self):
        # note bugsssss
        # shouldnt we use idx = int(self.detectors[0])
        # etc instead of  calib[0]
        en = self.array_len*[[]]
        if self.total:
            channels = Num.arange(len(self.data[0]))
            en[0] = self.calibration[0].channel_to_energy(channels)
        elif self.align:
            channels = Num.arange(len(self.data[0]))
            ref_en = self.calibration[0].channel_to_energy(channels)
            for i in range(self.array_len):
                en[i] = ref_en
        else:
            for i in range(self.array_len):
                channels = Num.arange(len(self.data[i]))
                en[i] = self.calibration[i].channel_to_energy(channels)
        return en

    #########################################################################
    def get_count_totals(self):
        lst  = []
        for mca in self.med.mcas:
            rt = mca.elapsed.real_time
            lt = mca.elapsed.live_time
            OCR = mca.elapsed.total_counts/mca.elapsed.live_time
            ICR = mca.elapsed.input_counts/mca.elapsed.live_time
            lst.append({'rt':lt,'lt':lt,'OCR':OCR,'ICR':ICR})
        return lst
            

    #########################################################################
    def fit(self,fwhm_flag=1,energy_flag=1,chi_exp=0.0):

        for i in range(self.array_len):
            self._fit_bgr(i)
            self._fit_peaks(i,fwhm_flag=fwhm_flag,energy_flag=energy_flag,chi_exp=chi_exp)
        return 

    #########################################################################
    def fit_bgr(self):
        """
        Remove the background from a spectrum
        """
        for i in range(self.array_len):
            self._fit_bgr(i)
        return

    def _fit_bgr(self,idx):

        # note idx refers to corrected/processed arrays not to the detector number

        if not idx in range(self.array_len): return
        
        self.bgr_params[idx].slope = self.calibration[idx].slope
        self.bgr[idx] = fitBgr.fit_background(self.data[idx], self.bgr_params[idx])

        return

    #########################################################################
    def fit_peaks(self,fwhm_flag=1,energy_flag=1,chi_exp=0.0):

        for i in range(self.array_len):
            self._fit_peaks(i,fwhm_flag=fwhm_flag,energy_flag=energy_flag,chi_exp=chi_exp)
        return 

    def _fit_peaks(self,idx,fwhm_flag=1,energy_flag=1,chi_exp=0.0):

        # note idx refers to corrected/processed arrays not to the detector number
        
        if not idx in range(self.array_len): return
        
        observed     = self.data[idx]
        background   = self.bgr[idx]
        # in case background was not subtracted
        if len(background) == 0:
            background = Num.zeros(len(observed))
        calibration  = self.calibration[idx]
        peaks        = self.peak_params[idx]
        fit          = self.fit_params[idx]

        # Copy parameters to fit
        fit.npeaks = len(peaks)
        fit.initial_energy_offset = calibration.offset
        fit.initial_energy_slope  = calibration.slope
        fit.nchans = len(observed)
        fit.last_chan = fit.nchans-1
        fit.fwhm_flag = fwhm_flag
        fit.energy_flag = energy_flag
        fit.chi_exp   = chi_exp
        
        t0 = time.time()
        [fit, peaks, fit_counts] = fitPeaks.fitPeaks(fit, peaks, observed - background)
        t1 = time.time()
        
        self.predicted[idx]   = fit_counts + background
        self.fit_params[idx]  = fit
        self.peak_params[idx] = peaks
        self.calibration[idx].offset = fit.energy_offset
        self.calibration[idx].slope  = fit.energy_slope

        return
    
    #########################################################################
    def calc_peaks(self):

        for i in range(self.array_len):
            self._calc_peaks(i)
        return 

    def _calc_peaks(self,idx):

        # note idx refers to corrected/processed arrays not to the detector number
        
        if not idx in range(self.array_len): return

        observed     = self.data[idx]
        background   = self.bgr[idx]
        # in case background was not subtracted
        if len(background) == 0:
            background = Num.zeros(len(observed))
        calibration  = self.calibration[idx]
        peaks        = self.peak_params[idx]
        fit          = self.fit_params[idx]

        # Copy parameters to fit
        fit.npeaks = len(peaks)
        fit.energy_offset = calibration.offset
        fit.energy_slope  = calibration.slope
        fit.nchans = len(observed)
        fit.last_chan = fit.nchans-1
        
        fit_counts = fitPeaks.predict_gaussian_spectrum(fit, peaks)
        self.predicted[idx] = fit_counts + background
        
        return
