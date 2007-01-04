
from Num import Num, num_version
import types
import time
import copy
import Med
import Mca
import fitPeaks
import fitBgr
import xrf_lookup
import CarsMcaFile


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
def read_xrf_file(file=None,bad_mca_idx=[],total=True,align=True,tau=None):
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

    if file:
        med = CarsMcaFile.read_med(file=file)
    else:
        med = Med.Med()
    xrf = XRF(med=med,bad_mca_idx=bad_mca_idx,total=total,align=align,tau=tau)
    return xrf


#####################################################
class XRF:
    def __init__(self,med=None,bad_mca_idx=[],total=True,align=True,tau=None):

        self.med         = med
        self.ndet        = 0
        self.total       = total
        self.align       = align
        self.bad_mca_idx = bad_mca_idx

        # lists for corrected/processed data
        self.sum         = []
        self.data        = []
        self.bgr         = []
        self.predicted   = []
        
        # lists for fitting parameters
        self.bgr_params  = []
        self.fit_params  = []
        self.peak_params = []

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
        lout = lout + '  Bad Detectors = %s\n' % str(self.bad_mca_idx)
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
    def init_data(self,bad_mca_idx=None,total=None,align=None,correct=True,tau=None,init_params=True):

        if self.med == None: return

        # 
        if total is not None:
            self.total = total
        if align is not None:
            self.align = align

        # update bad list
        if bad_mca_idx is not None:
            self.bad_mca_idx = bad_mca_idx

        # update deadtime correction factors
        self.med.update_correction(tau=tau)

        # recalc rois
        self.med.update_rois(correct=correct)
 
        # clear arrays
        self.data = []

        #get data 
        self.data = self.med.get_data(bad_mca_idx=self.bad_mca_idx,
                                      total=self.total,
                                      align=self.align,
                                      correct=correct)

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

        # reset all peak and bgr parameters to defaults
        self.fit_params  = self.ndet * [[]]
        self.bgr_params  = self.ndet * [[]]
        self.peak_params = self.ndet * [[]]
        self.bgr         = self.ndet * [[]]
        self.predicted   = self.ndet * [[]]

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
                self.init_peak_line(line,det_idx=i)

        return

    #########################################################################
    def get_params(self):
        return (self.fit_params, self.bgr_params, self.peak_params)

    #########################################################################
    def get_data(self):
        return self.data

    #########################################################################
    def get_energy(self):
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
        if self.total or self.align:
            idx = self.med.get_align_idx(self.bad_mca_idx)
            return self.med.mcas[idx].calibration
        else:
            return self.med.mcas[det_idx].calibration

    #########################################################################
    def get_calibration_idx(self,det_idx):
        "get an mca's calibration"
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
    def init_peak_line(self,line, det_idx=0):
        """
        create a new peak_param given an xrf line
        """
        if not det_idx in range(self.ndet): return

        # if line == None, blow away all peaks for this detector
        if line == None:
            self.peak_params[det_idx] = []
            return

        if type(line)==types.StringType:
            line = line.strip()
            en = xrf_lookup.lookup_xrf_line(line)
        else:
            return 

        # see if its a duplicate 
        dup = self.peak_idx(det_idx=det_idx,label=line)
        if dup > -1:
            self.peak_params[det_idx].pop(dup)

        # create a new McaPeak structure
        tmp = fitPeaks.McaPeak()
        tmp.label          = line
        tmp.initial_energy = en
        tmp.energy         = en
        self.peak_params[det_idx].append(tmp)

        # return pk_idx of new peak_param
        return (len(self.peak_params[det_idx]) - 1)

    def init_peak_line_all(self,line):
        for det_idx in range(self.ndet):
            self.init_peak_line(line, det_idx=det_idx)
        return

    #######################################################################
    def init_peak_en(self,label=None, energy=0.0, det_idx=0):
        """
        create a new peak_param given a name and energy
        """
        
        if not det_idx in range(self.ndet): return

        if label == None:
            label = str(energy)
        label = label.strip()

        # see if its a duplicate 
        dup = self.peak_idx(det_idx=det_idx,label=label)
        if dup > -1:
            self.peak_params[det_idx].pop(dup)

        # create a new McaPeak structure
        tmp = fitPeaks.McaPeak()
        tmp.label          = label
        tmp.initial_energy = energy
        tmp.energy         = energy
        self.peak_params[idx].append(tmp)

        # return pk_idx of new peak_param
        return (len(self.peak_params[det_idx]) - 1)

    def init_peak_en_all(self,label=None, energy=0.0):

        for det_idx in range(self.ndet):
            self.init_peak_en(label=label, energy=energy, det_idx=det_idx)
        return

    #########################################################################
    def peak_idx(self,label=None,det_idx=0):
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
    def set_peak(self,pk_idx=0,energy=None,ampl=None,fwhm=None, energy_flag=None,
                 fwhm_flag=None,ampl_factor=None,ignore=False,det_idx=0):

        if not det_idx in range(self.ndet): return

        if type(pk_idx) == types.StringType:
            pk_idx = self.peak_idx(label=pk_idx.strip(),det_idx=det_idx)
        
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

    def set_peak_all(self,pk_idx=0,energy=None,ampl=None,fwhm=None, energy_flag=None,
                     fwhm_flag=None,ampl_factor=None,ignore=False):

        for det_idx in range(self.ndet):
            set_peak(pk_idx=pk_idx,energy=energy,ampl=ampl,fwhm=fwhm, energy_flag=energy_flag,
                     fwhm_flag=fwhm_flag,ampl_factor=ampl_factor,ignore=ignore,det_idx=det_idx)
        return

    #########################################################################
    def set_bgr(self,slope=None,exponent=None,top_width=None, bottom_width=None,
                tangent=None, compress=None,det_idx=0):

        if not det_idx in range(self.ndet): return
        if slope        is not None: self.bgr_params[det_idx].slope        = slope
        if exponent     is not None: self.bgr_params[det_idx].exponent     = exponent
        if top_width    is not None: self.bgr_params[det_idx].top_width    = top_width
        if bottom_width is not None: self.bgr_params[det_idx].bottom_width = bottom_width
        if tangent      is not None: self.bgr_params[det_idx].tangent      = tangent
        if compress     is not None: self.bgr_params[det_idx].compress     = compress

        return

    def set_bgr_all(self,slope=None,exponent=None,top_width=None, bottom_width=None,
                    tangent=None, compress=None):

        for det_idx in range(self.ndet):
            self.set_bgr(slope=slope,exponent=exponent,top_width=top_width, bottom_width=bottom_width,
                    tangent=tangent, compress=compress, det_idx=det_idx)
        return



    #########################################################################
    def fit(self,fwhm_flag=1,energy_flag=1,chi_exp=0.0,fit_bgr=True):

        for i in range(self.ndet):
            if fit_bgr: self._fit_bgr(i)
            self._fit_peaks(i,fwhm_flag=fwhm_flag,energy_flag=energy_flag,chi_exp=chi_exp)
        return 

    #########################################################################
    def fit_bgr(self):
        """
        Remove the background from a spectrum
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
            self.bgr[det_idx] = fitBgr.fit_background(self.data[det_idx], self.bgr_params[det_idx])

        return

    #########################################################################
    def fit_peaks(self,fwhm_flag=1,energy_flag=1,chi_exp=0.0):

        for i in range(self.ndet):
            self._fit_peaks(i,fwhm_flag=fwhm_flag,energy_flag=energy_flag,chi_exp=chi_exp)
        return 

    def _fit_peaks(self,det_idx,fwhm_flag=1,energy_flag=1,chi_exp=0.0):

        if not det_idx in range(self.ndet): return
        if det_idx in self.bad_mca_idx:
            self.predicted[det_idx] = Num.zeros(len(self.data[det_idx]))
            return
        
        observed     = self.data[det_idx]
        background   = self.bgr[det_idx]
        # in case background was not subtracted
        if len(background) == 0:
            background = Num.zeros(len(observed))
        calib_idx    = self.get_calibration_idx(det_idx)
        peaks        = self.peak_params[det_idx]
        fit          = self.fit_params[det_idx]

        # Copy parameters to fit
        fit.npeaks = len(peaks)
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
    def calc_peaks(self):

        for i in range(self.ndet):
            self._calc_peaks(i)
        return 

    #########################################################################
    def _calc_peaks(self,det_idx):

        if not det_idx in range(self.ndet): return

        if det_idx in self.bad_mca_idx:
            self.predicted[det_idx] = Num.zeros(len(self.data[det_idx]))
            return

        background   = self.bgr[det_idx]

        # in case background was not subtracted
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
        if not det_idx in range(self.ndet): return

        if type(pk_idx) == types.StringType:
            pk_idx = self.peak_idx(pk_idx)
            
        if not pk_idx in range(len(self.peak_params[det_idx])): return

        if det_idx in self.bad_mca_idx:
            self.predicted[det_idx] = Num.zeros(len(self.data[det_idx]))
            return

        background   = self.bgr[det_idx]

        # in case background was not subtracted
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
        dict = {'energy':[],'counts':[],'background':[],'predicted':[],'peak_areas':[],'peaks':[]}
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
                calc = fitPeaks.predict_gaussian_spectrum(self.fit_params[mca_idx], [peak_params])
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
            return self.med.get_roi_counts_lbl(background_width=background_width,correct=correct)

    #########################################################################
    def xrf_set_roi(self,label,lrn=[],mcas=[],units='keV'):

        if mcas == []: mcas = range(self.med.n_detectors)

        # delete rois
        for i in mcas:
            d = int(i)
            if type(label) == types.IntType:
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

