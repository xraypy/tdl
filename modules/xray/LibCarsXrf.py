# T. Trainor (2006)
# Wrapper functions for CARS MED/MCA library
# T. Trainor, 8-11-2006
#
# --------------
# Modifications
# --------------
#
#
# --------------
# Todo
# --------------
#
# * add read/fit spectra command (ie loop)
#
# * modify things so that detectors = [-1,-2] means use all
#   except the negative values.
#
# * currently rois are directly from med.  should change to use corrected data in xrf class
#   and work on show rois etc....
#
# * test fitting etc when total = False
#
# * work on __repr__, ie nicer format, include fitting results
#
# * add Emin/Emax to read/set_detectors
#
#
# Note for read and set commands, should have an mcas arg instead of detectors
# For others use a det=integer syntax.  (or index=integer) ie data, energy
# fits etc should work on individual detectors or all of them at once
#
#
# -> mca_inc_array = [1,2,3,-4,5,6,7,-9]
# -> can use bad = [1,2]
# -> mca_inc_array = [1,2,3,-4,5,6,7,-9]
#
# -> put plot cmd in here and add all fancy options....
#
#
#################################################################################################
HelpXRF = """

tdl xrf functions:
------------------
xrf.read()
xrf.set_detectors() 

xrf.data()
xrf.energy()

set_roi()
xrf.get_rois()

xrf.set_bgr()
xrf.set_peak()

xrf.fit_bgr(substract=False)
xrf.fit_peaks(fit=True)
xrf.calc_peaks()

xrf.get_peaks()
xrf.data_dict()

xrf.plot()

------------------

Note functions that take a detectors argument:  detectors refer to the MED detector number,
starting from zero.  If detectors have been skipped (eg in set_data) then returned arrays
are indexed from 0 - len(detectors) -1 (unless total=True, then just one value)

Example read/plot/fit:
>>xrf = xrf.read('test/test2.xrf')
>>det = [1,2,4,5,6,7,8,9,10,11,12,14]
>>xrf.set_data(xrf,det)
>>d = xrf.data(xrf)
>>e = xrf.energy(xrf)
>>plot(e[0],d[0])
>>xrf.init_peaks(xrf,['Ca ka','Ca kb','Ti ka', 'Fe ka', 'Fe kb'])
>>pred = xrf.fit_peaks(xrf)
>>plot(e[0],pred[0])

"""

##########################################################################
from Num import Num, num_version

import os
import sys
import types
import copy
from Util import datalen

import CarsMcaFile
import Mca
import Med
import fitPeaks
import fitBgr
import XRF

title = "CARS XRF Library "

##############################################################################
def xrf_read(file=None,detectors=[],total=True,align=True,tau=None,tdl=None,**kws):
    """
    Read detector files
    >>m = xrf.read("file_name",detectors=[],total=True,align=True,tau=None)

    Returns an xrf object.  This function is appropriate for reading single and
    multi-element detectors.  (We always assume that the detector may be a
    multi-element detector)

    Inputs:
        file:
            File name
            
        detectors:
            A list of detectors to be used.  An empty
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
    xrf = XRF.XRF(med=med,detectors=detectors,total=total,align=align,tau=tau)
    return xrf

def xrf_read_cmd(val,file=None,tdl=None,**kws):
    """use command syntax, put ret into xrf.file"""
    name = val.med.name
    if '.' in name:
        #name = name.split('.',1)[0]
        name = name.replace('.','_')
    name = 'xrf_data.%s' % name
    tdl.setVariable(name,val=val)
    return

#############################################################################
def xrf_set_data(xrf,detectors=[],total=None,align=None,tau=None):
    """
    Reset detectors used for xrf analysis
    >>xrf.set_data(xrf,detectors=[],total=False,align=True,tau=[])

    Inputs:
        xrf:
            An instance of the xrf class or an xrf file name. See xrf.read
        
        detectors:
            A list of detectors to be used.  An empty list (default) means
            use all the detectors. Note detector indexing starts at zero!

        total:
            Set this keyword to work on the sum of the spectra from all
            of the detectors as a 1-D Numeric array.
            
        align:
            Set this keyword to return spectra which have been shifted and
            and stretched to match the energy calibration parameters of the
            first detector.  This permits doing arithmetic on a
            "channel-by-channel" basis. This keyword can be used alone
            or together with the TOTAL keyword, in which case the data
            are aligned before summing.

        tau:
            None or a list of deadtime factors for each detector (dimesion should
            be the same as the total number of detectors)
                     ocr = icr * exp(-icr*tau)
                     cor = (icr/ocr)*(rt/lt)
                     max_icr = 1/tau
                     max_ocr = max_icr*exp(-1)

    Notes:
        Setting any of the above keywords toggles thier values for
        any additional analysis.  This routine also clears all previous
        arrays and paramterss relating to analysis (background and peaks)!  
    """
    if type(xrf) == types.StringType:
        med = CarsMcaFile.read_med(file=xrf)
        xrf = XRF.XRF(med=med,detectors=detectors,total=total,align=align,tau=tau)
    else:
        xrf.init_data(detectors=detectors,total=total,align=align,tau=tau)
    return

#############################################################################
def xrf_data(xrf,detectors=[]):
    """
    Returns the data 
    >>data = xrf.data(xrf,detectors=[])

    Inputs:
        xrf:
            An instance of the xrf class. See xrf.read
        
        detectors:
            A list of detectors to return.  An empty list (default) means
            return all the detectors.  Note detector indexing
            starts at zero!  This argument is ignored if xrf.total = True
        
    Outputs:
        By default this function returns [counts] where
         counts = [ [cnts0],[cnts1], ...]

        These counts array is length = len(detectors), the first entry corresponds to
        the first detector in the detector list etc..

        If the align flag is set to True the first detector is used as the reference
        for alignment (all energy arrays have the same values).
        
        If the "total" keyword is set the counts array has only one entry.
        The counts are the sum for all detectors

    Example:
        >>data = xrf.data(xrf)  # xrf.total = True
        >>total_counts = data[0]

        >>data = xrf.data(xrf)  # xrf.total = False
        >>cnts0 = data[0]
        >>cnts1 = data[1]
        
    """
    cnts = xrf.get_data()
    if xrf.total == True or detectors == []:
        return cnts
    else:
        ret_cnts = []
        for d in detectors:
            idx = xrf.detector_idx(d)
            if idx >= 0:
                ret_cnts.append(cnts[idx])
        return ret_cnts

#############################################################################
def xrf_energy(xrf,detectors=[]):
    """
    Returns the energy arrays
    >>data = xrf.energy(xrf,detectors=[])

    Inputs:
        xrf:
            An instance of the xrf class or an xrf file name.

        detectors:
            A list of detectors to return.  An empty list (default) means
            return all the detectors.  Note detector indexing
            starts at zero!  This argument is ignored if xrf.total = True
        
    Outputs:
        By default this function returns [energy,counts] where
         energy = [ [e0],[e1], ...]

        The energy array length = len(detectors), the first entry corresponds to
        the first detector in the detector list etc..

        If the align flag is set to True the first detector is used as the reference
        for alignment (all energy arrays have the same values).
        
        If the "total" keyword is set the energy and counts arrays only have one enetry each.
        The counts are the sum for all detectors and then the energy for the first detector

    Example:
        >>energy = xrf.data(xrf)
        >>en0 = data[0]
        >>en1 = data[1]
        
    """
    energy = xrf.get_energy()
    if xrf.total == True or detectors == []:
        return energy
    else:
        ret_energy = []
        for d in detectors:
            idx = xrf.detector_idx(d)
            if idx >= 0:
                ret_energy.append(energy[idx])
        return ret_energy

#########################################################################
def xrf_set_bgr(xrf,exponent=None,top_width=None, bottom_width=None,
                tangent=None,compress=None,detectors=[]):
    """
    Set background parameters for detector(s)
    >>xrf.set_bgr(xrf,exponent=None,top_width=None,bottom_width=None,
                  tangent=None, compress=None,detectors=[])

    * Inputs
        xrf:
            An instance of the xrf class
            
        bottom_width:
            Specifies the width of the polynomials which are concave downward.
            The bottom_width is the full width in energy units at which the
            magnitude of the polynomial is 100 counts. The default is 4.

        top_width:
            The width of the polynomials which are concave upward.
            The top_width is the full width in energy units at which the
            magnitude of the polynomial is 100 counts. The default is 0, which
            means that concave upward polynomials are not used.

        exponent:
            Specifies the power of polynomial which is used. The power must be
            an integer. The default is 2, i.e. parabolas. Higher exponents,
            for example EXPONENT=4, results in polynomials with flatter tops
            and steeper sides, which can better fit spectra with steeply
            sloping backgrounds.

        tangent: flag 0/1
           Specifies that the polynomials are to be tangent to the slope of the
           spectrum. The default is 0 = vertical polynomials. This option works
           best on steeply sloping spectra. It has trouble in spectra with
           big peaks because the polynomials are very tilted up inside the
           peaks.

        compress:
           Compression factor to apply before fitting the background.
           Default=4, which means, for example, that a 2048 channel spectrum
           will be rebinned to 512 channels before fitting.

        detectors:
            A list of detectors to set parameters for.  An empty list (default) means
            set for all the detectors.  Note detector indexing starts at zero!
            Note if xrf.total == True this argument is ignored since the background is
            substracted after summing the detectors (only one set of parameters needed)
    """
    if xrf.total == True:
        xrf.set_bgr(exponent=exponent,top_width=top_width,
                    bottom_width=bottom_width, tangent=tangent, compress=compress)
    else:
        if detectors == []: detectors = xrf.detectors
        for d in detectors:
            idx = xrf.detector_idx(d)
            if idx >= 0:
                xrf.set_bgr(idx=idx,exponent=exponent,top_width=top_width,
                            bottom_width=bottom_width, tangent=tangent, compress=compress)

    return

#########################################################################
def xrf_fit_bgr(xrf,subtract=False,detectors=[]):
    """
    Remove the background from a spectrum
    >>data = xrf.fit_bgr(xrf,subtract=True,detectors=[])

    * Inputs
        xrf:
            An instance of the xrf class

        subtract: True/False
            Set this keyword to return background subtracted data.  False
            returns the background functions.  Defaults is False.
                    
        detectors:
            A list of detectors to fit.  An empty list (default) means
            fit all the detectors.  Note detector indexing starts at zero!
            Note if xrf.total == True this argument is ignored since the fit is
            performed after summing the detectors.

    * Outputs
        Returns an array [bgr0,bgr1,bgr2] 
        If the "total" keyword is set then the return only has a single entry
        If subtract == True the background substracted data is returned
    """
    if xrf.total == True:
        xrf.fit_bgr()
        if subtract:
            bgr = xrf.data[0] - xrf.bgr[0]
        else:
            bgr = xrf.bgr
    else:
        if detectors == []: detectors = xrf.detectors
        bgr = []
        for d in detectors:
            idx = xrf.detector_idx(d)
            if idx >= 0:
                xrf._fit_bgr(idx)
                if subtract:
                    bgr.append(xrf.data[idx] - xrf.bgr[idx])
                else:
                    bgr.append(xrf.bgr[idx])
    return bgr


#########################################################################
def xrf_init_peaks(xrf,peaks,detectors=[]):
    """
    Init peak parameters for detector(s)
    >>xrf.init_peaks(xrf,peaks,detectors=[])

    * Inputs
        xrf:
            An instance of the xrf class

        peaks:
            A list of xrf peaks or peak energies, eg ['Fe Ka', 'Fe Kb']
            If peaks = [] or None, then this blows away the previous
            values

        detectors:
            A list of detectors to set parameters for.  An empty list (default) means
            set for all the detectors.  Note detector indexing starts at zero!
            Note if xrf.total == True this argument is ignored since the fit is
            performed after summing the detectors (only one set of peak parameters needed)
    """
    # note if peaks = None or [] should blow away all previous peaks

    if xrf.total == True:
        if peaks == [] or peaks == None:
            xrf.init_peak(line=None)
        else:
            for line in peaks:
                xrf.init_peak(line=line)
    else:
        if detectors == []: detectors = xrf.detectors
        for d in detectors:
            idx = xrf.detector_idx(d)
            if idx >= 0:
                if peaks == [] or peaks == None:
                    xrf.init_peak(idx=idx,line=None)
                else:
                    for line in peaks:
                        xrf.init_peak(idx=idx,line=line)
    return

#########################################################################
def xrf_set_peak(xrf,label,energy=None,ampl=None,fwhm=None,energy_flag=None,
                  fwhm_flag=None,ampl_factor=None,ignore=False,detectors=[]):
    """
    Set peak parameters for detector(s)
    >>xrf.set_peak(xrf,label,energy=None,ampl=None,fwhm=None,energy_flag=None,
                  fwhm_flag=None,ampl_factor=None,ignore=False,detectors=[])

    * Inputs
        xrf:
            An instance of the xrf class

        label:
            A string describing the peak, use similiar syntax to init_peak

        energy:
            Peak energy

        ampl:
            Peak amplitude.  If ampl_factor is 0.0 then the fitting routine 
            will automaticailly determine a value for the initial_ampl

        fwhm:
            Peak full width at half max.
            This can be zero (Default) if fwhm_flag is 0
            
        energy_flag: 
            Flag for fitting energy of this peak
             0 = Fix energy (Default)
             1 = Optimize energy

        fwhm_flag:
            Flag for fitting FWHM of this peak
             0 = Fix FWHM to global curve (Default)
             1 = Optimize FWHM
             2 = Fix FWHM to input value

        ampl_factor:
            Flag for fitting amplitude of this peak
             0.0  = Optimize amplitude of this peak (Default)
             >0.0 = Fix amplitude to this value relative to amplitude of 
                    previous unconstrained peak

        ignore:
            Flag to ignore peak in the fit (Default is False)

        detectors:
            A list of detectors to set parameters for.  An empty list (default) means
            set for all the detectors.  Note detector indexing starts at zero!
            Note if xrf.total == True this argument is ignored since the fit is
            performed after summing the detectors (only one set of peak parameters needed)
    """
    if xrf.total == True:
        pk_idx = xrf.peak_idx(idx=0,label=label)
        if pk_idx == -1:
            pk_idx =xrf.init_peak(idx=0,line=label)

        xrf.set_peak(idx=0,pk_idx=pk_idx,energy=energy,ampl=ampl,fwhm=fwhm,
                     energy_flag=energy_flag,fwhm_flag=fwhm_flag,
                     ampl_factor=ampl_factor,ignore=ignore)

    else:
        if detectors == []: detectors = xrf.detectors
        for d in detectors:
            idx = xrf.detector_idx(d)
            if idx >= 0:
                pk_idx = xrf.peak_idx(idx=idx,label=label)
                if pk_idx == -1:
                    pk_idx = xrf.init_peak(idx=idx,line=label)

                xrf.set_peak(idx=0,pk_idx=pk_idx,energy=energy,ampl=ampl,fwhm=fwhm,
                             energy_flag=energy_flag,fwhm_flag=fwhm_flag,
                             ampl_factor=ampl_factor,ignore=ignore)
    return

#########################################################################
def xrf_fit_peaks(xrf,fit_bgr=True,fwhm_flag=1,energy_flag=1,chi_exp=0.0,detectors=[]):
    """
    Fit peaks
    >>xrf.fit_peaks(fit_bgr=True,fwhm_flag=1,energy_flag=1,chi_exp=0.0,detectors=[])

    * Inputs
        xrf:
            An instance of the xrf class

        fit_bgr: (True/False)
            Flag to indicate if the background should be fit (and removed)
            before fitting peaks

        fwhm_flag:
            0 = Fix FWHM coefficients
            1 = Optimize FWHM coefficients (Default)

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

        detectors:
            A list of detectors to fit.  An empty list (default) means
            fit all the detectors.  Note detector indexing starts at zero!
            Note if xrf.total == True this argument is ignored since the fit is
            performed after summing the detectors.

    * Outputs
        Returns an array [predicted0,predicted1,predicted2] 
        If the "total" keyword is set then the return only has a single entry
        The predicted function includes the background

    """
    if xrf.total == True:
        if fit_bgr == True:
            xrf.fit(fwhm_flag=fwhm_flag,energy_flag=energy_flag,chi_exp=chi_exp)
        else:
            xrf.fit_peaks(fwhm_flag=fwhm_flag,energy_flag=energy_flag,chi_exp=chi_exp)
        predicted = xrf.predicted
    else:
        if detectors == []: detectors = xrf.detectors
        predicted = []
        for d in detectors:
            idx = xrf.detector_idx(d)
            if idx >= 0:
                if fit_bgr: xrf._fit_bgr(idx)
                xrf._fit_peaks(idx,fwhm_flag=fwhm_flag,energy_flag=energy_flag,chi_exp=chi_exp)
                predicted.append(xrf.predicted[idx])
    return predicted

#########################################################################
def xrf_calc_peaks(xrf,detectors=[]):
    """
    Calculate peaks
    >>xrf.calc_peaks(xrf,detectors=[])

    * Inputs
        xrf:
            An instance of the xrf class

        detectors:
            A list of detectors to fit.  An empty list (default) means
            fit all the detectors.  Note detector indexing starts at zero!
            Note if xrf.total == True this argument is ignored since the fit is
            performed after summing the detectors.

    * Outputs
        Returns an array [predicted0,predicted1,predicted2] 
        If the "total" keyword is set then the return only has a single entry
        The predicted function includes the background
        
    """
    if xrf.total == True:
        xrf.calc_peaks()
        predicted = xrf.predicted
    else:
        predicted = []
        if detectors == []: detectors = xrf.detectors
        for d in detectors:
            idx = xrf.detector_idx(d)
            if idx >= 0:
                xrf._calc_peaks(idx)
                predicted.append(xrf.predicted[idx])
    return predicted

##############################################################################
def xrf_data_dict(xrf,detectors=[]):
    """ Dictionary of xrf data
    >>xrf.data_dict(xrf,detectors=[])

    * Inputs
        xrf:
            An instance of the xrf class

        detectors:
            A list of detectors to fit.  An empty list (default) means
            fit all the detectors.  Note detector indexing starts at zero!
            Note if xrf.total == True this argument is ignored since the fit is
            performed after summing the detectors.

    * Outputs
        Returns a dictionary 
          data = {'energy':[],'counts':[],'background':[],'predicted':[],'peak_areas':[],'peaks':[]}

        If the "total" keyword is set then each array has a single entry

    """
    data = {'energy':[],'counts':[],'background':[],'predicted':[],'peak_areas':[],'peaks':[]}
    counts = xrf.get_data()
    energy = xrf.get_energy()
    if xrf.total == True:
        data['energy']     = energy
        data['counts']     = counts
        data['background'] = xrf.bgr
        data['predicted']  = xrf.predicted

        # peak areas
        peak_areas = {}
        for peak in xrf.peak_params[0]:
            peak_areas[peak.label] = peak.area
        data['peak_areas'].append(peak_areas)

        # get each peak function
        peaks = {}
        for peak_params in xrf.peak_params[0]:
            calc = fitPeaks.predict_gaussian_spectrum(xrf.fit_params[0], [peak_params])
            peaks[peak_params.label] = calc
        data['peaks'].append(peaks)
        
    else:
        if detectors == []: detectors = xrf.detectors
        for d in detectors:
            idx = xrf.detector_idx(d)
            if idx >= 0:
                data['energy'].append(energy[idx])
                data['counts'].append(counts[idx])
                data['background'].append(xrf.bgr[idx])
                data['predicted'].append(xrf.predicted[idx])

                # peak areas
                peak_areas = {}
                for peak in xrf.peak_params[idx]:
                    peak_areas[peak.label] = peak.area
                data['peak_areas'].append(peak_areas)

                # get each peak function
                peaks = {}
                for peak_params in xrf.peak_params[idx]:
                    calc = fitPeaks.predict_gaussian_spectrum(xrf.fit_params[idx], [peak_params])
                    peaks[peak_params.label] = calc
                data['peaks'].append(peaks)

    return data

##############################################################################
def xrf_get_peaks(xrf,detectors=[]):
    """ Dictionary of peak fitting results
    >>xrf.get_peaks(xrf,detectors=[])

    Inputs
      xrf:
          An instance of the xrf class

      detectors:
          A list of detectors to fit.  An empty list (default) means
          fit all the detectors.  Note detector indexing starts at zero!
          Note if xrf.total == True this argument is ignored since the fit is
          performed after summing the detectors.

    Outputs
      Returns a dictionary 
        data = [{'label':{energy:0,ampl:0,fwhm:0,area:0},'label':{energy:0,ampl:0,fwhm:0,area:0}...}, {}]

      If the "total" keyword is set then each array has a single entry

    """
    results = []
    if xrf.total == True:
        # peak results
        peak_results = {}
        for peak in xrf.peak_params[0]:
            peak_results[peak.label] = {'energy':0.0,'ampl':0.0,'fwhm':0.0,'area':0.0}
            peak_results[peak.label]['energy'] = peak.energy
            peak_results[peak.label]['ampl'] = peak.ampl
            peak_results[peak.label]['fwhm'] = peak.fwhm
            peak_results[peak.label]['area'] = peak.area
        results.append(peak_results)
        
    else:
        if detectors == []: detectors = xrf.detectors
        for d in detectors:
            idx = xrf.detector_idx(d)
            if idx >= 0:
                peak_results = {}
                for peak in xrf.peak_params[idx]:
                    peak_results[peak.label] = {'energy':0.0,'ampl':0.0,'fwhm':0.0,'area':0.0}
                    peak_results[peak.label]['energy'] = peak.energy
                    peak_results[peak.label]['ampl'] = peak.ampl
                    peak_results[peak.label]['fwhm'] = peak.fwhm
                    peak_results[peak.label]['area'] = peak.area
                results.append(peak_results)

    return results

###############################################################################
def xrf_get_rois(xrf,detectors=[],total=False,net=False):
    """
    Return the roi values
    r = xrf.get_rois(xrf,detectors=[],total=False,net=False)

    Inputs:
        xrf:
            an instance of the xrf class or an xrf file name. See med.read

    Keywords:
        detectors:
            A list of detectors to be passed back.  An empty list (default)
            return all the detectors.  Note detector indexing
            starts at zero!

        total:
            If the total flag is True then the roi's from each detector
            will be summed

        net:
            If net is true then return the net (bgr subtracted) counts rather than the total counts
            
    Outputs:    
        Returns an array of roi values,
        r[0] = {roi1:value,roi2:value} --> for first det (or sum if total is true)...
    """

    if type(xrf) == types.StringType:
        med = CarsMcaFile.read_med(file=xrf)
        if detectors == []:
            detectors = range(med.n_detectors)
        else:
            detectors = map(int,detectors)        
    else:
        med = xrf.med
        if detectors == []:
            detectors = xrf.detectors
        else:
            detectors = map(int,detectors)        

    # note could just use net,total = roi[d].get_counts()
    # but that doesnt ret the label
    rois = med.get_rois()
    if not rois: return None

    ret = []
    for i in detectors:
        temp = {}
        for roi in rois[i]:
            if net==True:
                temp[roi.label]=roi.net
            else:
                temp[roi.label]=roi.total
        ret.append(temp)

    if total:
        temp = {}
        for i in range(len(ret)):
            if i == 0:
                temp = ret[i]
            else:
                for label in ret[i].keys():
                    temp[label] = temp[label] + ret[i][label]
        ret = [temp]

    return ret

def xrf_get_rois_cmd(val,**kws):
    for j in range(len(val)):
        print "detector %i:" % j
        print val[j]
    return

def xrf_show_rois(xrf):
    """Show detector counts/rois
    >>xrf.totals(xrf)

    Inputs:
        xrf:
            An instance of the xrf class. See xrf.read

    Output:
        A list of dictionaries.  Each entry in the list corresponds to a detector mca.
        The dictionary contains entries for real time (rt), live time (lt), output count rate (OCR)
        and input count rate (icr).  

    """
    print xrf.show_rois()

##########################################################################
def xrf_set_roi(med,label,lrn=[],detectors=[],units='keV'):
    """Set the roi values for a detector
    >>xrf.set_roi(xrf,label,lrn=[],detectors=[],units='keV')

    Inputs:
        xrf:
            An instance of the xrf class or an xrf file name. See xrf.read

        label:
            String label for the roi or an integer index (starting at zero)

    Keywords:
        lrn:
            lrn=[left,right,nbgr]
            if lrn=[] then delete it if it exists
            nbgr is optional.
             
        detectors:
            A list of detectors to set roi for.  An empty
            list (default) means use all the detectors.  Note detector indexing
            starts at zero!

        units:
            string specifying units for left and right. default is 'keV'
            Valid strings are "channel","keV" and "eV"
    """
    if type(xrf) == types.StringType:
        med = CarsMcaFile.read_med(file=xrf)
        if detectors == []:
            detectors = range(med.n_detectors)
        else:
            detectors = map(int,detectors)
    else:
        med = xrf.med
        if detectors == []:
            detectors = xrf.detectors
        else:
            detectors = map(int,detectors)

    # delete rois
    for i in detectors:
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
    else:
        # add roi
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
        for i in detectors:
            d = int(i)
            med.mcas[d].add_roi(roi)

    return


def xrf_get_count_totals(xrf):
    """Get detector count totals
    >>xrf.show_rois(xrf)

    Inputs:
        xrf:
            An instance of the xrf class. See xrf.read
    """
    return xrf.get_count_totals()

def xrf_get_count_totals_cmd(val,**kws):
    lout = ''
    print 'Number of MCAs = ', len(val)
    for j in range(len(val)):
        lout = lout + "MCA %3d, rt= %f,  lt= %f,  OCR= %6.0f,  ICR= %6.0f\n" % \
                      (j,val[j]['rt'],val[j]['lt'],val[j]['OCR'],val[j]['ICR'])
    print lout
    return 

##################################################################################

_help_  = {'xrf': HelpXRF}

# tdl functions
_func_ = {"xrf.read":(xrf_read,xrf_read_cmd),
          "xrf.set_data":xrf_set_data,
          "xrf.data":xrf_data,
          "xrf.energy":xrf_energy,
          "xrf.set_bgr":xrf_set_bgr,
          "xrf.fit_bgr":xrf_fit_bgr,
          "xrf.init_peaks":xrf_init_peaks,
          "xrf.set_peak":xrf_set_peak,
          "xrf.fit_peaks":xrf_fit_peaks,
          "xrf.calc_peaks":xrf_calc_peaks,
          "xrf.get_peaks":xrf_get_peaks,
          "xrf.data_dict":xrf_data_dict,
          "xrf.rois":(xrf_get_rois,xrf_get_rois_cmd),
          "xrf.show_rois":xrf_show_rois,
          "xrf.set_roi":xrf_set_roi,
          "xrf.totals":(xrf_get_count_totals,xrf_get_count_totals_cmd)}

_scripts_ = ['xrf_scripts.tdl']

