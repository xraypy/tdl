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
# * test fitting etc when total = False
#
# * work on __repr__, ie nicer format, include fitting results
#
# * add Emin/Emax to read/set_detectors
#
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
>>bad = [3,13]
>>xrf.set_data(xrf,bad_mca_idx=bad)
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
def xrf_read(file=None,bad_mca_idx=[],total=True,align=True,tau=None,tdl=None,**kws):
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

    xrf = XRF.read_xrf_file(file=file,bad_mca_idx=bad_mca_idx,total=total,align=align,tau=tau)
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
def xrf_set_data(xrf,bad_mca_idx=None,total=None,align=None,correct=None,tau=None):
    """
    Reset detectors used for xrf analysis
    >>xrf.set_data(xrf,detectors=[],total=False,align=True,correct=True,tau=[])

    Inputs:
        xrf:
            An instance of the xrf class or an xrf file name. See xrf.read
        
        bad_mca_idx:
            A list of bad detectors to be used.  An empty
            list (default) means use all the detectors.  Note detector indexing
            starts at zero!

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

        correct:
            Flag to determine if deadtime corrections are apprlied to the data
            
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
        xrf = XRF.read_xrf_file(file=file,bad_mca_idx=bad_mca_idx,total=total,align=align,tau=tau)
    else:
        xrf.init_data(bad_mca_idx=bad_mca_idx,total=total,align=align,correct=correct,tau=tau)
    return

#############################################################################
def xrf_data(xrf):
    """
    Returns the data 
    >>data = xrf.data(xrf)

    Inputs:
        xrf:
            An instance of the xrf class. See xrf.read
                
    Outputs:
        This function returns [counts] where
         counts = [ [cnts0],[cnts1], ...]

        These counts array is length = len(detectors), the first entry corresponds to
        the first detector in the detector list etc..

        If the align flag is set to True the first good detector is used as the reference
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
    return xrf.get_data()

#############################################################################
def xrf_energy(xrf):
    """
    Returns the energy arrays
    >>data = xrf.energy(xrf,detectors=[])

    Inputs:
        xrf:
            An instance of the xrf class or an xrf file name.
            
    Outputs:
        This function returns [energy,counts] where
         energy = [ [e0],[e1], ...]

        The energy array length = len(detectors), the first entry corresponds to
        the first detector in the detector list etc..

        If the align flag is set to True the first good detector is used as the reference
        for alignment (all energy arrays have the same values).
        
        If the "total" keyword is set the energy and counts arrays only have one enetry each.
        The counts are the sum for all detectors and then the energy for the first detector

    Example:
        >>energy = xrf.data(xrf)
        >>en0 = data[0]
        >>en1 = data[1]
        
    """
    return xrf.get_energy()


#########################################################################
def xrf_set_bgr(xrf,exponent=None,top_width=None, bottom_width=None,
                tangent=None,compress=None,det_idx=0):
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

        det_idx:
            Index of the detector to set bgr parameters for.  Note detector indexing starts at zero!
            Note if xrf.total == True then det_idx should be zero (default) 
    """
    xrf.set_bgr(exponent=exponent,top_width=top_width,
                bottom_width=bottom_width, tangent=tangent, compress=compress,det_idx=det_idx)
    return

#########################################################################
def xrf_fit_bgr(xrf,subtract=False):
    """
    Remove the background from a spectrum
    >>data = xrf.fit_bgr(xrf,subtract=True)

    * Inputs
        xrf:
            An instance of the xrf class

        subtract: True/False
            Set this keyword to return background subtracted data.  False
            returns the background functions.  Defaults is False.

    * Outputs
        Returns an array [bgr0,bgr1,bgr2] 
        If the "total" keyword is set then the return only has a single entry
        If subtract == True the background substracted data is returned
    """
    xrf.fit_bgr()
    if subtract:
        bgr = xrf.data - xrf.bgr
    else:
        bgr = xrf.bgr
    return bgr


#########################################################################
def xrf_init_peaks(xrf,peaks):
    """
    Init peak parameters for detector(s)
    >>xrf.init_peaks(xrf,peaks)

    * Inputs
        xrf:
            An instance of the xrf class

        peaks:
            A list of xrf peaks or peak energies, eg ['Fe Ka', 'Fe Kb']
            If peaks = None or [], then this blows away the previous
            values

    """
    # note if peaks = None blow away all previous peaks
    if peaks == [] or peaks == None:
        xrf.init_peak_line_all(None)
    else:
        for line in peaks:
            xrf.init_peak_line_all(line)

    return

#########################################################################
def xrf_set_peak(xrf,label,energy=None,ampl=None,fwhm=None,energy_flag=None,
                 fwhm_flag=None,ampl_factor=None,ignore=False,det_idx=0):
    """
    Set peak parameters for detector(s)
    >>xrf.set_peak(xrf,label,energy=None,ampl=None,fwhm=None,energy_flag=None,
                  fwhm_flag=None,ampl_factor=None,ignore=False,det_idx=0)

    * Inputs
        xrf:
            An instance of the xrf class

        label:
            A string describing the peak (may also be an integer corresponding to peak index)

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

        det_idx:
            Index for detector.  If det_idx == 'All' then the parameters are set for all detectors
            Note if xrf.total == then this should be zero (default)
    """

    if type(det_idx) == types.StringType:
        if det_idx.lower() == 'all':
            xrf.set_peak_all(pk_idx=pk_idx,energy=energy,ampl=ampl,fwhm=fwhm,
                             energy_flag=energy_flag,fwhm_flag=fwhm_flag,
                             ampl_factor=ampl_factor,ignore=ignore)
            return
        else:
            det_idx = int(det_idx)

    xrf.set_peak(pk_idx=pk_idx,energy=energy,ampl=ampl,fwhm=fwhm,
                 energy_flag=energy_flag,fwhm_flag=fwhm_flag,
                 ampl_factor=ampl_factor,ignore=ignore,det_idx=det_idx)

    return

#########################################################################
def xrf_fit_peaks(xrf,fwhm_flag=1,energy_flag=1,chi_exp=0.0,fit_bgr=True,):
    """
    Fit peaks
    >>xrf.fit_peaks(fit_bgr=True,fwhm_flag=1,energy_flag=1,chi_exp=0.0)

    * Inputs
        xrf:
            An instance of the xrf class

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

        fit_bgr: (True/False)
            Flag to indicate if the background should be fit (and removed)
            before fitting peaks (default = True)

    * Outputs
        Returns an array [predicted0,predicted1,predicted2] 
        If the "total" keyword is set then the return only has a single entry
        The predicted function includes the background

    """
    xrf.fit(fwhm_flag=fwhm_flag,energy_flag=energy_flag,chi_exp=chi_exp,fit_bgr=fit_bgr)
    predicted = xrf.predicted

    return predicted

#########################################################################
def xrf_calc_peaks(xrf):
    """
    Calculate peaks
    >>xrf.calc_peaks(xrf)

    * Inputs
        xrf:
            An instance of the xrf class

    * Outputs
        Returns an array [predicted0,predicted1,predicted2] 
        If the "total" keyword is set then the return only has a single entry
        The predicted function includes the background
        
    """
    xrf.calc_peaks()
    predicted = xrf.predicted
    return predicted

##############################################################################
def xrf_data_dict(xrf):
    """ Dictionary of xrf data
    >>xrf.data_dict(xrf)

    * Inputs
        xrf:
            An instance of the xrf class

    * Outputs
        Returns a dictionary 
          data = {'energy':[],'counts':[],'background':[],'predicted':[],'peak_areas':[],'peaks':[]}

        If the "total" keyword is set then each array has a single entry

    """
    dict = xrf.data_dict()
    return dict

##############################################################################
def xrf_get_peaks(xrf):
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
    results = xrf.get_peaks()
    return results

###############################################################################
def xrf_get_rois(xrf,background_width=1):
    """
    Return the roi values
    r = xrf.get_rois(xrf)

    Inputs:
        xrf:
            an instance of the xrf class or an xrf file name. See med.read
            
    Outputs:    
        Returns a list of dictionaries.  The list is of length num detectors
        each entry in the list holds a dictionary of {'lbl':(total, net),...}

    """
    return xrf.get_rois(background_width=background_width)

def xrf_get_rois_cmd(val,**kws):
    for j in range(len(val)):
        print "detector %i:" % j
        print val[j]
    return

###############################################################################
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
def xrf_set_roi(label,lrn=[],mcas=[],units='keV'):
    """Set the roi values for a detector
    >>xrf.set_roi(xrf,label,lrn=[],mcas=[],units='keV')

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
             
        mcas:
            A list of detectors to set roi for.  An empty
            list (default) means use all the detectors.  Note detector indexing
            starts at zero!

        units:
            string specifying units for left and right. default is 'keV'
            Valid strings are "channel","keV" and "eV"
    """

    xrf.xrf_set_roi(label,lrn=lrn,mcas=mcas,units=units)
    return

##########################################################################
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
_groups_ = [('xrf',True)]

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

