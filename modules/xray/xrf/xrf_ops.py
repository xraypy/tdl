#######################################################################
"""
T. Trainor (fftpt@uaf.edu)
Xrf methods for data handling and fitting

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

import sys
import types
import string
import copy
import numpy as Num

from  detector import medfile_cars
from  detector import medfile_emsa
from  detector import mca_calib as calib
from  xrf_model import Xrf

##############################################################################
def read(file,bad_mca_idx=[],total=True,align=True,correct=True,tau=None,
         det_idx=0,emin=-1.0,emax=-1.0,fmt='CARS',xrf_params={},lines=None):
    """
    Read detector files
    >>m = xrf.read(file="file_name",bad_mca_idx=[],
                   total=True,align=True,tau=None)

    if xrf_params == None: returns a list of med objects
    otherwise: returns a list of xrf objects (default)

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
            tmp = rd(file=f, bad_mca_idx=bad_mca_idx,
                     total=total,align=align,correct=correct,tau=tau)
            med.append(tmp)
    else:
        return None

    if xrf_params == None:
        return med
    else:
        # generate xrf from med 
        xrf = med2xrf(med,xrf_params=xrf_params,lines=lines,
                      det_idx=det_idx,emin=emin,emax=emax)
        return xrf

##############################################################################
def read_files(prefix,start=0,end=100,nfmt=3,bad_mca_idx=[],
               total=True,align=True,correct=True,tau=None,det_idx=0,
               emin=-1.0,emax=-1.0,fmt='CARS',xrf_params={},lines=None):
    """
    Read multiple files
    if xrf_params == None: returns a list of med objects
    otherwise: returns a list of xrf objects (default)
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

    if xrf_params == None:
        return med
    else:
        # generate xrf from med 
        xrf = med2xrf(med,xrf_params=xrf_params,lines=lines,
                      det_idx=det_idx,emin=emin,emax=emax)
        return xrf

##############################################################################
def med2xrf(med,xrf_params={},lines=None,det_idx=0,emin=-1.,emax=-1.):
    """
    Given an med object (or list of med objects) return a
    new xrf object (or list of xrf objects).
    
    Note xrf_params should have the same format as returned
    by Xrf.get_params
    """
    if type(med) == types.ListType:
        xrf = []
        for m in med:
            tmp = _med2xrf(m,xrf_params=xrf_params,lines=lines,
                           det_idx=det_idx,emin=emin,emax=emax)
            xrf.append(tmp)
    else:
        xrf = _med2xrf(med,xrf_params=xrf_params,lines=lines,
                       det_idx=det_idx,emin=emin,emax=emax)
    return xrf

def _med2xrf(med,xrf_params={},lines=None,det_idx=0,emin=-1.,emax=-1.):
    # Xrf just fits one data/energy array    
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
    
    # Create an XrfSpectrum object
    x =  Xrf(data=data,chans=chans,params=xrf_params)

    #Add bgr    
    if xrf_params.has_key('bgr'):
        x.init_bgr(params=xrf_params['bgr'])
    else:
        x.init_bgr()

    # if lines add lines
    if lines:
        x.init_lines(lines)

    return x

#################################################################################
def fit(xrf,xrf_params={},use_prev_fit=False,fit_init=-1,guess=False,verbose=True):
    """
    Given a (list of) xrf objects fit them all

    fit_init: the index of the first scan in the list to fit
              it will be used as the seed value for fitting
              all the xrf objects restarting at index zero.
    use_prev_fit: if True then use index-1 as seed parameters

    If fit_init >-1 and use_prev_fit=True then fit_init will only
    be used as the seed for fitting the first index only.
    """
    if type(xrf) != types.ListType:
        xrf = [xrf]

    if fit_init >=0:
        if verbose: sys.__stdout__.write("Fitting index = %d\n" % fit_init)
        xrf[fit_init].fit()
        params = xrf[fit_init].get_params()
        init = True
    elif len(xrf_params) > 0:
        params = xrf_params
        init = True
    else:
        init = False
    
    for j in range(len(xrf)):
        if verbose: sys.__stdout__.write("Fitting index = %d\n" % j)
        if (j > 0) and (use_prev_fit == True):
            params = xrf[j-1].get_params()
            xrf[j].init(params=params)
        elif init:
            xrf[j].init(params=params,guess=guess)
        xrf[j].fit()
    return

#################################################################################
def fit_scan_files(prefix='xrf',start=0,end=100,nfmt=3,
                   bad_mca_idx=[],total=True,align=True,tau=None,
                   det_idx=0,emin=-1.0,emax=-1.0,fmt='CARS',xrf_params={},
                   use_prev_fit=False,fit_init=-1):
    """
    Fit a bunch of files.  Dont really need this...
    """
    xrf = read_files(prefix=prefix,start=start,end=end,nfmt=nfmt,
                     bad_mca_idx=bad_mca_idx,total=total,align=align,tau=tau,
                     det_idx=det_idx,emin=emin,emax=emax,fmt=fmt,
                     xrf_params=xrf_params)
    
    fit(xrf,xrf_params=xrf_params,use_prev_fit=use_prev_fit,fit_init=fit_init)
    return xrf

#################################################################################
def peak_areas(xrf,line):
    """
    get peak areas for given line
    """
    if type(xrf) != types.ListType:
        xrf = [xrf]
    results = []
    for x in xrf:
        for pk in x.peaks:
            if string.lower(pk.label) == string.lower(line):
                results.append(pk.area)
                break
    return Num.array(results)

#################################################################################
def xrf_plot(xrf,d='Data',f='Fit',p=None,ylog=True,xlog=False,hold=False):
    """
    Plot options:
    d = 'Data', 'Data-Bgr'
    f = 'Fit', 'Bgr', 'Fit-Bgr', 'Fit and Bgr'
    p = 'Peaks', 'Peaks+Bgr'
    """
    tiny = 1.e-6

    try:
        import pylab
    except:
        return
    if hold == False: pylab.clf()
    
    en  = xrf.get_energy() + tiny
    da  = xrf.get_data() + tiny
    fit = copy.copy(xrf.predicted)
    fit = Num.array(fit) + tiny
    bgr = copy.copy(xrf.bgr.bgr)
    bgr = Num.array(bgr)  + tiny  

    if d == 'Data':
        pylab.plot(en,da,'k.-',label='Data')
    if d == 'Data-Bgr':
        if len(bgr) == len(da):
            pylab.plot(en,da-bgr,'k.-',label='Data')
        else:
            pylab.plot(en,da,'k.-',label='Data')

    if f == 'Fit':
        pylab.plot(en,fit,'r-',label='Fit')
    elif f == 'Bgr':
        pylab.plot(en,bgr,'r-',label='Bgr')
    elif f == 'Fit-Bgr':
        if len(bgr) == len(da):
            pylab.plot(en,fit-bgr,'r-',label='Fit-Bgr')
        else:
            pylab.plot(en,fit,'r-',label='Fit')
    elif f == 'Fit and Bgr':
        pylab.plot(en,fit,'r-',label='Fit')
        pylab.plot(en,bgr,'g-',label='Bgr')

    if p == 'Peaks':
        npks = len(xrf.peaks)
        lbls = []
        for pk in xrf.peaks:
            lbls.append(pk.label)
        pk_fit = xrf.calc_peaks()
        for j in range(npks):
            pylab.plot(en,pk_fit[j],label=lbls[j])
    elif p == 'Peaks+Bgr':
        npks = len(xrf.peaks)
        lbls = []
        for pk in xrf.peaks:
            lbls.append(pk.label)
        pk_fit = xrf.calc_peaks()
        for j in range(npks):
            if len(bgr) == len(pk_fit[j]):
                pylab.plot(en,pk_fit[j]+bgr,label=lbls[j])
            else:
                pylab.plot(en,pk_fit[j],label=lbls[j])

    if ylog:
        pylab.semilogy()
        pylab.ylim(ymin=1.)
    if xlog:
        pylab.semilogx()

    # annotation
    pylab.legend()
    pylab.xlabel('keV')
    pylab.ylabel('counts')


##############################################################################
##############################################################################
if __name__ == "__main__":
    import pylab
    xrf = read(file='_test.xrf',bad_mca_idx=[0,2,13],emin=4.,emax=9.)
    pylab.plot(xrf.get_energy(),xrf.get_data(),'-k')
    #pylab.show()
    #
    xrf.init_lines(['Fe ka',7.12])
    #xrf.init_bgr()
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
    
