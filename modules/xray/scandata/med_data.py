#######################################################################
"""
T. Trainor (fftpt@uaf.edu)
Functions for operating on / working with
med's in scandata

Modifications:
--------------

"""
#######################################################################

import types
import numpy as num
from matplotlib import pyplot
import deadtime

########################################################################
def fit_deadtime(data,x='io',y='Med',norm='Seconds',offset=True,display=True):
    """
    Do a deadtime fit to data

      x = linear axis.  ie this should be proportional
          to the real input count.  default = 'io'
      y = 'Med'
      norm='Seconds'
      offset=True
      display=True

    returns the deadtime tau values.  

    Note make sure normalization is checked carefully
    e.g. ocr from the med is cps, to compare vs io should
    normalize io/count_time
    
    """
    if norm != None:
        norm = data[norm].astype(float)
        norm = 1./norm

    xfit = data[x]
    if norm != None:
        xfit = xfit * norm
    
    if y == 'Med':
        # this is cps, therefore dont apply norm
        yfit_arr = data.med_ocr()
    else:
        yfit_arr = data[y]
        if norm != None:
            yfit_arr = yfit_arr * norm
        yfit_arr = [yfit_arr]

    # do the fits
    tau = []
    a   = []
    if offset: off = []
    else: off = None
    for yfit in yfit_arr:
        params = deadtime.fit(xfit,yfit,offset=offset)
        tau.append(params[0])
        a.append(params[1])
        if offset:
            off.append(params[2])
        if display:  print params

    # update med taus...
    if y == 'Med':
        data.med_update_tau(tau)
    else:
        # for scaler calc correction
        # and post as y_c
        pass

    # show fits if wanted
    if display:
        (ndet,npts) = yfit_arr.shape
        ycorr_arr   = num.zeros((ndet,npts))
        if y == 'Med':
            for j in range(npts):
                for k in range(ndet):
                    cts = data.med[j].mca[k].get_data(correct=True)
                    ycorr_arr[k][j] = cts.sum()
            if norm != None:
                for k in range(ndet):
                    ycorr_arr[k] = ycorr_arr[k] * norm
        else:
            # need to apply correction to scaler above,
            # then get it here to plot
            pass
        
        _display_deadtime_fit(xfit,yfit_arr,ycorr_arr,tau,a,off)
        
    return tau

def _display_deadtime_fit(xfit,yfit_arr,ycorr_arr,tau,a,off):
    """
    plot x,y
    plot fits to x,y
    compute a corrected y and plot
    (for med compute by summing corrected data )
    """
    pyplot.clf()
    pyplot.subplot(2,1,1)

    # data
    for yfit in yfit_arr:
        pyplot.plot(xfit,yfit,'k.')
        
    # plot fit
    for j in range(len(tau)):
        if off != None:
            params = (tau[j],a[j],off[j])
            offset = True
        else:
            params = (tau[j],a[j])
            offset = False
        ycal = deadtime.calc_ocr(params,xfit,offset)
        pyplot.plot(xfit,ycal,'r-')
        
    # plot corrected data
    pyplot.subplot(2,1,2)
    for j in range(len(yfit_arr)):
        pyplot.plot(xfit,yfit_arr[j],'k.')
        pyplot.plot(xfit,ycorr_arr[j],'r-')
        if offset:
            pyplot.plot(xfit,a[j]*xfit+off[j],'k--')
        else:
            pyplot.plot(xfit,a[j]*xfit,'k--')

########################################################################
def med_plot(data,scan_pnt=0,hold=False,ylog=True):
    """
    plot med spectra
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
class ImageScan:
    """
    Class to hold a collection of meds associated with a scan
    ie one med per scan point
    """
    def __init__(self,med=[]):
        self.med = med

    ################################################################
    def med_ocr(self):
        """
        return outgoing count rate (ocr) values from meds
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
    def med_update_tau(self,tau):
        """
        update med tau factors
        """
        npnt = len(self.med)
        for j in range(npnt):
            self.med[j].update_correction(tau)
        
    ################################################################
    def med2xrf(self,xrf_params={},det_idx=0,emin=-1.,emax=-1.):
        """
        convert med objects to xrf objects
        """
        xrf = xrf_data.med2xrf(self.med,xrf_params=xrf_params,
                              lines = self.xrf_lines,
                              det_idx=det_idx,emin=emin,emax=emax)
        if xrf: self.xrf = xrf

########################################################################
########################################################################
if __name__ == '__main__':
    pass
