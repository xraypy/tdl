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
import pylab
import deadtime

########################################################################
def fit_deadtime(data,x='io',y='Med',norm='Seconds',offset=True,display=True):
    """
    do a deadtime fit to data

    Note problem - we arent normalizing Io by time
    but ocr from the med is cps
    Also in correction plot, need yarr_corr/time etc
    
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
    pylab.clf()
    pylab.subplot(2,1,1)

    # data
    for yfit in yfit_arr:
        pylab.plot(xfit,yfit,'k.')
        
    # plot fit
    for j in range(len(tau)):
        if off != None:
            params = (tau[j],a[j],off[j])
            offset = True
        else:
            params = (tau[j],a[j])
            offset = False
        ycal = deadtime.calc_ocr(params,xfit,offset)
        pylab.plot(xfit,ycal,'r-')
        
    # plot corrected data
    pylab.subplot(2,1,2)
    for j in range(len(yfit_arr)):
        pylab.plot(xfit,yfit_arr[j],'k.')
        pylab.plot(xfit,ycorr_arr[j],'r-')
        if offset:
            pylab.plot(xfit,a[j]*xfit+off[j],'k--')
        else:
            pylab.plot(xfit,a[j]*xfit,'k--')

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
        pylab.clf()
    if nchan == 1:
        pylab.plot(en[0],ch[0])
    else:
        for j in range(nchan):
            if j not in bad:
                label = str(j)
                pylab.plot(en[j],ch[j]+tiny,label=label)
        pylab.legend()
    if ylog:
        pylab.semilogy()
        pylab.ylim(ymin=1.)

########################################################################

