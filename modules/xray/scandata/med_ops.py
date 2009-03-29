#######################################################################
"""
T. Trainor (fftpt@uaf.edu)
Functions for operating on / working with
med's in ScanData

Modifications:
--------------

"""
#######################################################################

import types
import numpy as Num

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
        ycorr_arr   = Num.zeros((ndet,npts))
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
    import pylab
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
def med_inspect(data):
    """
    Interactively inspect med's
    """
    options = """
###################
Options:
(Number of meds = %s)
1.  Display med
2.  Toggle 'total':    %s
3.  Toggle 'align':    %s
4.  Toggle 'correct':  %s
5.  Set bad_mca's:     %s
6.  Select scan point: %i
7.  Toggle ylog:       %s
8.  Toggle hold:       %s
9.  Plot
10. Med2Xrf
11. Break
###################
"""
    import pylab
    prompt = 'Select option (1)>'
    scan_pnt = 0
    max_pnt = len(data.med)
    ylog = True
    hold = False
    ret = 0
    while ret != 11:
        op = options % (str(data.med[scan_pnt].n_detectors),
                        str(data.med[scan_pnt].total),
                        str(data.med[scan_pnt].align),
                        str(data.med[scan_pnt].correct),
                        str(data.med[scan_pnt].bad_mca_idx),
                        scan_pnt, str(ylog),str(hold))
        print op
        ret = raw_input(prompt)
        if not ret:
            ret = 1
        else:
            ret = int(ret)
        ####
        if ret == 1:
            print "###################"
            print data.med[scan_pnt]
        elif ret == 2:
            for j in range(max_pnt):
                data.med[j].total = not data.med[j].total 
        elif ret == 3:
            for j in range(max_pnt):
                data.med[j].align = not data.med[j].align 
        elif ret == 4:
            for j in range(max_pnt):
                data.med[j].correct = not data.med[j].correct 
        elif ret == 5:
            default = str(data.med[scan_pnt].bad_mca_idx)
            p2 = "Enter bad idx(%s)>" % default
            bad = raw_input(p2)
            if len(bad.strip()) == 0:
                bad = default
            bad = eval(bad)
            for j in range(max_pnt):
                data.med[j].bad_mca_idx = bad
        elif ret == 6:
            p2 = "Enter scan point, max = %i (%i)>" % (max_pnt-1,scan_pnt)
            pnt = raw_input(p2)
            if len(pnt.strip()) == 0:
                pnt = scan_pnt
            else:
                pnt = int(pnt)
            if pnt > max_pnt-1: pnt = max_pnt-1
            if pnt < 0: pnt = 0
            scan_pnt = pnt
        elif ret == 7:
            ylog = not ylog
            if ylog:
                pylab.semilogy()
                pylab.ylim(ymin=1.)
        elif ret == 8:
            hold = not hold
        elif ret == 9:
            med_plot(data,scan_pnt=scan_pnt,hold=hold,ylog=ylog)
        elif ret == 10:
            data.med2xrf()
        elif ret == 11:
            print "All done"
        else:
            print "Unknown option %s" % ret
    
########################################################################
def med_plot(data,scan_pnt=0,hold=False,ylog=True):
    """
    plot med spectra
    """
    import pylab
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

