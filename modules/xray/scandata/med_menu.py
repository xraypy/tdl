#######################################################################
"""
T. Trainor (fftpt@uaf.edu)
Functions for interactive med 

Modifications:
--------------

"""
#######################################################################

import types
import numpy as Num
import deadtime


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

