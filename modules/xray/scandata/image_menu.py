#######################################################################
"""
Tom Trainor (fftpt@uaf.edu)

Menu function to handle interactive
image processing

Modifications:
--------------


"""
#######################################################################

import types
import copy
import numpy as Num
import pylab

from   plotter import cursor
import image_data

########################################################################
def image_menu(data):
    """
    Interactively inspect images in scan data object
    """
    options = """
###################
Options:
Number of images  = %s
Current image     = %s
Current image roi = %s
1.  Display image
2.  Set roi from image zoom (Figure 1)
3.  Plot row/column sums (Figure 2)
4.  Select roi from sum plots (Figure 2)
5.  Apply current roi to all images
6.  Set num bgr points: %s
7.  Integrate current image
8.  Integrate all images
9.  Select scan point
10. Select next point 
11. Select previous point 
12. Remove current data point
13. Done
###################
"""
    prompt   = 'Select option (1)>'
    npts     = len(data.image)
    scan_pnt = 0
    roi      = []
    nbgr     = 4
    ret      = 0
    
    # check rois.
    if data.image_rois == None:
        data.image_rois = []
    if len(data.image_rois) != npts:
        data.image_rois = []
        for j in range(npts):
            data.image_rois.append([])

    # loop
    while ret != 13:
        ###
        try:
            roi = data.image_rois[scan_pnt]
        except:
            print "** error getting roi **"
            roi = []
        op = options % (str(npts),str(scan_pnt),str(roi),
                        str(nbgr))
        print op
        rret = raw_input(prompt)
        if not rret:
            ret = 1
        else:
            try:
                ret = int(rret)
            except:
                print "** Invalid option '%s' **" % rret
                ret = -1
        ####
        if ret == -1:
            pass
        elif ret == 1:
            figtitle = "Scan Point = %i, L = %6.3f" % (scan_pnt, data['L'][scan_pnt])
            image_data.image_plot(data.image[scan_pnt],fig=1,verbose=True,figtitle=figtitle)
        elif ret == 2:
            pylab.figure(1)
            (x1,x2,y1,y2) = pylab.axis()
            roi  = [int(x1),int(y1),int(x2),int(y2)]
            data.image_rois[scan_pnt] = roi
        elif ret == 3:
            roi = data.image_rois[scan_pnt]
            image = image_data.clip_image(data.image[scan_pnt],roi)
            image_data.sum_plot(image,nbgr=nbgr,fig=2)
        elif ret == 4:
            #roi = data.image_rois[scan_pnt]
            #image = image_data.clip_image(data.image[scan_pnt],roi)
            image = data.image[scan_pnt]
            image_data.sum_plot(image,nbgr=-1,fig=2)
            #c = pylab.cursor(fig=2)
            c = cursor(fig=2)
            (c1,y) = c.get_click(msg="Select left col sum")
            (c2,y) = c.get_click(msg="Select right col sum")
            (r1,y) = c.get_click(msg="Select left row sum")
            (r2,y) = c.get_click(msg="Select right row sum")
            roi = [int(c1),int(r1),int(c2),int(r2)]
            data.image_rois[scan_pnt] = roi
        elif ret == 5:
            roi = data.image_rois[scan_pnt]
            if len(roi) == 4:
                data.image_rois = []
                for j in range(npts):
                    data.image_rois.append(copy.copy(roi))
        elif ret == 6:
            p2 = "Num bgr points (%s)>" % str(nbgr)
            print p2
            n = raw_input('>>')
            if len(n) > 0:
                try: nbgr = int(n)
                except: pass
        elif ret == 7:
            data.integrate_image(idx=[scan_pnt],nbgr=nbgr,plot=True,fig=3)
        elif ret == 8:
            p2 = "Plot all images (False)>"
            print p2
            yn = raw_input('>>')
            if len(yn.strip()) == 0:
                yn = False
            else:
                yn = eval(yn)
            data.integrate_image(nbgr=nbgr,plot=yn)
            #
            pylab.figure(5, figsize = [5,4])
            pylab.clf()
            pylab.plot(data['L'],data['I_c'],'b',label='col sum')
            pylab.errorbar(data['L'],data['I_c'],data['Ierr_c'],fmt='o')
            pylab.plot(data['L'],data['I_r'],'g',label='row sum')
            pylab.errorbar(data['L'],data['I_r'],data['Ierr_r'],fmt='o')
            pylab.semilogy()
            pylab.legend()
            pylab.xlabel('L')
            pylab.ylabel('Integrated Intensity')
        elif ret == 9:
            p2 = "Enter scan point, max = %i (%i)>" % (npts-1,scan_pnt)
            print p2
            pnt = raw_input('>>')
            if len(pnt.strip()) == 0:
                pnt = scan_pnt
            else:
                pnt = int(pnt)
            if pnt > npts-1: pnt = npts-1
            if pnt < 0: pnt = 0
            scan_pnt = pnt
            figtitle = "Scan Point = %i, L = %6.3f" % (scan_pnt, data['L'][scan_pnt])
            image_data.image_plot(data.image[scan_pnt],fig=1,verbose=True,figtitle=figtitle)
        elif ret == 10:
            if scan_pnt + 1 < npts: 
                scan_pnt = scan_pnt + 1
                figtitle = "Scan Point = %i, L = %6.3f" % (scan_pnt, data['L'][scan_pnt])
                image_data.image_plot(data.image[scan_pnt],fig=1,verbose=True,figtitle=figtitle)
        elif ret == 11:
            if scan_pnt - 1 > -1: 
                scan_pnt = scan_pnt - 1
                figtitle = "Scan Point = %i, L = %6.3f" % (scan_pnt, data['L'][scan_pnt])
                image_data.image_plot(data.image[scan_pnt],fig=1,verbose=True,figtitle=figtitle)
        elif ret ==12:
            if data.image_peaks['I_c']== None:
                print "Remove not possible before first integration"
                pass
            else:
                p2 = "Are you sure you want to remove this Datapoint? (y/n)"
                print p2
                pnt = raw_input('>>')
                if pnt.lower() == 'y' or pnt.lower() =='yes':
                    data.image_peaks['I_c'][scan_pnt] = -1.
                    data.image_peaks['I_r'][scan_pnt] = -1.
                    data.image_peaks['Ierr_c'][scan_pnt] = 0.
                    data.image_peaks['Ierr_r'][scan_pnt] = 0.
                else:
                    print "Datapoint not removed" 
                    pass
        elif ret == 13:                
            print "All done"
        else:
            print "** Invalid option '%s' **" % ret
            pass

        # integrate all and plot
        """
        data.integrate_image(nbgr=nbgr,plot=False)
        pylab.figure(5, figsize = [5,4])
        pylab.clf()
        pylab.plot(data['L'],data['I_c'],'b',label='col sum')
        pylab.errorbar(data['L'],data['I_c'],data['Ierr_c'],fmt='o')
        pylab.plot(data['L'],data['I_r'],'g',label='row sum')
        pylab.errorbar(data['L'],data['I_r'],data['Ierr_r'],fmt='o')
        pylab.semilogy()
        """

################################################################################

