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
OPTIONS = """
###################
Options:
Number of images  = %s
Current image     = %s
Current image roi = %s
1.  Display image
2.  Set roi from image zoom (Figure 1)
3.  Plot row/column sums (Figure 2)
4.  Select roi from sum plots (Figure 2)
5.  Set background parameters
6.  Apply current roi and background params to all images
7.  Integrate current image
8.  Integrate all images
9.  Select scan point
10. Select next point 
11. Select previous point 
12. Flag as bad point
13. Set max image intensity value
14. Done
###################
"""

########################################################################
def image_menu(data):
    """
    Interactively inspect images in scan data object
    """
    prompt   = 'Select option (1)>'
    npts     = len(data.image)
    scan_pnt = 0
    roi      = []
    ret = 0

    # check init
    if len(data.image_peaks) != npts:
        data._init_image()

    # loop
    while ret != 14:
        ###
        roi  = data.image_rois[scan_pnt]
        op = OPTIONS % (str(npts),str(scan_pnt),str(roi))
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
            image_data.sum_plot(image,nbgr=bgr_params['nbgr'],fig=2)
        elif ret == 4:
            image = data.image[scan_pnt]
            image_data.sum_plot(image,nbgr=0,fig=2)
            c = cursor(fig=2)
            (c1,y) = c.get_click(msg="Select left col sum")
            (c2,y) = c.get_click(msg="Select right col sum")
            (r1,y) = c.get_click(msg="Select left row sum")
            (r2,y) = c.get_click(msg="Select right row sum")
            roi = [int(c1),int(r1),int(c2),int(r2)]
            data.image_rois[scan_pnt] = roi
        elif ret == 5:
            # should make this a seperate menu function
            # the text could then go into a 'h'elp menu
            bgr_params = data.image_bgrpar[scan_pnt]

            print "Note that linear backgrounds are applied to row"
            print "and column sum integrals. These sums apply to the"
            print "non-background substracted image roi"
            p2 = "Num bgr points for linear background (0 for no bgr) (%s)>" % str(bgr_params['nbgr'])
            print p2
            tmp = raw_input('>>')
            if len(tmp.strip()) > 0:
                data.image_bgrpar[scan_pnt]['nbgr'] = int(tmp)
                
            print "The below widths are used for determining the image"
            print "roi background.  If both are zero no background is "
            print "applied.  If only one is non-zero the background is"
            print "only calculated in the given direction (row or column)"
            print "If both are non-zero the background is taken as the"
            print "average of the background determined for each direction"
            p2 = "Enter peak width in column direction, cwidth (%s)>" % bgr_params['cwidth']
            print p2
            tmp = raw_input('>>')
            if len(tmp.strip()) > 0:
                data.image_bgrpar[scan_pnt]['cwidth'] = int(tmp)
            p2 = "Enter peak width in row direction, rwidth (%s)>" % bgr_params['rwidth']
            print p2
            tmp = raw_input('>>')
            if len(tmp.strip()) > 0:
                data.image_bgrpar[scan_pnt]['rwidth'] = int(tmp)
            print "bgr params:", data.image_bgrpar[scan_pnt]
            
            #plot bgr
            if  (data.image_bgrpar[scan_pnt]['rwidth'] > 0) or (data.image_bgrpar[scan_pnt]['cwidth'] > 0):
                img = image_data.clip_image(data.image[scan_pnt],
                                            roi=data.image_rois[scan_pnt])
                image_data.image_bgr(img,
                                     cwidth=data.image_bgrpar[scan_pnt]['cwidth'],
                                     rwidth=data.image_bgrpar[scan_pnt]['rwidth'],
                                     plot=True)
        elif ret == 6:
            roi = data.image_rois[scan_pnt]
            data.image_rois = []
            for j in range(npts):
                data.image_rois.append(copy.copy(roi))
            
            bgr_par = data.image_bgrpar[scan_pnt]
            data.image_bgrpar = []
            for j in range(npts):
                data.image_bgrpar.append(copy.copy(bgr_par))
        elif ret == 7:
            data.integrate_image(idx=[scan_pnt],
                                 plot=True,fig=3)
        elif ret == 8:
            p2 = "Plot all images (False)>"
            print p2
            yn = raw_input('>>')
            if len(yn.strip()) == 0:
                yn = False
            else:
                yn = eval(yn)
            data.integrate_image(plot=yn)
            #
            pylab.figure(5, figsize = [5,4])
            pylab.clf()
            #
            pylab.plot(data['L'],data['I'],'b',label='image sum')
            pylab.errorbar(data['L'],data['I'],data['Ierr'],fmt='bo')
            #
            pylab.plot(data['L'],data['I_c'],'r',label='col sum')
            pylab.errorbar(data['L'],data['I_c'],data['Ierr_c'],fmt='ro')
            #
            pylab.plot(data['L'],data['I_r'],'g',label='row sum')
            pylab.errorbar(data['L'],data['I_r'],data['Ierr_r'],fmt='go')
            #
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
            if int(scan_pnt) in data.bad_points:
                data.bad_points.remove(scan_pnt)
                print "Data point removed from bad list"
            else:
                data.bad_points.append(int(scan_pnt))
                print "Data point added to bad list"
        elif ret == 13:
            print "Enter maximum intensity value to apply to image"
            Im_max = raw_input('>>')
            if len(Im_max.strip()) == 0:
                Im_max = None
            else:
                Im_max = int(Im_max)
            figtitle = "Scan Point = %i, L = %6.3f" % (scan_pnt, data['L'][scan_pnt])
            image_data.image_plot(data.image[scan_pnt],fig=1,verbose=True,figtitle=figtitle, Im_max = Im_max)
        elif ret == 14:                
            print "All done"
        else:
            print "** Invalid option '%s' **" % ret
            pass

################################################################################

