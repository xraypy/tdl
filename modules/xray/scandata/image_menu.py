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
import numpy as num
import pylab

from   pds.shellutil import Menu, get_yn, get_int, get_flt, 
from   pds.shellutil import get_tf, show_more
from   plotter import cursor
import image_data

########################################################################
IMG_HEADER = """
Number of images  = %s
Current image     = %s
Current image roi = %s
Current image rotation angle = %s
"""

IMG_LABELS = ['display','imax','rotangle','zoomroi','plotsums',
              'selectroi','bgr','copyall','integrate','intall',
              'point','next','previous','flag','help','done']
IMG_DESCR = ["Display image",
             "Set max image intensity value",
             "Set image rotation angle (deg ccw)",
             "Set roi from image zoom (Figure 1)",
             "Plot row/column sums (Figure 2)",
             "Select roi from sum plots (Figure 2)",
             "Set background parameters",
             "Apply current roi and background params to all images",
             "Integrate current image",
             "Integrate all images",
             "Select scan point",
             "Select next point ",
             "Select previous point", 
             "Flag as bad point",
             "Show options",
             "Done"]

########################################################################
def image_menu(data):
    """
    Interactively inspect images in scan data object
    """
    prompt   = 'Select option >'
    npts     = len(data.image)
    scan_pnt = 0
    roi      = []
    ret      = ''
    im_max   = None

    # local plot fun
    def _implot(data,scan_pnt,im_max):
        rotangle = data.image_rotangle[scan_pnt]
        figtitle = "Scan Point = %i, L = %6.3f" % (scan_pnt, data['L'][scan_pnt])
        image_data.image_plot(data.image[scan_pnt],fig=1,verbose=True,figtitle=figtitle,
                              im_max=im_max,rotangle=rotangle)

    # check init and plot first
    if len(data.image_peaks) != npts:
        data._init_image()
    _implot(data,scan_pnt,im_max)
    
    # make menu
    m = Menu(labels=IMG_LABELS,descr=IMG_DESCR,sort=False,matchidx=True)
    
    # loop
    while ret != 'done':
        roi      = data.image_rois[scan_pnt]
        rotangle = data.image_rotangle[scan_pnt]
        header   = IMG_HEADER % (str(npts),str(scan_pnt),str(roi),str(rotangle))
        m.header = header
        ret      = m.prompt(prompt)

        if ret == 'display':
            _implot(data,scan_pnt,im_max)
        elif ret == 'imax':
            print "Enter maximum intensity value to apply to image"
            im_max = raw_input('>>')
            if len(im_max.strip()) == 0:
                im_max = None
            else:
                im_max = int(im_max)
            _implot(data,scan_pnt,im_max)
        elif ret == "rotangle":
            rotangle = get_flt(p='Enter rotation angle in degrees ccw>',
                               d=0.0,min=-360.,max=360.)
            data.image_rotangle[scan_pnt] = rotangle
            _implot(data,scan_pnt,im_max)
        elif ret == 'zoomroi':
            pylab.figure(1)
            (x1,x2,y1,y2) = pylab.axis()
            roi  = [int(x1),int(y1),int(x2),int(y2)]
            data.image_rois[scan_pnt] = roi
        elif ret == 'plotsums':
            roi      = data.image_rois[scan_pnt]
            rotangle = data.image_rotangle[scan_pnt]
            image    = image_data.clip_image(data.image[scan_pnt],roi,rotangle=rotangle)
            bgr_par  = data.image_bgrpar[scan_pnt]
            image_data.sum_plot(image,fig=2,**bgr_par)
        elif ret == 'selectroi':
            image = data.image[scan_pnt]
            bgr_par = data.image_bgrpar[scan_pnt]
            image_data.sum_plot(image,fig=2,**bgr_par)
            c = cursor(fig=2)
            (c1,y) = c.get_click(msg="Select left col sum")
            (c2,y) = c.get_click(msg="Select right col sum")
            (r1,y) = c.get_click(msg="Select left row sum")
            (r2,y) = c.get_click(msg="Select right row sum")
            roi = [int(c1),int(r1),int(c2),int(r2)]
            data.image_rois[scan_pnt] = roi
        elif ret == 'bgr':
            bgr_par = data.image_bgrpar[scan_pnt]
            bgr_par = bgr_menu(bgr_par)
            data.image_bgrpar[scan_pnt] = bgr_par
        elif ret == 'copyall':
            roi      = data.image_rois[scan_pnt]
            rotangle = data.image_rotangle[scan_pnt]
            bgr_par  = data.image_bgrpar[scan_pnt]
            #
            data.image_rois = []
            data.image_rotangle = []
            data.image_bgrpar = []
            #
            for j in range(npts):
                data.image_rois.append(copy.copy(roi))
                data.image_rotangle.append(rotangle)
                data.image_bgrpar.append(copy.copy(bgr_par))
        elif ret == 'integrate':
            data.integrate_image(idx=[scan_pnt],
                                 plot=True,fig=3)
        elif ret == 'intall':
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
            pylab.legend(loc = 9)
            pylab.xlabel('L')
            pylab.ylabel('Integrated Intensity')
        elif ret == 'point':
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
            _implot(data,scan_pnt,im_max)
        elif ret == 'next':
            if scan_pnt + 1 < npts: 
                scan_pnt = scan_pnt + 1
                _implot(data,scan_pnt,im_max)
        elif ret == 'previous':
            if scan_pnt - 1 > -1: 
                scan_pnt = scan_pnt - 1
                _implot(data,scan_pnt,im_max)
        elif ret == 'flag':
            if int(scan_pnt) in data.bad_points:
                data.bad_points.remove(scan_pnt)
                print "Data point removed from bad list"
            else:
                data.bad_points.append(int(scan_pnt))
                print "Data point added to bad list"
        else:
            pass


########################################################################
BGR_INFO = """
################################################
* bgrflag is flag for how to do backgrounds:
   = 0 determine row and column backgrounds after summation
   = 1 determine 2D background using 'c'olumn direction 
   = 2 determine 2D background using 'r'ow direction 
   = 3 determine 2D background using average of the
       'r'ow and 'c'olumn directions 

-> below params are for 'c'olumn and 'r'ow directions

* c/rnbgr = number of end points to use in linear background determination
  (see background.background)
      
* c/rwidth should correspond roughly to the actual peak
  widths. The background funciton should fit
  features that are in general broader than these values
  Note estimate cwidth using width of peak in row sum
  and rwidth using the width of the peak in the col sum.
  Note that width = 0 corresponds to no polynomial bgr
  
* c/rpow is the power of the polynomial used in background determination
  (see background.background)
  
* c/rtangent is a flag to indicate if local slope of the data should be fitted 
  (see background.background)
################################################
"""

BGR_HEADER = """
* bgrflag=%i,
* cnbgr=%i, cwidth=%g,cpow=%g,ctan=%s
* rnbgr=%i, rwidth=%g,rpow=%g,rtan=%s
"""

BGR_LABELS = ['help','info','bgrflag',
              'cnbgr','cwidth','cpow','ctan',
              'rnbgr','rwidth','rpow','rtan',
              'done']

BGR_DESCR = ["Show options","Get more info on parameter defintions", "Set background flag",
             "Set column sum num bgr for linear background",
             "Set column sum peak width for non-linear background",
             "Set column sum polynomial power for non-linear background",
             "Set column sum tangent flag (True or False)",
             "Set column row num bgr for linear background",
             "Set column row peak width for non-linear background",
             "Set column row polynomial power for non-linear background",
             "Set column row tangent flag (True or False)",
             "All done"]

IMG_BGR_PARAMS = {'bgrflag':0,
                  'cnbgr':5,'cwidth':0,'cpow':2.,'ctan':False,
                  'rnbgr':5,'rwidth':0,'rpow':2.,'rtan':False}

def bgr_menu(bgr_params=IMG_BGR_PARAMS):
    """
    Get bgr options
    """
    prompt   = 'Select option >'

    # make menu
    m   = Menu(labels=BGR_LABELS,descr=BGR_DESCR,sort=False,matchidx=True)
    ret = ''
    
    while ret != 'done':
        header = BGR_HEADER % (bgr_params['bgrflag'],
                               bgr_params['cnbgr'],bgr_params['cwidth'],
                               bgr_params['cpow'],str(bgr_params['ctan']),
                               bgr_params['rnbgr'],bgr_params['rwidth'],
                               bgr_params['rpow'],str(bgr_params['rtan']))
        m.header = header
        ret      = m.prompt(prompt)
        #
        if ret == 'bgrflag':
            bgr_params['bgrflag'] = get_int(p='Enter bgrflag>',
                                            d=bgr_params['bgrflag'],o=[0,1,2,3])
        #
        elif ret == 'cnbgr':
            bgr_params['cnbgr'] = get_int(p='Enter col nbgr>',
                                          d=bgr_params['cnbgr'])
        elif ret == 'cwidth':
            bgr_params['cwidth'] = get_int(p='Enter col width>',
                                           d=bgr_params['cwidth'])
        elif ret == 'cpow':
            bgr_params['cpow'] = get_flt(p='Enter col pow>',
                                         d=bgr_params['cpow'],min=0.,max=5.)
        elif ret == 'ctan':
            bgr_params['ctan'] = get_tf(p='Enter col tan flag>',
                                           d=bgr_params['ctan'])
        #
        elif ret == 'rnbgr':
            bgr_params['rnbgr'] = get_int(p='Enter row nbgr>',
                                          d=bgr_params['rnbgr'])
        elif ret == 'rwidth':
            bgr_params['rwidth'] = get_int(p='Enter row width>',
                                           d=bgr_params['rwidth'])
        elif ret == 'rpow':
            bgr_params['rpow'] = get_flt(p='Enter row pow>',
                                         d=bgr_params['rpow'],min=0.,max=5.)
        elif ret == 'rtan':
            bgr_params['rtan'] = get_tf(p='Enter row tan flag>',
                                           d=bgr_params['rtan'])
        elif ret == 'info':
            show_more(BGR_INFO)
    return bgr_params
        

################################################################################

