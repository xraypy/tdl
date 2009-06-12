#######################################################################
"""
Tom Trainor (fftpt@uaf.edu)
Frank Heberling (Frank.Heberling@ine.fzk.de)
Matt Newville (newville@cars.uchicago.edu)

Methods to handle image processing
Reads in data for the pilatus detector
Simple integrations and plotting

Modifications:
--------------


"""
#######################################################################
"""
Todo
 - Need to improve show image / options
 - Improve plotting in integrate
 - Need non-linear backgrounds in integrate?


"""
#######################################################################
"""
Note in older versions of PIL:
 To read the pilatus, you must add this line:
    (1, 1, 1, (32,), ()): ("F", "F;32F"),
 to the OPEN_INFO dict in TiffImagePlugin.py (part of the Image module)

 starting at line 130 of TiffImagePlugin.py (from PIL 1.1.6):
   (1, 1, 2, (8,), ()): ("L", "L;R"),
   (1, 1, 1, (16,), ()): ("I;16", "I;16"),
   (1, 1, 1, (32,), ()): ("F", "F;32F"),
   (1, 2, 1, (16,), ()): ("I;16S", "I;16S"),
"""
#######################################################################

import types
import copy
import numpy as Num
import Image
import pylab

from  mathutil import LinReg
from  plotter import cursor

#######################################################################

def read_files(file_prefix,start=0,end=100,nfmt=3):
    """
    read files
    """
    images  = []
    format = '%' + str(nfmt) + '.' + str(nfmt) + 'd'
    for j in range(start,end+1):
        ext  = format % j
        file = file_prefix + '_' + ext + '.tif'
        arr  = read_file(file)
        images.append(arr)
    return images

def read_file(tiff_file):
    """
    read file
    """
    try:
        im  = Image.open(tiff_file)
    except:
        print "Error reading file: %s" % tiff_file
        return Num.array([[0]])
    arr = Num.fromstring(im.tostring(), dtype='int32')
    arr.shape = im.size[1],im.size[0]
    #return arr.transpose()
    return arr

def image_plot(img,fig=None,figtitle='',cmap=None,verbose=False,):
    """
    show image
    
    * arguments:
       img              # the image array to be displayed
       fig = None       # Figure to plot to
       figtitle = ''    # Title
       cmap = None      # Colormap.  None uses default
                        # you can pass a string name if its in pylab.cm.colormaps
                        # or you can pass explicitly the colormap
       verbose = False  # Print some fig statistics

    * examples:
       >>>image_plot(im,fig=1,figtitle='Image',cmap='hot')
       >>>image_plot(im,fig=1,figtitle='Image',cmap=pylab.cm.Spectral)
       
    """
    if verbose:
        print '###################'
        print 'Some statistics for plotted image'
        print 'Image total= ', img.sum()
        print 'Max value = ',  img.max()
        print 'Min value = ',  img.min()
    if fig != None:
        pylab.figure(fig)
        pylab.clf()

    # pylab.imshow(img, cmap = pylab.cm.hot)
    if cmap != None:
        if type(cmap) == types.StringType:
            if cmap in pylab.cm.cmapnames:
                cmap = getattr(pylab.cm,cmap)
            else:
                cmap = None
    pylab.imshow(img,cmap=cmap)
    pylab.colorbar(orientation='horizontal')

    if figtitle:
        pylab.title(figtitle, fontsize = 12)

############################################################################
def clip_image(image,roi=[]):
    """
    roi = [c1,r1,c2,r2]
    """
    if len(roi) != 4:
        roi = [0,0,image.shape[1], image.shape[0]]
    if roi[0] < roi[2]:
        c1 = roi[0]
        c2 = roi[2]
    else:
        c1 = roi[2]
        c2 = roi[0]
    if roi[1] < roi[3]:
        r1 = roi[1]
        r2 = roi[3]
    else:
        r1 = roi[3]
        r2 = roi[1]
    #(c1,r1,c2,r2) = roi
    return image[r1:r2, c1:c2]

def line_sum(image,sumflag='c',nbgr=5):
    """
    sum down 'c'olumns or across 'r'ows
    """
    #
    if sumflag == 'c':
        data = image.sum(axis=0)
    elif sumflag == 'r':
        data = image.sum(axis=1)
    npts     = len(data)
    data_idx = Num.arange(npts,dtype=float)
    #data_err  = data**(0.5)

    ### compute linear background
    bgr = []
    if nbgr > 0:
        xbgr = Num.arange(0,nbgr,1,dtype=float)
        xbgr = Num.append(xbgr,Num.arange(npts-nbgr,npts,1))
        ybgr = Num.array(data[0:nbgr],dtype=float)
        ybgr = Num.append(ybgr,data[npts-nbgr:])
        lr   = LinReg(xbgr,ybgr,plot=False)
        bgr  = lr.calc_y(data_idx)
    return (data, data_idx, bgr)

def sum_plot(image,nbgr=5,fig=None):
    #
    if fig != None:
        pylab.figure(fig)
        pylab.clf()
    else:
        pylab.figure()
    #
    pylab.subplot(211)
    pylab.title('Column Sum')
    (data, data_idx, bgr) = line_sum(image,sumflag='c',nbgr=nbgr)
    pylab.plot(data_idx, data, 'r')
    if nbgr > 0:
        pylab.plot(data_idx, bgr, 'b')
    #
    pylab.subplot(212)
    pylab.title('Row Sum')
    (data, data_idx, bgr) = line_sum(image,sumflag='r',nbgr=nbgr)
    pylab.plot(data_idx, data, 'r')
    if nbgr > 0:
        pylab.plot(data_idx, bgr, 'b')

################################################################################
class ImageAna:
    def __init__(self,image,roi=[],nbgr=5,plot=True,fig=None,figtitle=''):
        """
        Note take x as the horizontal image axis index (col index), 
        and y as vertical axis index (row index).
        Therefore, in image indexing:
          image[y1:y2,x1:x2] ==> rows y1 to y2 and cols x1 to x2

        Parameters:
          * roi  = [x1,y1,x2,y2]
          * nbgr = number of bgr points
          * show_plot = show fancy plot 

        """
        ### roi = [x1,y1,x2,y2] = [c1,r1,c2,r2]
        if len(roi) < 4:
            roi = [0,0,image.shape[1], image.shape[0]]
        if roi[0] < roi[2]:
            c1 = roi[0]
            c2 = roi[2]
        else:
            c1 = roi[2]
            c2 = roi[0]
        if roi[1] < roi[3]:
            r1 = roi[1]
            r2 = roi[3]
        else:
            r1 = roi[3]
            r2 = roi[1]
        self.roi    = (int(c1),int(r1),int(c2),int(r2))
        #
        self.image  = image
        self.nbgr   = nbgr
        self.title  = figtitle
        #
        self.I_c    = 0.0
        self.I_r    = 0.0
        self.Ibgr_c = 0.0
        self.Ibgr_r = 0.0
        self.Ierr_c = 0.0
        self.Ierr_r = 0.0
        #
        self.integrate_sums()
        # plot
        if plot == True:
            self.plot(fig=fig)
            
    ############################################################################
    def integrate_sums(self):
        """
        Convert image to line plots and integrate.

        Note approx error by assuming data and bgr errors are:
           sig_i = sqrt(I_i)
        Therefore:
          (sig_Itot)**2 = Sum(I_i)
          (sig_Itot)**2 = Sum(Ibgr_i)
          Ierr = sqrt((sig_Itot)**2 + (sig_Itot)**2)
        """
        def _integrate(data, data_idx, bgr):
            ### integrate
            #Itot = data.sum()
            #Ibgr = bgr.sum()
            Itot = 0.0
            Ibgr = 0.0
            Itot = Num.trapz(data)
            Ibgr = Num.trapz(bgr)
            I    = Itot - Ibgr
            ### compute errors
            #Ierr = (data.sum() + bgr.sum())**(0.5)
            Ierr = (Itot + Ibgr)**(0.5)
            return(I,Ierr,Ibgr)

        # clip image
        image = clip_image(self.image,self.roi)
        nbgr = self.nbgr

        # integrate col sum
        (data, data_idx, bgr) = line_sum(image,sumflag='c',nbgr=nbgr)
        (I,Ierr,Ibgr) = _integrate(data, data_idx, bgr)
        self.I_c    = I
        self.Ierr_c = Ierr
        self.Ibgr_c = Ibgr
        
        # integrate row sum
        (data, data_idx, bgr) = line_sum(image,sumflag='r',nbgr=nbgr)
        (I,Ierr,Ibgr) = _integrate(data, data_idx, bgr)
        self.I_r    = I
        self.Ierr_r = Ierr
        self.Ibgr_r = Ibgr

    ############################################################################
    def plot(self,fig=None):
        """
        make fancy 4-panel plot
        """
        colormap = pylab.cm.hot
        if fig != None:
            pylab.figure(fig)
            pylab.clf()
            pylab.figure(fig,figsize=[12,8])
        else:
            pylab.figure(fig,figsize=[12,8])
        title_c = 'I = %g, Ierr = %g, Ibgr = %g' % (self.I_c,self.Ierr_c,self.Ibgr_c)
        title_r = 'I = %g, Ierr = %g, Ibgr = %g' % (self.I_r,self.Ierr_r,self.Ibgr_r)
        #
        nbgr = self.nbgr
        
        # clipped image
        image = clip_image(self.image,self.roi)

        # full image with roi box
        (c1,r1,c2,r2) = self.roi
        bild= copy.copy(self.image)
        bild[r1-1:r1,c1:c2] = bild.max()
        bild[r2:r2+1,c1:c2] = bild.max()
        bild[r1:r2,c1-1:c1] = bild.max()
        bild[r1:r2,c2:c2+1] = bild.max()

        # plot column sum
        pylab.subplot(221)
        pylab.title(title_c, fontsize = 12)
        (data, data_idx, bgr) = line_sum(image,sumflag='c',nbgr=nbgr)
        pylab.plot(data_idx, data, 'r')
        #pylab.plot(data_idx, bgr, 'b')
        pylab.fill_between(data_idx,bgr,data,where= data>=bgr,facecolor='green')
        pylab.fill_between(data_idx,0,bgr,facecolor='blue')
        pylab.axis([0, data_idx.max(), 0, data.max()*1.05])

        # plot full image with ROI
        pylab.subplot(222)
        pylab.title(self.title, fontsize = 12)
        pylab.imshow(bild,cmap=colormap)
        pylab.colorbar(orientation='horizontal')

        # plot zoom on image 
        pylab.subplot(223)
        pylab.imshow(image, cmap=colormap, aspect='auto')

        # plot row sum
        pylab.subplot(224)
        pylab.title(title_r, fontsize = 12)
        (data, data_idx, bgr) = line_sum(image,sumflag='r',nbgr=nbgr)
        pylab.plot(data, data_idx, 'r')
        pylab.plot(bgr, data_idx, 'b')
        #pylab.fill_between(data, data_idx,bgr, where=bgr <= data, facecolor = 'green')
        pylab.axis([0,data.max()*1.05, data_idx.max(), 0])
        pylab.xticks(rotation=-45)

########################################################################
def image_integrate(data):
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
4.  Select roi from sum plots (Figure 3)
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
            image_plot(data.image[scan_pnt],fig=1,verbose=True,figtitle=figtitle)
        elif ret == 2:
            pylab.figure(1)
            (x1,x2,y1,y2) = pylab.axis()
            roi  = [int(x1),int(y1),int(x2),int(y2)]
            data.image_rois[scan_pnt] = roi
        elif ret == 3:
            roi = data.image_rois[scan_pnt]
            image = clip_image(data.image[scan_pnt],roi)
            sum_plot(image,nbgr=nbgr,fig=2)
        elif ret == 4:
            #roi = data.image_rois[scan_pnt]
            #image = clip_image(data.image[scan_pnt],roi)
            image = data.image[scan_pnt]
            sum_plot(image,nbgr=-1,fig=2)
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
            image_plot(data.image[scan_pnt],fig=1,verbose=True,figtitle=figtitle)
        elif ret == 10:
            if scan_pnt + 1 < npts: 
                scan_pnt = scan_pnt + 1
                figtitle = "Scan Point = %i, L = %6.3f" % (scan_pnt, data['L'][scan_pnt])
                image_plot(data.image[scan_pnt],fig=1,verbose=True,figtitle=figtitle)
        elif ret == 11:
            if scan_pnt - 1 > -1: 
                scan_pnt = scan_pnt - 1
                figtitle = "Scan Point = %i, L = %6.3f" % (scan_pnt, data['L'][scan_pnt])
                image_plot(data.image[scan_pnt],fig=1,verbose=True,figtitle=figtitle)
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
################################################################################
if __name__ == '__main__':
    import sys
    try:
        fname = sys.argv[1]
    except:
        print 'read_pilatus  tiff file'
        sys.exit()
    a = read_file(fname)
    image_show(a)
