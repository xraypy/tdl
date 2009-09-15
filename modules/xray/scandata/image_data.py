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
 - Improve image background determination
 - In ImageAna.plot should we plot row and column sums
   after bgr subtraction?
   
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
import numpy as num
import Image
import pylab

from  mathutil import LinReg
import xrf_bgr

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
        return num.array([[0]])
    arr = num.fromstring(im.tostring(), dtype='int32')
    arr.shape = im.size[1],im.size[0]
    #return arr.transpose()
    return arr

############################################################################
def image_plot(img,fig=None,figtitle='',cmap=None,verbose=False,Im_max = None):
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
    pylab.imshow(img,cmap=cmap, vmax = Im_max)
    pylab.colorbar(orientation='horizontal')

    if figtitle:
        pylab.title(figtitle, fontsize = 12)

############################################################################
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

############################################################################
def clip_image(image,roi=[],cp=False):
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
    if cp == True:
        return copy.copy(image[r1:r2, c1:c2])
    else:
        return image[r1:r2, c1:c2]

############################################################################
def line_sum(image,sumflag='c',nbgr=3):
    """
    sum down 'c'olumns or across 'r'ows
    this returns the summed data and linear background
    """
    if sumflag == 'c':
        data = image.sum(axis=0)
    elif sumflag == 'r':
        data = image.sum(axis=1)
    npts     = len(data)
    data_idx = num.arange(npts,dtype=float)
    #data_err  = data**(0.5)

    ### compute linear background
    bgr = num.zeros(npts)
    if nbgr > 0:
        xbgr = num.arange(0,nbgr,1,dtype=float)
        xbgr = num.append(xbgr,num.arange(npts-nbgr,npts,1))
        ybgr = num.array(data[0:nbgr],dtype=float)
        ybgr = num.append(ybgr,data[npts-nbgr:])
        lr   = LinReg(xbgr,ybgr,plot=False)
        bgr  = lr.calc_y(data_idx)
    return (data, data_idx, bgr)

############################################################################
def line_sum_integral(image,sumflag='c',nbgr=3):
    """
    calc the integral after image is summed down 'c'olumns
    or across 'r'ows.  This uses a linear background
    """
    # get line sum
    (data, data_idx, bgr) = line_sum(image,sumflag=sumflag,nbgr=nbgr)

    ### integrate
    #Itot = data.sum()
    #Ibgr = bgr.sum()
    Itot = 0.0
    Ibgr = 0.0
    Itot = num.trapz(data)
    Ibgr = num.trapz(bgr)
    I    = Itot - Ibgr
    ### compute errors
    #Ierr = (data.sum() + bgr.sum())**(0.5)
    Ierr = (Itot + Ibgr)**(0.5)
    return(I,Ierr,Ibgr)

############################################################################
def image_bgr(image,cwidth=100,rwidth=100,plot=False):
    """
    calculate a background for the image.
    we assume that image is already clipped
    and that it consists of a 'single peak'
    and a non-lin bgr.

    cwidth and rwidth are the 'widths' to be used
    for column and row sums respectivley. The cwidth
    value should correspond roughly to the actual peak
    width along the column direction (or the peak width in
    the row sum).  Visa versa for rrwidth.  

    Therefore the background should fit
    features that are in general broader than these values

    initially try using the xrf bgr algorithm
    if this works then we can optimize for images

    """
    bgr_arr_r = num.zeros(image.shape)
    bgr_arr_c = num.zeros(image.shape)

    # fit to rows
    if rwidth > 0:
        b = xrf_bgr.Background(bottom_width=rwidth,compress=1)
        for j in range(image.shape[0]):
            b.calc(image[j])
            bgr_arr_r[j,:] = copy.copy(b.bgr)

    # fit to cols
    if cwidth > 0:
        b = xrf_bgr.Background(bottom_width=cwidth,compress=1)
        for j in range(image.shape[1]):
            b.calc(image[:,j])
            bgr_arr_c[:,j] = copy.copy(b.bgr)

    # combine the two bgrs 
    """
    # combine by taking the minimum of the
    # two at each point ??
    bgr = bgr_arr_r
    idx = num.where(bgr > bgr_arr_c)
    bgr[idx] = bgr_arr_c[idx]
    """
    # take the average of the 2 bgrs
    if (rwidth > 0) and (cwidth > 0):
        bgr = (bgr_arr_r + bgr_arr_c)/2.
    elif rwidth > 0:
        bgr = bgr_arr_r
    elif cwidth > 0:
        bgr = bgr_arr_c
    else:
        bgr = num.zeros(image.shape)

    #show
    if plot:
        pylab.figure(3)
        pylab.clf()
        pylab.subplot(3,1,1)
        pylab.imshow(image)
        pylab.title("image")
        pylab.colorbar()

        pylab.subplot(3,1,2)
        pylab.imshow(bgr)
        pylab.title("background")
        pylab.colorbar()

        pylab.subplot(3,1,3)
        pylab.imshow(image-bgr)
        pylab.title("image - background")
        pylab.colorbar()

    return bgr

def _bgr(y,width):
    """
    calc background for a given line
    simplified form of the algorithm used for
    XRF bgr
    """

    # first smooth and rebin y
    # at end restore length with spline
    npts = len(y)
    
    # here calc polynomial 
    pow     = 0.5
    npoly   = int(2*width) + 1
    pdelx   = num.array(range(npoly),dtype=float) - float((npoly-1)/2)
    poly    = (width**2. - pdelx**2)**pow  - (width**2.)**pow

    # loop through each point
    delta = num.zeros(len(poly))
    n = (npoly-1)/2
    for j in range(npts):
        lidx  = num.max((0,j-n))
        ridx  = num.min((npts,j+n))
        print n, lidx, ridx
        delta = y[lidx:ridx] - (y[j] + poly[])
    return delta

################################################################################
class ImageAna:
    def __init__(self,image,roi=[],nbgr=5,cwidth=0,rwidth=0,
                 plot=True,fig=None,figtitle=''):
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
        self.I      = 0.0
        self.Ibgr   = 0.0
        self.Ierr   = 0.0
        #
        self.I_c    = 0.0
        self.I_r    = 0.0
        self.Ibgr_c = 0.0
        self.Ibgr_r = 0.0
        self.Ierr_c = 0.0
        self.Ierr_r = 0.0
        #
        self.bgr_cwidth = cwidth
        self.bgr_rwidth = rwidth
        #
        self.integrate()
        # plot
        if plot == True: self.plot(fig=fig)
            
    ############################################################################
    def integrate(self):
        """
        Integrate Image.

        Note approx error by assuming data and bgr errors are:
           sig_i = sqrt(I_i)
        Therefore:
          (sig_Itot)**2 = Sum(I_i)
          (sig_Ibgr)**2 = Sum(Ibgr_i)
          Ierr = sqrt((sig_Itot)**2 + (sig_Ibgr)**2)
        """

        # clip image
        image = clip_image(self.image,self.roi)

        # integrate the roi
        self.I = num.sum(num.trapz(image))

        # subtract image background
        if (self.bgr_cwidth > 0) or (self.bgr_rwidth > 0):
            bgr = image_bgr(image,cwidth=self.bgr_cwidth,
                            rwidth=self.bgr_rwidth)
            #image = image - bgr
            #self.I = num.sum(num.trapz(image))
            self.Ibgr = num.sum(num.trapz(bgr))
            self.I    = self.I - self.Ibgr

        # error
        self.Ierr = (self.I + self.Ibgr)**0.5

        # the below integrations should only be done
        # if we have a flag set to do so??
        # note these are done on the non background
        # subtracted image
        
        # integrate col sum
        (I,Ierr,Ibgr) = line_sum_integral(image,sumflag='c',nbgr=self.nbgr)
        self.I_c      = I
        self.Ierr_c   = Ierr
        self.Ibgr_c   = Ibgr
        
        # integrate row sum
        (I,Ierr,Ibgr) = line_sum_integral(image,sumflag='r',nbgr=self.nbgr)
        self.I_r      = I
        self.Ierr_r   = Ierr
        self.Ibgr_r   = Ibgr

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
        title_c = 'I_c = %g, Ierr_c = %g, Ibgr_c = %g' % (self.I_c,self.Ierr_c,self.Ibgr_c)
        title_r = 'I_r = %g, Ierr_r = %g, Ibgr_r = %g' % (self.I_r,self.Ierr_r,self.Ibgr_r)
        title_roi = 'I = %g, Ierr = %g, Ibgr = %g' % (self.I,self.Ierr,self.Ibgr)
        #
        nbgr = self.nbgr
        
        # clipped image
        image = clip_image(self.image,self.roi)

        # background
        if (self.bgr_cwidth > 0) or (self.bgr_rwidth > 0):
            ibgr = image_bgr(image,cwidth=self.bgr_cwidth,
                             rwidth=self.bgr_rwidth)
            #image = image - bgr
            title_roi = title_roi + '\n(background subtracted)'
        else:
            ibgr = None

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
        pylab.plot(data_idx, bgr, 'b')
        pylab.axis([0, data_idx.max(), 0, data.max()*1.05])
        #pylab.fill_between(data_idx,bgr,data,where= data>=bgr,facecolor='green')
        #pylab.fill_between(data_idx,0,bgr,facecolor='blue')
        if ibgr != None:
            (data, data_idx, bgr) = line_sum(image-ibgr,sumflag='c',nbgr=0)
            pylab.plot(data_idx, data, 'k-.')
            (data, data_idx, bgr) = line_sum(ibgr,sumflag='c',nbgr=0)
            pylab.plot(data_idx, data, 'g-.')

        # plot full image with ROI
        pylab.subplot(222)
        pylab.title(self.title, fontsize = 12)
        pylab.imshow(bild,cmap=colormap)
        pylab.colorbar(orientation='horizontal')

        # plot zoom on image 
        pylab.subplot(223)
        pylab.title(title_roi,fontsize = 12)
        if ibgr != None:
            pylab.imshow(image - ibgr, cmap=colormap, aspect='auto')
        else:
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
        if ibgr != None:
            (data, data_idx, bgr) = line_sum(image-ibgr,sumflag='r',nbgr=0)
            pylab.plot(data, data_idx, 'k-.')
            (data, data_idx, bgr) = line_sum(ibgr,sumflag='r',nbgr=0)
            pylab.plot(data, data_idx, 'g-.')

################################################################################
################################################################################
if __name__ == '__main__':
    """
    import sys
    try:
        fname = sys.argv[1]
    except:
        print 'read_pilatus  tiff file'
        sys.exit()
    a = read_file(fname)
    image_show(a)
    """
    import pylab
    x = num.array(range(30))
    y = 20.+num.exp( -1*((x-15.)**2.)/20.)
    #pylab.plot(x,y)
    #
    p = _bgr(y,10)
    print p
    pylab.newplot(p,'o-')

    