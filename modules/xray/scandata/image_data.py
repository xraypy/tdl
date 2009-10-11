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
import pylab
#import Image
from Image import open as imopen
from scipy import ndimage

from mpcutils.mathutil import LinReg
from background import background

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
        #im  = Image.open(tiff_file)
        im  = imopen(tiff_file)
    except:
        print "Error reading file: %s" % tiff_file
        return num.array([[0]])
    arr = num.fromstring(im.tostring(), dtype='int32')
    arr.shape = im.size[1],im.size[0]
    #return arr.transpose()
    return arr

############################################################################
def clip_image(image,roi=[],rotangle=0.0,cp=False):
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
        
    # rotate
    if rotangle != 0:
        image = ndimage.rotate(image,rotangle)    

    if cp == True:
        return copy.copy(image[r1:r2, c1:c2])
    else:
        return image[r1:r2, c1:c2]

############################################################################
def image_plot(img,fig=None,figtitle='',cmap=None,verbose=False,
               im_max=None,rotangle=0.0):
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
       im_max  = None   # Max intensity value
       rotangle = 0.0   # Rotation angle in degrees ccw

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
        print '###################'
    if fig != None:
        pylab.figure(fig)
        pylab.clf()

    # rotate
    if rotangle != 0:
        img = ndimage.rotate(img,rotangle)

    # pylab.imshow(img, cmap = pylab.cm.hot)
    if cmap != None:
        if type(cmap) == types.StringType:
            if cmap in pylab.cm.cmapnames:
                cmap = getattr(pylab.cm,cmap)
            else:
                cmap = None
    pylab.imshow(img,cmap=cmap, vmax = im_max)
    pylab.colorbar(orientation='horizontal')

    if figtitle:
        pylab.title(figtitle, fontsize = 12)

############################################################################
def sum_plot(image,bgrflag=0,
             cnbgr=5,cwidth=0,cpow=2.,ctan=False,
             rnbgr=5,rwidth=0,rpow=2.,rtan=False,
             fig=None):
    """
    plot sums with background.

    Note should we calc the bgr according
    to bgrflag???
    """
    #
    if fig != None:
        pylab.figure(fig)
        pylab.clf()
    else:
        pylab.figure()

    # col sum
    pylab.subplot(211)
    pylab.title('Column Sum')
    (data, data_idx, bgr) = line_sum(image,sumflag='c',nbgr=cnbgr,
                                     width=cwidth,pow=cpow,tangent=ctan)
    pylab.plot(data_idx, data, 'r')
    pylab.plot(data_idx, bgr, 'b')

    # row sum
    pylab.subplot(212)
    pylab.title('Row Sum')
    (data, data_idx, bgr) = line_sum(image,sumflag='r',nbgr=rnbgr,
                                     width=rwidth,pow=rpow,tangent=rtan)
    pylab.plot(data_idx, data, 'r')
    pylab.plot(data_idx, bgr, 'b')

############################################################################
def line_sum(image,sumflag='c',nbgr=0,width=0,pow=2.,tangent=False):
    """
    sum down 'c'olumns or across 'r'ows
    this returns the summed data and background
    """
    if sumflag == 'c':
        data = image.sum(axis=0)
    elif sumflag == 'r':
        data = image.sum(axis=1)
    npts     = len(data)
    data_idx = num.arange(npts,dtype=float)
    #data_err  = data**(0.5)

    ### compute background
    bgr = background(data,nbgr=nbgr,width=width,pow=pow,tangent=tangent)
    
    return (data, data_idx, bgr)

############################################################################
def line_sum_integral(image,sumflag='c',nbgr=0,width=0,pow=2.,
                      tangent=False):
    """
    Calc the integral after image is summed down 'c'olumns
    or across 'r'ows.  
    """
    # get line sum
    (data, data_idx, bgr) = line_sum(image,sumflag=sumflag,nbgr=nbgr,width=width,
                                     pow=pow,tangent=tangent)
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

##############################################################################
def image_bgr(image,lineflag='c',nbgr=3,width=100,pow=2.,
              tangent=False,plot=False):
    """
    Calculate a 2D background for the image.

    * image is the (hopefully clipped) image data

    * lineflag ('c' or 'r') corresponds to to direction which is
      used to generate the background 

    * nbgr = number of end points to use in linear background determination
      (see background.background)
          
    * width should correspond roughly to the actual peak
      width in the lineflag direction. The background should fit
      features that are in general broader than this value
      Note that width = 0 corresponds to no polynomial bgr
      
    * pow is the power of the polynomial used in background determination
      (see background.background)
      
    * tangent is a flag to indicate if local slope of the data should be fitted 
      (see background.background)

     * plot is a flag to indicate if a 'plot' should be made
    """
    bgr_arr = num.zeros(image.shape)

    # fit to rows
    if lineflag=='r':
        for j in range(image.shape[0]):
            bgr_arr[j,:] = background(image[j],nbgr=nbgr,width=width,pow=pow,
                                      tangent=tangent)

    # fit to cols
    if lineflag=='c':
        for j in range(image.shape[1]):
            bgr_arr[:,j] = background(image[:,j],nbgr=nbgr,width=width,pow=pow,
                                      tangent=tangent)

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

    return bgr_arr

################################################################################
class ImageAna:
    def __init__(self,image,roi=[],rotangle=0.0,
                 bgrflag=2,
                 cnbgr=5,cwidth=0,cpow=2.,ctan=False,
                 rnbgr=5,rwidth=0,rpow=2.,rtan=False,
                 plot=True,fig=None,figtitle=''):
        """
        Analyze images
        
        * image is the image data

        * roi is the 'peak' roi = [x1,y1,x2,y2]
          Note take x as the horizontal image axis index (col index), 
          and y as vertical axis index (row index).
          Therefore, in image indexing:
            image[y1:y2,x1:x2] ==> rows y1 to y2 and cols x1 to x2

        * rotangle is the image rotation angle in degrees ccw

        * bgrflag is flag for how to do backgrounds:
           = 0 determine row and column backgrounds after summation
           = 1 determine 2D background using 'c'olumn direction 
           = 2 determine 2D background using 'r'ow direction 
           = 3 determine 2D background using both 'r'ow and 'c'olumn direction 

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
        
        * plot = show fancy plot 

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
        self.rotangle   = rotangle
        self.image      = image
        self.clpimg     = None
        self.bgrimg     = None
        self.integrated = False
        #
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
        self.bgrflag = bgrflag
        self.cbgr = {'nbgr':cnbgr,'width':cwidth,'pow':cpow,'tan':ctan}
        self.rbgr = {'nbgr':rnbgr,'width':rwidth,'pow':rpow,'tan':rtan}
        self.plotflag = plot
        
        self.integrate()

        # plot
        if self.plotflag == True: self.plot(fig=fig)
            
    ############################################################################
    def integrate(self):
        """
        Integrate image.

        Note approx error by assuming data and bgr std deviation are:
           sig_i = sqrt(I_i)
           
        Therefore:
          (sig_Itot)**2 = Sum(I_i)
          (sig_Ibgr)**2 = Sum(Ibgr_i)
          Ierr = sqrt((sig_Itot)**2 + (sig_Ibgr)**2)

        """
        # clip image
        self.clpimg = clip_image(self.image,self.roi,rotangle=self.rotangle)

        # integrate the roi
        self.I = num.sum(num.trapz(self.clpimg))

        # calculate and subtract image background
        self.bgrimg = None
        self.Ibgr = 0.0
        if self.bgrflag > 0:
            if self.bgrflag == 1:
                self.bgrimg = image_bgr(self.clpimg,lineflag='c',nbgr=self.cbgr['nbgr'],
                                        width=self.cbgr['width'],pow=self.cbgr['pow'],
                                        tangent=self.cbgr['tan'],plot=False)
            elif self.bgrflag == 2:
                self.bgrimg = image_bgr(self.clpimg,lineflag='r',nbgr=self.rbgr['nbgr'],
                                        width=self.rbgr['width'],pow=self.rbgr['pow'],
                                        tangent=self.rbgr['tan'],plot=False)
            else:
                bgr_r = image_bgr(self.clpimg,lineflag='c',nbgr=self.cbgr['nbgr'],
                                  width=self.cbgr['width'],pow=self.cbgr['pow'],
                                  tangent=self.cbgr['tan'],plot=False)
                bgr_c = image_bgr(self.clpimg,lineflag='r',nbgr=self.rbgr['nbgr'],
                                  width=self.rbgr['width'],pow=self.rbgr['pow'],
                                  tangent=self.rbgr['tan'],plot=False)
                # combine the two bgrs by taking avg
                self.bgrimg = (bgr_r + bgr_c)/2.
                
            # correct for 2D bgr
            self.Ibgr = num.sum(num.trapz(self.bgrimg))
            self.I    = self.I - self.Ibgr
        
        # error
        self.Ierr = (self.I + self.Ibgr)**0.5
        
        # integrate col sum
        if self.bgrimg != None:
            (I,Ierr,Ibgr) = line_sum_integral(self.clpimg-self.bgrimg,sumflag='c',nbgr=0)
        else:
            (I,Ierr,Ibgr) = line_sum_integral(self.clpimg,sumflag='c',
                                              nbgr=self.rbgr['nbgr'],
                                              width=self.rbgr['width'],
                                              pow=self.rbgr['pow'],
                                              tangent=self.rbgr['tan'])
        self.I_c      = I
        self.Ierr_c   = Ierr
        self.Ibgr_c   = Ibgr
        
        # integrate row sum
        if self.bgrimg != None:
            (I,Ierr,Ibgr) = line_sum_integral(self.clpimg-self.bgrimg,sumflag='r',nbgr=0)
        else:
            (I,Ierr,Ibgr) = line_sum_integral(self.clpimg,sumflag='r',
                                              nbgr=self.cbgr['nbgr'],
                                              width=self.cbgr['width'],
                                              pow=self.cbgr['pow'],
                                              tangent=self.cbgr['tan'])
        self.I_r      = I
        self.Ierr_r   = Ierr
        self.Ibgr_r   = Ibgr

    ############################################################################
    def plot(self,fig=None):
        """
        make fancy 4-panel plot
        """
        if self.integrated == False:
            self.integrate()
        
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
        if self.bgrimg != None:
            title_roi = title_roi + '\n(background subtracted)'

        # calc full image with an roi box
        (c1,r1,c2,r2) = self.roi
        if self.rotangle != 0.0:
            bild = ndimage.rotate(self.image,self.rotangle)
        else:
            bild = copy.copy(self.image)
        bild[r1-1:r1,c1:c2] = bild.max()
        bild[r2:r2+1,c1:c2] = bild.max()
        bild[r1:r2,c1-1:c1] = bild.max()
        bild[r1:r2,c2:c2+1] = bild.max()

        ###################################
        #### plot column sum
        pylab.subplot(221)
        pylab.title(title_c, fontsize = 12)
        # plot raw sum
        (data, data_idx, bgr) = line_sum(self.clpimg,sumflag='c',nbgr=0)
        rawmax = data.max()
        pylab.plot(data_idx, data, 'k',label='raw sum')
        # get bgr and data-bgr
        if self.bgrimg != None:
            # here data is automatically bgr subracted
            (data, data_idx, xx) = line_sum(self.clpimg-self.bgrimg,sumflag='c',nbgr=0)
            bgr = self.bgrimg.sum(axis=0)
        else:
            # here data is data and bgr is correct, therefore data = data-bgr
            (data, data_idx, bgr) = line_sum(self.clpimg,sumflag='c',
                                             nbgr=self.rbgr['nbgr'],
                                             width=self.rbgr['width'],
                                             pow=self.rbgr['pow'],
                                             tangent=self.rbgr['tan'])
            data = data-bgr
        # plot bgr and bgr subtracted data
        pylab.plot(data_idx, bgr, 'r',label='bgr')
        pylab.plot(data_idx, data, 'b',label='data-bgr')
        pylab.axis([0, data_idx.max(), 0, rawmax*1.25])
        pylab.legend()

        ####################################
        # plot full image with ROI
        pylab.subplot(222)
        pylab.title(self.title, fontsize = 12)
        pylab.imshow(bild,cmap=colormap)
        pylab.colorbar(orientation='horizontal')

        ####################################
        # plot zoom on image 
        pylab.subplot(223)
        pylab.title(title_roi,fontsize = 12)
        if self.bgrimg != None:
            pylab.imshow(self.clpimg-self.bgrimg, cmap=colormap, aspect='auto')
        else:
            pylab.imshow(self.clpimg, cmap=colormap, aspect='auto')
            
        ####################################
        # plot row sum
        pylab.subplot(224)
        pylab.title(title_r, fontsize = 12)
        # plot raw sum
        (data, data_idx, bgr) = line_sum(self.clpimg,sumflag='r',nbgr=0)
        rawmax = data.max()
        pylab.plot(data, data_idx, 'k',label='raw sum')
        # get bgr and data-bgr
        if self.bgrimg != None:
            # here data is automatically bgr subracted
            (data, data_idx, xx) = line_sum(self.clpimg-self.bgrimg,sumflag='r',nbgr=0)
            bgr = self.bgrimg.sum(axis=1)
        else:
            # here data is data and bgr is correct, therefore data = data-bgr
            (data, data_idx, bgr) = line_sum(self.clpimg,sumflag='r',
                                             nbgr=self.cbgr['nbgr'],
                                             width=self.cbgr['width'],
                                             pow=self.cbgr['pow'],
                                             tangent=self.cbgr['tan'])
            data = data-bgr
        # plot bgr and bgr subtracted data
        pylab.plot(bgr, data_idx, 'r',label='bgr')
        pylab.plot(data, data_idx, 'b',label='data-bgr')
        pylab.axis([0,rawmax*1.25, data_idx.max(), 0])
        pylab.xticks(rotation=-45)
        pylab.legend()

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
    
 
    
    