#######################################################################
"""
Tom Trainor (fftpt@uaf.edu)
Class to store and operate on generic 'scan data'

Modifications:
--------------

"""
#######################################################################
"""
Todo

"""
########################################################################

import types
from   glob import glob
import os, copy
import numpy as Num

from   specfile import SpecFile
import image_data
import xrf_ops

########################################################################
class ScanData:
    """
    ScanData Class holds:
    name = scan name
    dims ->[npts1] for a 1-d scan
           [npts1,npts2] for a 2-d scan
            where npts1 is the outer loop
            and npts2 is for the inner loop
    scalers = {'io':[],'i1':[]}  -> each array has dims = dims
    positioners  = {'E':[]}      -> these may be single values,
                                    or same size as scalers
    primary_axis = ['E']         -> if multi-dim, list of outer, inner
    primary_det  = 'i1'          -> use for default plotting
    med = []                     -> array has dims = dims or is [] for No med
    xrf = []                     -> array has dims = dims or is [] for No xrf
    images = []                  -> array has dims = dims or is [] for No images
    state = {}                   -> dictionary of additional state information
    """
    ################################################################
    def __init__(self,name='',dims=[],scalers={},positioners={},
                 primary_axis=[],primary_det=None,state={},
                 med=[],xrf=[],xrf_lines=None,image=[],image_rois=[]):
        
        self.name         = name
        self.dims         = dims
        self.primary_axis = primary_axis
        self.primary_det  = primary_det
        self.scalers      = scalers
        self.positioners  = positioners
        self.state        = state
        #
        self.med          = med
        self.xrf          = xrf
        self.xrf_lines    = xrf_lines
        self.xrf_peaks    = {}
        self.image        = image
        self.image_rois   = image_rois
        self.image_peaks  = {}
        #
        self.bad_points   = []

    ########################################################################
    """
    def __setitem__(self,arg,val):
        # note this currenlty craps for sub arrays.
        # ie x[0][0].  In such a case arg is a tuple
        # first is the name, the rest are the subscripts
        #print arg, type(arg)
        #print val, type(val)
        if type(arg) == types.StringType:
            if arg.startswith("_"):
                raise TypeError, "Private attribute"
        #if hasattr(self,arg):
        #    setattr(self,arg,val)
        setattr(self,arg,val)
        return
    """
    ########################################################################
    def __getitem__(self,arg):
        """
        Get items.  Note that order matters:
        spectra, images, scalers, positioners, state
        """
        if type(arg) == types.StringType:
            #if arg.startswith("_"):
            #    raise TypeError, "Private attribute"
            #return getattr(self,arg)
            if arg == 'med':
                return self.med
            if arg == 'xrf':
                return self.xrf
            if arg == 'image':
                return self.image
            #
            for key in self.xrf_peaks.keys():
                if key == arg:
                    return self.xrf_peaks[key]
            #
            for key in self.image_peaks.keys():
                if key == arg:
                    return self.image_peaks[key]
            #
            for key in self.scalers.keys():
                if key == arg:
                    return self.scalers[key]
            for key in self.positioners.keys():
                if key == arg:
                    return self.positioners[key]
            for key in self.state.keys():
                if key == arg:
                    return self.state[key]
        return None

    ########################################################################
    def __repr__(self):
        lout = "* Scan Name = %s\n" % self.name
        lout = lout + "* Scan Dimensions = %s\n" % str(self.dims)

        lout = lout + "* Scalers:\n"
        ct = 0
        for sc in self.scalers.keys():
            lout = lout + "%12s" % sc
            ct = ct + 1
            if ct > 7:
                lout = lout + '\n'
                ct = 0

        lout = lout + "\n* Positioners:\n"
        ct = 0
        for p in self.positioners.keys():
            lout = lout + "%12s" % p
            ct = ct + 1
            if ct > 7:
                lout = lout + '\n'
                ct = 0

        if len(self.state) > 0:
            lout = lout + "\n* Additional State Variables:\n"  
            ct = 0
            for p in self.state.keys():
                lout = lout + "%12s" % p
                ct = ct + 1
                if ct > 7:
                    lout = lout + '\n'
                    ct = 0

        lout = lout + "\n* Primary scan axis = %s\n"  % str(self.primary_axis)
        lout = lout + "* Primary detector  = %s\n"  % self.primary_det

        spectra = False
        if self.med != []:
            spectra = True
            lout = lout + "* Scan includes %i med spectra\n" % len(self.med)
        if self.xrf != []:
            spectra = True
            lout = lout + "* Scan includes %i xrf spectra\n" % len(self.xrf)
            if self.xrf_lines != None:
                lout = lout + "  -> Xrf lines = %s\n" % str(self.xrf_lines)
        if spectra == False:
            lout = lout + "* Scan does not include spectra\n"

        if self.image != []:
            lout = lout + "* Scan includes %i images\n" % len(self.image)
        else:
            lout = lout + "* Scan does not include images\n"
        
        return lout
    
    ################################################################
    def get_scaler(self,label=None):
        """
        return scaler
        """
        if label == None:
            label = self.primary_det[0]
        return self.scalers.get(label)
        
    ################################################################
    def get_positioner(self,label=None):
        """
        return positioner
        """
        if label == None:
            label = self.primary_axis[0]
        return self.positioners.get(label)

    ################################################################
    def med_ocr(self):
        """
        return outgoing count rate (ocr) values from meds
        """
        ndet = self.med[0].n_detectors
        npnt = len(self.med)
        ocr  = Num.zeros((npnt,ndet))
        for j in range(npnt):
            for k in range(ndet):
                tot = self.med[j].mca[k].total_counts
                lt  = self.med[j].mca[k].live_time
                ocr[j][k] = float(tot)/float(lt)
        return Num.transpose(ocr)

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
        xrf = xrf_ops.med2xrf(self.med,xrf_params=xrf_params,
                              lines = self.xrf_lines,
                              det_idx=det_idx,emin=emin,emax=emax)
        if xrf: self.xrf = xrf

    ################################################################
    ################################################################
    def get_xrf(self,pnt=0):
        """
        get xrf
        """
        pnt = int(pnt)
        if self.xrf == []:
            return None
        if pnt not in range(self.dims[0]):
            return None
        return self.xrf[pnt]

    ################################################################
    def init_xrf_lines(self,lines=None):
        """
        init xrf lines
        """
        if lines == None:
            lines = self.xrf_lines
        if type(lines) != types.ListType:
            lines = [lines]
        for x in self.xrf:
            x.init_lines(lines)

    ################################################################
    def xrf_calc(self):
        """
        calc xrf
        """
        for x in self.xrf:
            x.calc()
        self.update_xrf_peaks()

    ################################################################
    def xrf_fit(self,xrf_params={},use_prev_fit=False,fit_init=-1,
                guess=False,verbose=True):
        """
        fit xrf
        """
        xrf_ops.fit(self.xrf,xrf_params=xrf_params,use_prev_fit=use_prev_fit,
                    fit_init=fit_init,guess=guess,verbose=verbose)
        self.update_xrf_peaks()

    ################################################################
    def update_xrf_peaks(self,):
        """
        update xrf peaks
        """
        lines = []
        for pk in self.xrf[0].peaks:
            lines.append(pk.label)
        self.xrf_lines = lines
        for l in lines:
            p = xrf_ops.peak_areas(self.xrf,l)
            self.xrf_peaks[l] = p

    ################################################################
    ################################################################
    def integrate_image(self,idx=[],roi=[],bgr_params={},plot=True,fig=None):
        """
        integrate images
        roi  = [x1,y1,x2,y2]
        """
        # make sure arrays exist:
        init = False
        if len(self.image_peaks) > 0:
            if len(self.image_peaks['I_c']) != len(self.image):
                init = True
        else:
            init = True
        if init:
            self._init_image()

        # idx of images to integrate
        if len(idx) == 0:
            idx = Num.arange(len(self.image))
        
        # update roi
        if len(roi) == 4:
            for j in idx:
                self.image_rois[j] = roi
        elif len(roi) == len(idx):
            for j in idx:
                self.image_rois[j] = roi[j]
        
        # do integrations
        for j in idx:
            if self.image_peaks['I_c'][j] != -1:
                self._integrate_image(idx=j,roi=self.image_rois[j],
                                      bgr_params=bgr_params,plot=plot,fig=fig)
    
    ################################################################
    def _init_image(self):
        npts = len(self.image)
        if npts == 0:
            self.image_rois = []
            self.image_peaks = {}

        if self.image_rois == None:
            self.image_rois = []
            
        if len(self.image_rois) != npts:
            self.image_rois = []
            for j in range(npts):
                self.image_rois.append([])

        # should we init all these or set based on an integrate flag?
        self.image_peaks  = {}
        self.image_peaks['I']      = Num.zeros(npts,dtype=float)
        self.image_peaks['Ierr']   = Num.zeros(npts,dtype=float)
        self.image_peaks['Ibgr']   = Num.zeros(npts,dtype=float)
        #
        self.image_peaks['I_c']    = Num.zeros(npts,dtype=float)
        self.image_peaks['Ierr_c'] = Num.zeros(npts,dtype=float)
        self.image_peaks['Ibgr_c'] = Num.zeros(npts,dtype=float)
        #
        self.image_peaks['I_r']    = Num.zeros(npts,dtype=float)
        self.image_peaks['Ierr_r'] = Num.zeros(npts,dtype=float)
        self.image_peaks['Ibgr_r'] = Num.zeros(npts,dtype=float)
        
    ################################################################
    def _integrate_image(self,idx=0,roi=[],bgr_params={},plot=True,fig=None):
        """
        integrate an image
        """
        if idx < 0 or idx > len(self.image): return None
        #
        figtitle = "Scan Point = %i, L = %6.3f" % (idx,self.scalers['L'][idx])
        nbgr   = bgr_params.get('nbgr')
        cwidth = bgr_params.get('cwidth')
        rwidth = bgr_params.get('rwidth')
        if  nbgr   == None: nbgr=0
        if  cwidth == None: cwidth=0
        if  rwidth == None: rwidth=0

        img_ana = image_data.ImageAna(self.image[idx],roi=roi,
                                      nbgr=nbgr,cwidth=cwidth,rwidth=rwidth,
                                      plot=plot,fig=fig,figtitle=figtitle)

        # results into image_peaks dictionary
        self.image_peaks['I'][idx]      = img_ana.I
        self.image_peaks['Ierr'][idx]   = img_ana.Ierr
        self.image_peaks['Ibgr'][idx]   = img_ana.Ibgr
        #
        self.image_peaks['I_c'][idx]    = img_ana.I_c
        self.image_peaks['Ierr_c'][idx] = img_ana.Ierr_c
        self.image_peaks['Ibgr_c'][idx] = img_ana.Ibgr_c
        #
        self.image_peaks['I_r'][idx]    = img_ana.I_r
        self.image_peaks['Ierr_r'][idx] = img_ana.Ierr_r
        self.image_peaks['Ibgr_r'][idx] = img_ana.Ibgr_r
        
##############################################################################

