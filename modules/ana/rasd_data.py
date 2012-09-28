"""
Function for extracting RASD data from spec "raxr" scans
mainly copied from ctr_data

Authors/Modifications:
----------------------
* Frank Heberling (Frank.Heberling@kit.edu)
* Tom Trainor

Notes:
------
To use the function you may type in these 4 variables with the values you need:

    spec_path = 'D:/Data/Surf_Diff/Aps_Data_ID' (path to your spec file)
    spec = 'caco3_july09h.spc'                  (name of your spec file)
    first_scan = 45 
    last_scan = 143 (first and last hklscan in your spec file, that belong to the raxr scan that is to be integrated)
    
and run it like this:

    A = ana.rasd_data.rasd_data(spec_path, spec, first_scan, last_scan)

"rasd_data" appends all the hklscans specified to one RasdData object 'A' (the two images and the 'io's in one
hklscan are summed up).
'A' can then be integrated in the CTRData application.
"""
#############################################################################

import types
import pylab
import numpy as num

from   tdl.modules.specfile.reader import Reader
from   tdl.modules.ana.scan_data import ScanData
from   tdl.modules.ana import image_data
from   tdl.modules.ana import ctr_data
from   tdl.modules.ana.ctr_data import *

#############################################################################
def rasd_data(spec_path,spec,first_scan,last_scan):
    """
    Merges all the hklscans of a spec "raxr" scan to one RASDData object
    containing images and info from the spec file, integrates the images
    interactively (as a ScanData object), calculates polarization + lorentz correction
    factors, plots the RASD wiggle, and writes the Information into a .rsd file.
    """
    spec_file = Reader(spec, spec_path)
    spec_file.image = False
    nscans = last_scan - first_scan + 1

    S = []
    for i in range(nscans):
        tmp = spec_file.spec_scan(first_scan + i, image = True)
        tmp.dims[0] = 1
        tmp.scalers['ROI3']= num.array([tmp.scalers['ROI3'][0]])
        tmp.scalers['i1']= num.array([tmp.scalers['i1'][0]])
        tmp.scalers['H']= num.array([tmp.scalers['H'][0]])
        tmp.scalers['K']= num.array([tmp.scalers['K'][0]])
        tmp.scalers['ROI2']= num.array([tmp.scalers['ROI2'][0]])
        tmp.scalers['Beta']= num.array([tmp.scalers['Beta'][0]])
        tmp.scalers['IROI']= num.array([tmp.scalers['IROI'][0]])
        tmp.scalers['Epoch']= num.array([tmp.scalers['Epoch'][0]])
        tmp.scalers['ROI1']= num.array([tmp.scalers['ROI1'][0]])
        tmp.scalers['io']= num.array([tmp.scalers['io'][0]+tmp.scalers['io'][1]])
        tmp.scalers['Alpha']= num.array([tmp.scalers['Alpha'][0]])
        tmp.scalers['AmpTek_sc']= num.array([tmp.scalers['AmpTek_sc'][0]])
        tmp.scalers['Seconds']= num.array([tmp.scalers['Seconds'][0]])
        tmp.scalers['Bicron']= num.array([tmp.scalers['Bicron'][0]])
        tmp.positioners['chi']= num.array([tmp.positioners['chi'][0]])
        tmp.positioners['eta']= num.array([tmp.positioners['eta'][0]])
        tmp.positioners['phi']= num.array([tmp.positioners['phi'][0]])
        tmp.positioners['nu']= num.array([tmp.positioners['nu'][0]])
        tmp.positioners['mu']= num.array([tmp.positioners['mu'][0]])
        tmp.positioners['del']= num.array([tmp.positioners['del'][0]])
        tmp.image.image = [tmp.image.image[0] + tmp.image.image[1]]
        tmp.image.bgrpar = [tmp.image.bgrpar[0]]
        tmp.image.im_max = [tmp.image.im_max[0]]
        tmp.image.rois = [tmp.image.rois[0]]
        tmp.image.rotangle = [tmp.image.rotangle[0]]
        for i in tmp.image.peaks.keys():
            tmp.image.peaks[i] = num.array([tmp.image.peaks[i][0]])
        S.append(tmp)

    Z = RasdData(S)


    return Z


##################################################################################
##################################################################################
class RasdData:
    """
    RASD data

    Attributes:
    -----------
    * bad is a list of index values of points flagged as 'bad'
    * scan is a list holding all the scan data ojects
    * scan_index is a list (possibly of tuples) that give the
      scan index corresponding to a data point. see the get_scan method
    * labels is a dictionary of lists of the string labels
      for 'I','Inorm','Ierr','Ibgr'
    * corr_params = list of correction parameter dictionaries
    * scan_type list of string flags for scan type
    * H     = array of H values
    * K     = array of K values
    * L     = array of L values
    * E     = array of Energy values
    * I     = array of integrated intensities
    * Inorm = array of normalization values
    * Ierr  = array of intensity errors
    * Ibgr  = array of intensity backgrounds
    * ctot  = array of total correction factors
    * F     = array of structure factor magnitudes
    * Ferr  = array of structure factor error bars
    """
    ##########################################################################
    def __init__(self,scans=[],I='I',Inorm='io',Ierr='Ierr',
                 Ibgr='Ibgr',corr_params={},scan_type='image'):
        """
        Initialize the object.

        Parameters:
        -----------
        * scans = list of scan data instances

        * I = string label corresponding to the intensity array, ie
          let y be the intensity, then y = scan[I]

        * Inorm=string label corresponding to the normalization intensity array, 
          ie, norm = scan[Inorm], where the normalized intensity is taken as
          ynorm= y/norm

        * Ierr = string label corresponding to the intensity error array, ie
          y_err = scan[Ierr].  We assume these are standard deviations 
          of the intensities (y).
        
          Note when the data are normalized we assume the statistics of norm 
          go as norm_err = sqrt(norm).  Therefore the error of normalized
          data is
            ynorm_err = ynorm * sqrt( (y_err/y)^2 + (norm_err/(norm))^2 )
                      = ynorm * sqrt( (y_err/y)^2 + 1/norm )
          If its assumed that y_err = sqrt(y) then the expression could be
          simplified futher, but we wont make that assumption since there
          may be other factors that determine how yerr was calculated.  

        * Ibgr = string label corresponding to background intensity array

        * corr_params = a dictionary containing the necessary information for
          data corrections.
          corr_params['geom'] = Type of goniometer geometry ('psic' is default)
          corr_params['beam_slits'] = dictionary describing the beam slits,e.g.
                                      {'horz':.6,'vert':.8}
          corr_params['det_slits'] = dictionary describing the beam slits,e.g.
                                     {'horz':.6,'vert':.8}
          corr_params['sample'] = a dictionary describing the sample shape.
                                  {'dia':0,'angles':{},'polygon':[]}
                                  if dia >=0 then assume a round sample.
                                  Otherwise use polygon/angles
                                  If all are None, then no sample correction is
                                  computed
          corr_params['scale'] = scale factor to multiply by all the intensity
                                 values. e.g.  if Io ~ 1million cps
                                 then using 1e6 as the scale makes the normalized
                                 intensity close to cps.  ie y = scale*y/norm

          See the Correction class for more info...

        * scan_type = Type of scans (e.g. 'image', 'phi', etc..)
        
        """
        self.fig    = None
        self.cursor = None
        self.bad    = []
        self.scan   = []
        self.scan_index  = []
        #
        self.labels      = {'I':[],'Inorm':[],'Ierr':[],'Ibgr':[]}
        self.corr_params = []
        self.scan_type   = []
        #
        self.H     = num.array([],dtype=float)
        self.K     = num.array([],dtype=float)
        self.L     = num.array([],dtype=float)
        self.E     = num.array([],dtype=float)
        self.I     = num.array([],dtype=float)
        self.Inorm = num.array([],dtype=float)
        self.Ierr  = num.array([],dtype=float)
        self.Ibgr  = num.array([],dtype=float)
        self.ctot  = num.array([],dtype=float)
        self.F     = num.array([],dtype=float)
        self.Ferr  = num.array([],dtype=float)
        #
        self.append_scans(scans,I=I,Inorm=Inorm,Ierr=Ierr,Ibgr=Ibgr,
                          corr_params=corr_params,
                          scan_type=scan_type)

    ##########################################################################
    def __repr__(self,):
        """ display """
        lout = "RASD DATA\n"
        lout = "%sNumber of scans = %i\n" % (lout,len(self.scan))
        lout = "%sNumber of structure factors = %i\n" % (lout,len(self.L))
        return lout

    ##########################################################################
    def __save__(self,):
        """
        this is called by the pickler so we can delete stuff
        we dont want to pickle
        """
        del self.cursor
        self.cursor = None

    ##########################################################################
    def append_scans(self,scans,I=None,Inorm=None,Ierr=None,Ibgr=None,
                     corr_params=None,scan_type=None):
        """
        Append new scan data

        Parameters:
        -----------
        * scans is a list of scan data objects.

        * The rest of the arguments (defined in __init__)
          should be the same for each scan in the list.

        For any argument with None passed we use previous defined
        values - based on the last exisiting data point.  
        """
        if type(scans) != types.ListType:
            scans = [scans]
        
        # if None passed use the last values
        if I == None:           I = self.labels['I'][-1]
        if Inorm == None:       Inorm = self.labels['Inorm'][-1]
        if Ierr == None:        Ierr = self.labels['Ierr'][-1]
        if Ibgr == None:        Ibgr = self.labels['Ibgr'][-1]
        if corr_params == None: corr_params = self.corr_params[-1]
        if scan_type == None:   scan_type = self.scan_type[-1]

        # get all the data parsed out of each scan and append
        for scan in scans:
            data = self._scan_data(scan,I,Inorm,Ierr,Ibgr,corr_params,scan_type)
            if data == None: return

            #self.scan.append([])
            self.scan.append(scan)
            #
            self.scan_index.extend(data['scan_index'])
            self.labels['I'].extend(data['I_lbl'])
            self.labels['Inorm'].extend(data['Inorm_lbl'])
            self.labels['Ierr'].extend(data['Ierr_lbl'])
            self.labels['Ibgr'].extend(data['Ibgr_lbl'])
            self.corr_params.extend(data['corr_params'])
            self.scan_type.extend(data['scan_type'])
            #
            self.H     = num.append(self.H,data['H'])
            self.K     = num.append(self.K,data['K'])
            self.L     = num.append(self.L,data['L'])
            self.I     = num.append(self.I,data['I'])
            self.Inorm = num.append(self.Inorm,data['Inorm'])
            self.Ierr  = num.append(self.Ierr,data['Ierr'])
            self.Ibgr  = num.append(self.Ibgr,data['Ibgr'])
            self.ctot  = num.append(self.ctot,data['ctot'])
            self.F     = num.append(self.F,data['F'])
            self.Ferr  = num.append(self.Ferr,data['Ferr'])
            self.E     = num.append(self.E,data['E'])

    ##########################################################################
    def _scan_data(self,scan,I,Inorm,Ierr,Ibgr,corr_params,scan_type):
        """
        Parse scan into data...
        """
        data = {'scan_index':[],'I_lbl':[],'Inorm_lbl':[],
                'Ierr_lbl':[],'Ibgr_lbl':[],'corr_params':[],'scan_type':[],
                'H':[],'K':[],'L':[],'I':[],'Inorm':[],'Ierr':[],'Ibgr':[],
                'ctot':[],'F':[],'Ferr':[],'E':[]}

        # compute a scan index
        scan_idx = len(self.scan)
        
        # image scan -> each scan point is a unique HKL
        if scan_type == 'image':
            if scan.image._is_integrated == False:
                scan.image.integrate()
            npts = int(scan.dims[0])
            for j in range(npts):
                data['scan_index'].append((scan_idx,j))
                data['I_lbl'].append(I)
                data['Inorm_lbl'].append(Inorm)
                data['Ierr_lbl'].append(Ierr)
                data['Ibgr_lbl'].append(Ibgr)
                data['corr_params'].append(corr_params)
                data['scan_type'].append(scan_type)
                #
                data['H'].append(scan['H'][j])
                data['K'].append(scan['K'][j])
                data['L'].append(scan['L'][j])
                data['E'].append(scan.state['ENERGY'][0])
                # get F
                d = image_point_F(scan,j,I=I,Inorm=Inorm,
                                  Ierr=Ierr,Ibgr=Ibgr,
                                  corr_params=corr_params)
                data['I'].append(d['I'])
                data['Inorm'].append(d['Inorm'])
                data['Ierr'].append(d['Ierr'])
                data['Ibgr'].append(d['Ibgr'])
                data['ctot'].append(d['ctot'])
                data['F'].append(d['F'])
                data['Ferr'].append(d['Ferr'])
        return data

    ##########################################################################
    def integrate_point(self,idx,**kw):
        """
        (Re)-integrate an individual structure factor point.

        Parameters:
        -----------
        * idx is the index number of the point

        If scan type is image use the following kw arguments
        (if the argument is not passed the existing value is used,
        ie just use these to update/change parameters)
        
        * bad        = True/False flags point as bad or not
        * roi        = image roi
        * rotangle   = image rotation angle
        * bgr_params = image background parameters
        * plot       = True/False to show integration plot 
        * fig        = Fig number for plot
        * I          = Intensity label
        * Inorm      = Intensity norm label
        * Ierr       = Intensity error label
        * Ibgr       = Intensity background label
        * corr_params = CTR correction parameters
        
        """
        if idx not in range(len(self.L)): return

        bad = kw.get('bad')
        if bad != None:
            if bad == True:
                if idx not in self.bad:
                    self.bad.append(idx)
            elif bad == False:
                if idx in self.bad:
                    self.bad.remove(idx)
            else:
                print "Warning: bad should be True/False"

        if self.scan_type[idx]=="image":
            (scan_idx,point) = self.scan_index[idx]
            scan = self.scan[scan_idx]
            if scan.image._is_init() == False:
                scan.image._init_image()
            
            # parse integration parameters
            roi        = kw.get('roi')
            rotangle   = kw.get('rotangle')
            bgr_params = kw.get('bgr_params')
            #
            plot       = kw.get('plot',False)
            fig        = kw.get('fig')
            
            if idx in self.bad:
                bad = [point]
            else:
                bad = []
            # integrate the scan.  
            scan.image.integrate(idx=[point],roi=roi,rotangle=rotangle,
                                 bgr_params=bgr_params,bad_points=bad,
                                 plot=plot,fig=fig)
            
            # parse all the correction info and re-compute 
            I           = kw.get('I',self.labels['I'][idx])
            Inorm       = kw.get('Inorm',self.labels['Inorm'][idx])
            Ierr        = kw.get('Ierr', self.labels['Ierr'][idx])
            Ibgr        = kw.get('Ibgr', self.labels['Ibgr'][idx])
            corr_params = kw.get('corr_params',self.corr_params[idx])
            d = image_point_F(scan,point,I=I,Inorm=Inorm,
                              Ierr=Ierr,Ibgr=Ibgr,
                              corr_params=corr_params)
            
            # store results
            self.labels['I'][idx]     = I
            self.labels['Inorm'][idx] = Inorm
            self.labels['Ierr'][idx]  = Ierr
            self.corr_params[idx]     = corr_params
            self.H[idx]               = scan['H'][point]
            self.K[idx]               = scan['K'][point]
            self.L[idx]               = scan['L'][point]
            self.I[idx]               = d['I']
            self.Inorm[idx]           = d['Inorm']
            self.Ierr[idx]            = d['Ierr']
            self.Ibgr[idx]            = d['Ibgr']
            self.ctot[idx]            = d['ctot']
            self.F[idx]               = d['F']
            self.Ferr[idx]            = d['Ferr']
        return 



    ##########################################################################
    def plot(self,fig=None,num_col=1,cursor=True,verbose=True,spnt=None):
        """
        Plot the structure factor data vs L

        Parameters:
        -----------
        * fig is the figure number to plot to
        * num_col is the number of plot columns
        * cursor is a flag to indicate cursor updates
        * verbose is a flag to indicate if cursor clicks
          should also print info
        * spnt is point index to plot in red
        """
        hksets  = sort_rasd(self)
        nset    = len(hksets)
        num_col = float(num_col)
        num_row = num.ceil(nset/num_col)
        pyplot.figure(fig)
        pyplot.clf()
        for j in range(nset):
            pyplot.subplot(num_row,num_col,j+1)
            d = hksets[j]
            title = 'H=%2.2f,K=%2.2f,L=%2.3f' % (d['H'][0],d['K'][0],d['L'][0])
            pyplot.title(title, fontsize = 12)
            pyplot.plot(d['E'],d['F'],'b.-')
            pyplot.errorbar(d['E'],d['F'],d['Ferr'], fmt ='o')
            if spnt != None:
                if spnt in d['point_idx']:
                    idx = num.where(d['point_idx']==spnt)
                    pyplot.plot(d['E'][idx],d['F'][idx],'ro')
            #
            min_L = num.floor(num.min(d['E']))
            max_L = num.ceil(num.max(d['E']))
            idx   = num.where(d['F'] > 0.)
            min_F = num.min(d['F'][idx])
            #min_F = 10.**(num.round(num.log10(min_F)) - 1)
            max_F = num.max(d['F'][idx])
            #max_F = 10.**(num.round(num.log10(max_F)) + 1)
            #pyplot.axis([min_L,max_L,min_F,max_F])
            #
            pyplot.xlabel('E')
            pyplot.ylabel('|F|')    
        fig = pyplot.gcf()
        self.fig = fig.number
        if self.cursor != None:
            self.cursor._disconnect()
            del self.cursor
            self.cursor = None
        if cursor == True:
            self.cursor = plotter.cursor(fig=self.fig,verbose=verbose)
    ##########################################################################
    def plot_I(self,fig=None,num_col=1,cursor=True,verbose=True,spnt=None):
        """
        Plot the raw intensities vs E

        Parameters:
        -----------
        * fig is the figure number to plot to
        * num_col is the number of plot columns
        * cursor is a flag to indicate cursor updates
        * verbose is a flag to indicate if cursor clicks
          should also print info
        * spnt is point index to plot in red
        
        """
        hksets  = sort_rasd(self)
        nset    = len(hksets)
        num_col = float(num_col)
        num_row = num.ceil(nset/num_col)
        pyplot.figure(fig)
        pyplot.clf()
        for j in range(nset):
            pyplot.subplot(num_row,num_col,j+1)
            d = hksets[j]
            title = 'H=%2.2f,K=%2.2f,L=%2.3f' % (d['H'][0],d['K'][0],d['L'][0])
            pyplot.title(title, fontsize = 12)
            I  = d['I']
            In = d['Inorm']
            Ib = d['Ibgr']
            Ie = d['Ierr']
            y  = I/In
            yb = Ib/In
            ye = Ie/In
            #
            pyplot.plot(d['E'],y,'b.-',label='I/Inorm')
            pyplot.errorbar(d['E'],y,ye, fmt ='o')
            pyplot.plot(d['E'],yb,'m.-',label='Ibgr/Inorm')
            #
            min_L = num.floor(num.min(d['E']))
            max_L = num.ceil(num.max(d['E']))
            #
            idx   = num.where(y > 0.)
            min_I = num.min(y[idx])
            min_I = 10.**(num.round(num.log10(min_I)) - 1)
            idxb  = num.where(yb > 0.)
            if len(idxb[0]) == 0:
                min_Ibgr = min_I
            else:
                min_Ibgr = num.min(yb[idxb])
                min_Ibgr = 10.**(num.round(num.log10(min_Ibgr)) - 1)
            min_I = min(min_I,min_Ibgr)
            #
            max_I = num.max(y[idx])
            max_I = 10.**(num.round(num.log10(max_I)) + 1)
            if len(idxb[0]) == 0:
                max_Ibgr = max_I
            else:
                max_Ibgr = num.max(yb[idxb])
                max_Ibgr = 10.**(num.round(num.log10(max_Ibgr)) + 1)
            max_I = max(max_I,max_Ibgr)
            #
            idx  = num.where(I <= 0.)
            tmp  = I[idx] * 0.0 + 1.1*min_I 
            pyplot.plot(d['E'][idx],tmp,'bo')
            #
            if spnt != None:
                if spnt in d['point_idx']:
                    idx = num.where(d['point_idx']==spnt)
                    if I[idx] <= 0.:
                        pyplot.plot(d['E'][idx],[1.1*min_I],'ro')
                    else:
                        pyplot.plot(d['E'][idx],y[idx],'ro')
            #
            pyplot.axis([min_L,max_L,min_I,max_I])
            pyplot.xlabel('E')
            pyplot.ylabel('Intensity')
            pyplot.legend()
        fig = pyplot.gcf()
        self.fig = fig.number
        if self.cursor != None:
            self.cursor._disconnect()
            del self.cursor
            self.cursor = None
        if cursor == True:
            self.cursor = plotter.cursor(fig=self.fig,verbose=verbose)

    ##########################################################################
    def plot_point(self,idx=None,fig=None,show_int=False,cmap=None):
        """
        Plot the raw data for a selected point

        Parameters:
        -----------
        * idx = point index.
          if idx = None, then uses last cursor click
        * fig = fig number
        * show_int is a flag for showing integration plot for images
        * cmap is the color map for images
        """
        if idx == None:
            idx = self.get_idx()
        if self.scan_type[idx] == 'image':
            if show_int:
                self.integrate_point(idx,plot=True,fig=fig)
            else:
                (scan_idx,point) = self.scan_index[idx]
                self.scan[scan_idx].image.plot(idx=point,fig=fig,cmap=cmap)
        else:
            # plot scan data
            pass

    ##########################################################################
    def get_idx(self):
        """
        Get point index from last plot selection

        Example:
        -------
        >>idx = ctr.get_idx()
        >>L = ctr.L[idx]
        """
        if self.cursor == None:
            return None
        if self.cursor.clicked == False:
            #self.cursor.get_click(msg="Select a data point")
            return None
        L = self.cursor.x
        subplot = self.cursor.subplot
        if subplot < 0:
            return None
        hksets  = sort_data(self)
        return self._get_idx(subplot,L,hksets)
    
    def _get_idx(self,subplot,L,hksets):
        d   = hksets[subplot]
        tmp = num.fabs(d['L']-L)
        idx = num.where(tmp==min(tmp))
        if len(idx) > 0:
            idx = idx[0]
        point_idx = d['point_idx'][idx]
        return point_idx

    ##########################################################################
    def get_scan(self,idx=None):
        """
        Get scan from point index (idx)

        If idx is None it uses the last plot selection        

        Returns:
        --------
        * tuple of (scan,point) where scan is the
          scan data object and point is the point
          in the scan (for image/stationary scans)
        """
        if idx == None:
            idx = self.get_idx()
            if idx == None: return None
        scan_idx = self.scan_index[idx]
        if len(scan_idx) > 1:
            point = scan_idx[1]
            idx = scan_idx[0]
        else:
            idx = scan_idx
            point = 0
        scan = self.scan[idx]
        return (scan,point)
    
    ##########################################################################
    def get_points(self,fig=None):
        """
        get index value of points from the plot zoom
        """
        if self.cursor == None:
            return None
        if self.cursor.clicked == False:
            #self.cursor.get_click(msg="Zoom on plot")
            return None
        self.cursor._zoom()
        z = self.cursor.zoom
        z = (z[0][0],z[1][0])
        Lmin = min(z)
        Lmax = max(z)
        subplot = self.cursor.subplot
        if subplot < 0:
            return None
        hksets  = sort_data(self)
        return self._get_idx_range(subplot,Lmin,Lmax,hksets)
        
    def _get_idx_range(self,subplot,Lmin,Lmax,hksets):
        d    = hksets[subplot]
        tmp  = d['L']
        idx1 = num.where(tmp<Lmin)
        idx2 = num.where(tmp>Lmax)
        tmp  = num.fabs(tmp)
        tmp[idx1] = 0.0
        tmp[idx2] = 0.0
        idx = num.where(tmp > 0 )
        point_idx = d['point_idx'][idx]
        return (point_idx)
    
    ##########################################################################
    def write_RSD(self,fname = None):
        """
        dump data file
        """
        if fname == None: fname = 'rasd'+str(int(round(self.H[0])))+str(int(round(self.K[0])))+'_'+str(round(self.L[0],3))+'.rsd'
        f = open(fname, 'w')
        header = "#%8s %8s %8s %8s %8s %8s %8s %8s\n" % ('Energy','H','K','L','F','Ferr','Alpha','Beta')
        f.write(header)
        for i in range(len(self.L)):
            if not i in self.bad:
                line = "%8.2f %8.2f %8.2f %8.3f %8.6g %8.6g %8.6g %8.6g\n" % (self.E[i],round(self.H[i]),
                                                                  round(self.K[i]),
                                                                  self.L[i],self.F[i],
                                                                  self.Ferr[i], self.scan[i].scalers['Alpha'][0],\
                                                                  self.scan[i].scalers['Beta'][0])
                f.write(line)
        f.close()

##########################################################################
##########################################################################
def sort_rasd(ctr,hkdecimal=2):
    """
    Return a dict of sorted data

    Assume H,K define a set with a range of L values
    All arrays should be of len npts. 

    Parameters:
    -----------
    * ctr is a ctr data object
    * hkdecimal is the number of precision to round H and K
      values to for sorting.  

    """
    # round H and K to sepcified precision
    H = num.around(ctr.H,decimals=hkdecimal)
    K = num.around(ctr.K,decimals=hkdecimal)
    L = num.around(ctr.L,decimals=hkdecimal)
    E = ctr.E
    F = ctr.F
    Ferr = ctr.Ferr
    #
    I     = ctr.I
    Inorm = ctr.Inorm
    Ierr  = ctr.Ierr
    Ibgr  = ctr.Ibgr
    #
    scan_idx = ctr.scan_index
    npts = len(F)

    #find all unique sets
    hkvals = []
    for j in range(npts):
        s = (H[j],K[j],L[j]) 
        if s not in hkvals:
            hkvals.append(s)

    # sort the hkvals
    # and stick data in correct set
    hkvals.sort()
    nsets = len(hkvals)
    #d = {'H':[],'K':[],'L':[],'F':[],'Ferr':[],'idx':[]}
    #hkset  = [copy.copy(d) for j in range(nsets)]
    hkset = []
    for j in range(nsets):
        hkset.append({'H':[],'K':[],'L':[],'E':[],'F':[],'Ferr':[],
                      'I':[],'Inorm':[],'Ierr':[],'Ibgr':[],
                      'point_idx':[],'scan_idx':[]})

    for j in range(npts):
        s      = (H[j],K[j],L[j])
        setidx = hkvals.index(s)
        hkset[setidx]['H'].append(H[j])
        hkset[setidx]['K'].append(K[j])
        hkset[setidx]['L'].append(L[j])
        hkset[setidx]['E'].append(E[j])
        hkset[setidx]['F'].append(F[j])
        hkset[setidx]['Ferr'].append(Ferr[j])
        hkset[setidx]['I'].append(I[j])
        hkset[setidx]['Inorm'].append(Inorm[j])
        hkset[setidx]['Ierr'].append(Ierr[j])
        hkset[setidx]['Ibgr'].append(Ibgr[j])
        hkset[setidx]['point_idx'].append(j)
        hkset[setidx]['scan_idx'].append(scan_idx[j])

    # make arrays num arrays
    for j in range(nsets):
        hkset[j]['H'] = num.array(hkset[j]['H'])
        hkset[j]['K'] = num.array(hkset[j]['K'])
        hkset[j]['L'] = num.array(hkset[j]['L'])
        hkset[j]['E'] = num.array(hkset[j]['E'])
        hkset[j]['F'] = num.array(hkset[j]['F'])
        hkset[j]['Ferr'] = num.array(hkset[j]['Ferr'])
        hkset[j]['I']     = num.array(hkset[j]['I'])
        hkset[j]['Inorm'] = num.array(hkset[j]['Inorm'])
        hkset[j]['Ierr']  = num.array(hkset[j]['Ierr'])
        hkset[j]['Ibgr']  = num.array(hkset[j]['Ibgr'])
        hkset[j]['point_idx']  = num.array(hkset[j]['point_idx'])
        hkset[j]['scan_idx']  = num.array(hkset[j]['scan_idx'])

    # now sort each set by E
    for j in range(nsets):
        lidx = num.argsort(hkset[j]['E'])
        hkset[j]['H'] = hkset[j]['H'][lidx]
        hkset[j]['K'] = hkset[j]['K'][lidx]
        hkset[j]['L'] = hkset[j]['L'][lidx]
        hkset[j]['E'] = hkset[j]['E'][lidx]
        hkset[j]['F'] = hkset[j]['F'][lidx]
        hkset[j]['Ferr'] = hkset[j]['Ferr'][lidx]
        hkset[j]['I'] = hkset[j]['I'][lidx]
        hkset[j]['Inorm'] = hkset[j]['Inorm'][lidx]
        hkset[j]['Ierr'] = hkset[j]['Ierr'][lidx]
        hkset[j]['Ibgr'] = hkset[j]['Ibgr'][lidx]
        hkset[j]['point_idx'] = hkset[j]['point_idx'][lidx]
        hkset[j]['scan_idx'] = hkset[j]['scan_idx'][lidx]

    return hkset
