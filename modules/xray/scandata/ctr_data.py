##############################################################################
"""
T. Trainor (tptrainor@alaska.edu)
Frank Heberling (Frank.Heberling@ine.fzk.de)
Functions for extracting ctr data from ScanData objects

Modifications:
--------------
"""
##############################################################################
"""
Notes
------

References:
Vlieg 1997
Schlepuetz 2005


Todo
----
- Test!
- sorting and specify HK for plots
- averaging/merging and merge statistics
- corrections for image and rocking scans
"""
#############################################################################

import types
import numpy as num
from matplotlib import pyplot

from mathutil import cosd, sind, tand
from mathutil import arccosd, arcsind, arctand
from xtal.active_area import active_area
import gonio_psic 
import active_area_psic

##############################################################################
class CtrData:
    def __init__(self,scans=[],I_lbl='I_c',Ierr_lbl='Ierr_c',AAparams=[]):
        #
        self.I_lbl    = I_lbl
        self.Ierr_lbl = Ierr_lbl
        self.AAparams = AAparams
        #
        self.H    = num.array([],float)
        self.K    = num.array([],float)
        self.L    = num.array([],float)
        self.I0   = num.array([],float)
        self.I    = num.array([],float)
        self.Ierr = num.array([],float)
        self.corr = num.array([],float)
        self.F    = num.array([],float)
        self.Ferr = num.array([],float)
        #
        self.append_data(scans)


    ##########################################################################
    def get_data_from_scans(self,scans,type='Image'):
        # data = get_data_from_scans(scans)
        # this should parse the scans into arrays
        # data['H'], data['K'],data['L'],data['Io'],
        # data['I'], data['Ierr'], data['corr']
        # and data['idx']
        # Note this should append scandata instances
        # to a list self.scandata (which should write images
        # to hdf then delete)
        # data['idx'] should point to this, ie
        # self.scan_data[idx] -> returns original scan data object
        # of a raw point.  
        pass
    

    ##########################################################################
    def append_data(self,scans):
        """
        scans is a list of scan data objects
        """
        if type(scans) != types.ListType:
            scans = [scans]

        # data = get_data_from_scans(scans)
        # this should parse the scans into arrays
        # data['H'], data['K'],data['L'],data['Io'],
        # data['I'], data['Ierr'], data['corr']
        # and data['idx']
        # then below appends are basically fine

        for scan in scans:
            #fill dictionary with data from all scans
            self.H    = num.append(self.H,scan.scalers['H'])
            self.K    = num.append(self.K,scan.scalers['K'])
            self.L    = num.append(self.L,scan.scalers['L'])
            self.I0   = num.append(self.I0,scan.scalers['io'])
            self.I    = num.append(self.I,scan.image_peaks[self.I_lbl])
            self.Ierr = num.append(self.Ierr,scan.image_peaks[self.Ierr_lbl])
            # the correction Factor is multiplied by 17000 to rescale to Counts/sec

            # work here!!!
            #-->norm by Io: scan.scalers['io']
            # then compute correction 
            corr = self.calc_correction(scan)*17000 
            F    = (corr*scan.image_peaks[self.I_lbl])**0.5
            Ferr = (corr*scan.image_peaks[self.Ierr_lbl])**0.5
            #
            self.corr = num.append(self.corr,corr)
            self.F    = num.append(self.F,F)
            self.Ferr = num.append(self.Ferr,Ferr)

    ##########################################################################
    def plot(self):
        """
        plot the rod.
        Need multi panel plot 
        """
        pyplot.figure()
        pyplot.plot(self.L, self.F,'b')
        pyplot.errorbar(self.L,self.F,self.Ferr, fmt ='o')
        pyplot.semilogy()
        pyplot.show()

    ##########################################################################
    def write_HKL(self,fname = 'ctr.hkl'):
        #write data file for row data
        f = open(fname, 'w')
        f.write('#idx        H            K            L            F           F_err \n')
        for i in range(len(self.L)):
            if self.I[i] > 0:
                line = "%i    %i    %i    %6.3f    %6.6g    %6.6g\n" % (i,
                                                                  round(self.H[i]),
                                                                  round(self.K[i]),
                                                                  self.L[i],
                                                                  self.F[i],
                                                                  self.Ferr[i])
                f.write(line)
        f.close()

##############################################################################
def append_ctr(scans,ctr=None,I_lbl='I_c',Ierr_lbl='Ierr_c',AAparams=[]):
    if ctr == None:
        return CtrData(scans=scans,I_lbl=I_lbl,Ierr_lbl=Ierr_lbl,AAparams=AAparams)
    else:
        return ctr.append_data(scans)

##############################################################################
##############################################################################
class CtrCorrectionPsic:
    """
    Data point operations / corrections for Psic geometry

    Note: All correction factors are defined such that the
    measured data is corrected by multiplying times
    the correction: 
      Ic  = Im*ct
    where
      Im = Idet/Io = uncorrected (measured) intensity

    In other words we use the following formalism:
      Im = (|F|**2)* prod_i(Xi)
    where Xi are various (geometric) factors that 
    influence the measured intensity.  To get the
    structure factor:
      |F| = sqrt(Im/prod_i(Xi)) = sqrt(Im* ct)
    and
      ct = prod_i(1/Xi) = prod_i(ci)
      ci = 1/Xi
      
    If there is an error or problem in the routine for a sepcific
    correction factor, (e.g. divide by zero), the routine should
    return a zero.  This way the corrected data is zero'd....

    The correction factors depend on the goniometer geometry
     gonio = gonio_psic.Psic instance

    The slits settings are needed.  Note if using a large area detector
    you may pass det_slits = None and just spill off will be computed
       beam_slits = {'horz':.6,'vert':.8}
       det_slits = {'horz':20.0,'vert':10.5}
      these are defined wrt psic phi-frame:
       horz = beam/detector horz width (total slit width in lab-z,
              or the horizontal scattering plane)
       vert = detector vert hieght (total slit width in lab-x,
              or the vertical scattering plane)

    A sample description is needed.
    If sample = a number then is is taken as the diameter
    of a round sample mounted on center.

    Otherwise is may describe a general sample shape:    
      sample = {}
        sample['polygon'] = [[1.,1.], [.5,1.5], [-1.,1.],
                             [-1.,-1.],[0.,.5],[1.,-1.]]
        sample['angles']  = {'phi':108.0007,'chi':0.4831}

        polygon = [[x,y,z],[x,y,z],[x,y,z],....]
                  is a list of vectors that describe the shape of
                  the sample.  They should be given in general lab
                  frame coordinates.

         angles = {'phi':0.,'chi':0.,'eta':0.,'mu':0.}
             are the instrument angles at which the sample
             vectors were determined.

      The lab frame coordinate systems is defined such that:
        x is vertical (perpendicular, pointing to the ceiling of the hutch)
        y is directed along the incident beam path
        z make the system right handed and lies in the horizontal scattering plane
          (i.e. z is parallel to the phi axis)

        The center (0,0,0) of the lab frame is the rotation center of the instrument.

        If the sample vectors are given at the flat phi and chi values and with
        the correct sample hieght (sample Z set so the sample surface is on the
        rotation center), then the z values of the sample vectors will be zero.
        If 2D vectors are passed we therefore assume these are [x,y,0].  If this
        is the case then make sure:
            angles = {'phi':flatphi,'chi':flatchi,'eta':0.,'mu':0.}

        The easiest way to determine the sample coordinate vectors is to take a picture
        of the sample with a camera mounted such that is looks directly down the omega
        axis and the gonio angles set at the sample flat phi and chi values and
        eta = mu = 0. Then find the sample rotation center and measure the position
        of each corner (in mm) with up being the +x direction, and downstream
        being the +y direction.  
    
    """
    def __init__(self,gonio=None,beam_slits={},det_slits={},sample={}):
        self.gonio      = gonio
        if self.gonio.calc_psuedo == False:
            self.gonio._update_psuedo()
        self.beam_slits = beam_slits
        self.det_slits  = det_slits
        self.sample     = sample
        # fraction horz polarization
        self.fh         = 1.0
        self.scale      = 100.

    ##########################################################################
    def ctot_stationary(self,plot=False):
        """
        correction factors for stationary measurements (e.g. images)
        """
        cp = self.polarization()
        cl = self.lorentz_stationary()
        ca = self.active_area(plot=plot)
        ct = (self.scale)*(cp)*(cl)*(ca)
        if plot == True:
            print "Polarization=%f" % cp
            print "Lorentz=%f" % cl
            print "Area=%f" % ca
            print "Scale=%f" % self.scale
            print "Total=%f" % ct
        return ct

    ##########################################################################
    def lorentz_stationary(self):
        """
        Compute the Lorentz factor for a stationary (image)
        measurement.  See Vlieg 1997

        Measured data is corrected for Lorentz factor as: 
          Ic  = Im * cl
        """
        beta  = self.gonio.pangles['beta']
        cl = sind(beta)
        return cl

    ##########################################################################
    def lorentz_scan(self):
        """
        Compute the Lorentz factor for a generic scan
        Note this is approximate. In general for bulk xrd
        with single crystals lorentz factor is defined as:
          L = 1/sin(2theta)
        We need L for specific scan types, e.g. phi, omega, etc..
        See Vlieg 1997

        Measured data is corrected for Lorentz factor as: 
          Ic  = Im * cl = Im/L
        """
        tth  = self.gonio.pangles['tth']
        cl = sind(tth)
        return cl

    ##########################################################################
    def rod_intercept(self,):
        """
        Compute the dl of the rod intercepted by the detector.
        (this correction only applies for rocking scans)
        This can be (roughly) approximated from:
          dl = dl_o * cos(beta)
          where dl_o is the detector acceptance at beta = 0

        Note this is an approximation for all but simple specular scans,
        Should use the detector acceptance polygon and determine
        the range of dl for the specific scan axis used.
        
        Measured data is corrected for the intercept as: 
          Ic  = Im * cr = Im/dl
        """
        beta  = self.gonio.pangles['beta']
        cr   = cosd(beta)
        return cr
    
    ##########################################################################
    def polarization(self,):
        """
        Compute polarization correction factor.
        
        For a horizontally polarized beam (polarization vector
        parrallel to the lab-frame z direction) the polarization
        factor is normally defined as:
           p = 1-(cos(del)*sin(nu))^2
        For a beam with mixed horizontal and vertical polarization:
           p = fh( 1-(cos(del)*sin(nu))^2 ) + (1-fh)(1-sin(del)^2)
        where fh is the fraction of horizontal polarization.

        Measured data is corrected for polarization as: 
          Ic  = Im * cp = Im/p
        """
        fh    = self.fh
        delta = self.gonio.angles['delta']
        nu    = self.gonio.angles['nu']
        p = 1. - ( cosd(delta) * sind(nu) )**2.
        if fh != 1.0:
            p = fh * c_p + (1.-fh)*(1.0 - (sind(delta))**2.)
        if p == 0.:
            cp = 0.
        else:
            cp = 1./p

        return cp

    ##########################################################################
    def active_area(self,plot=False):
        """
        Compute active area correction (c_a = A_beam/A_int)
        Use to correct scattering data for area effects,
        including spilloff, i.e.
           Ic = Im * ca = Im/A_ratio 
           A_ratio = A_int/A_beam 
        where
           A_int = intersection area (area of beam on sample
                   viewed by detector)
           A_beam = total beam area
        """
        alpha = self.gonio.pangles['alpha']
        beta  = self.gonio.pangles['beta']
        if plot == True:
            print 'alpha = ', alpha, ', beta = ', beta
        if alpha < 0.0:
            print 'alpha is less than 0.0'
            return 0.0
        elif beta < 0.0:
            print 'beta is less than 0.0'
            return 0.0
        
        # get beam vectors
        bh = self.beam_slits['horz']
        bv = self.beam_slits['vert']
        beam = gonio_psic.beam_vectors(h=bh,v=bv)
        # get det vectors
        dh = self.det_slits['horz']
        dv = self.det_slits['vert']
        det  = gonio_psic.det_vectors(h=dh,v=dv,
                                      nu=self.gonio.angles['nu'],
                                      delta=self.gonio.angles['delta'])
        # get sample poly
        if type(self.sample) == types.DictType:
            sample_vecs = self.sample['polygon']
            sample_angles = self.sample['angles']
            sample = gonio_psic.sample_vectors(sample_vecs,
                                               angles=sample_angles,
                                               gonio=self.gonio)
        else:
            sample = self.sample
        # compute active_area
        (A_beam,A_int) = active_area_psic.active_area(self.gonio.nm,
                                                      ki=self.gonio.ki,
                                                      kr=self.gonio.kr,
                                                      beam=beam,det=det,
                                                      sample=sample,plot=plot)
        if A_int == 0.:
            ca = 0.
        else:
            ca = A_beam/A_int

        return ca

##############################################################################
##############################################################################
def test1():
    import gonio_psic
    #psic = psic_from_spec(G,angles={})
    psic = gonio_psic.test2(show=False)
    psic.set_angles(phi=12.,chi=30.,eta=20.,
                    mu=25.,nu=75.,delta=20.)
    #print psic
    #
    beam_slits = {'horz':.6,'vert':.8}
    det_slits = {'horz':20.0,'vert':10.5}
    sample = {}
    sample['polygon'] = [[1.,1.], [.5,1.5], [-1.,1.], [-1.,-1.],[0.,.5],[1.,-1.]]
    sample['angles']  = {'phi':108.0007,'chi':0.4831}
    #
    cor = CtrCorrectionPsic(gonio=psic,beam_slits=beam_slits,
                            det_slits=det_slits,sample=sample)
    ct = cor.ctot_stationary(plot=True)

##########################################################################
if __name__ == "__main__":
    """
    test 
    """
    test1()
    #test2()


