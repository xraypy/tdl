"""
This module defines a device-independent MultiChannel Analyzer (MCA) class.
Original written by Mark Rivers, GSECARS
See http://cars9.uchicago.edu/software/python/index.html

--------------
 Modifications
--------------
- Modified for Tdl, tpt

"""
########################################################################

import numpy as Num
import math
import types
import DeadTime as deadtime
import McaCalib as calib

########################################################################
class Mca:
    """
    MultiChannel Analyzer (MCA) class 
        
    Fields:
        self.name        = 'mca'  # Name of the mca object
        self.nchans      = 2048   # number of mca channels
        self.data        = None   # MCA data 
        self.channels    = None   # MCA channels value 

        # Counting parameters
        self.start_time   = ''    # Start time and date, a string
        self.live_time    = 0.    # Elapsed live time in seconds
        self.real_time    = 0.    # Elapsed real time in seconds
        self.read_time    = 0.    # Time that the Mca was last read in seconds
        self.total_counts = 0.    # Total counts between the preset start and stop channels
        self.input_counts = 0.    # Actual total input counts (eg given by detector software)
                                  # Note total_counts and input counts ARE NOT time normalized
                                  #
        self.tau          = -1.0  # Factor for deadtime/detector saturation calculations, ie
                                  #     ocr = icr * exp(-icr*tau)
        self.icr_calc     = -1.0  # Calculated input count rate from above expression
        self.cor_factor   = 1.0   # Calculated correction factor based on icr,ocr,lt,rt
                                  # data_corrected = data * cor_factor
                                  #
        # Calibration parameters
        self.offset       = 0.    # Offset
        self.slope        = 1.0   # Slope
        self.quad         = 0.    # Quadratic
        self.units        = 'keV' # Calibration units, a string
        self.two_theta    = 0.    # 2-theta of this Mca for energy-dispersive diffraction

    Notes:
    If input_counts are read from the data file then there is no need for tau
    values (and hence icr_calc is not used) - ie we assume the detector
    provides the icr value. 

    The application of corrections uses the following logic:
       if tau > 0 this will be used in the correction factor calculation 
       if tau = 0 then we assume ocr = icr in the correction factor calculation
       if tau < 0 (or None):
          if input_counts > 0 this will be used for icr in the factor calculation
          if input_counts <= 0 we assume ocr = icr in the correction factor calculation

    Energy calibration is based on the following:
        energy = offset + slope*channel + quad*channel**2

    If channels are not explicitly given, we'll assume:
       channels = Num.arange(NCHAN,dtype=int)
    with NCHAN = 2048
    """
    ###############################################################################
    def __repr__(self):
        lout = "  **Detector %s, nchans = %i\n" % (self.name, self.nchans)
        lout = lout + "    >RealTime= %6.5f, LiveTime= %6.5f\n" % (self.real_time, self.live_time)
        lout = lout + "    >OCR     = %6.1f, ICR     = %6.1f\n" % (self.total_counts / self.live_time,
                                                                     self.input_counts / self.live_time)
        lout = lout + "    >Tau     = %6.5f, ICR_Calc= %6.1f\n" %  (self.tau, self.icr_calc)
        lout = lout + "    >Correct = %6.5f, Deadtime =%5.2f%%\n" %  (self.cor_factor,
                                                                                    100.0*(1 - 1/self.cor_factor))
        return lout

    ###############################################################################
    def __init__(self,data=None,channels=None, **kws):
        """ Initialize the MCA structure """
        self.det_type    = "MCA"
        self.name        = 'mca'  # Name of the mca object
        self.nchans      = 2048   # number of mca channels
        self.data        = []     # MCA data 
        self.channels    = []     # MCA channels value 

        # Counting parameters
        self.start_time   = ''    # Start time and date, a string
        self.read_time    = 0.    # Time that the Mca was last read in seconds
                                  # start_time and read_time are not used
                                  # for any corrections
        self.real_time    = 0.    # Elapsed real time in seconds (requested counting time)
        self.live_time    = 0.    # Elapsed live time in seconds (time detector is live)
        self.total_counts = 0.    # Total counts between the preset start and stop channels
        self.input_counts = -1.0  # Actual total input counts (eg given by detector software)
                                  # Note total_counts and input counts ARE NOT time normalized
                                  #
        self.tau          = -1.0  # Factor for deadtime/detector saturation calculations, ie
                                  #     ocr = icr * exp(-icr*tau)
        # Calculated correction values
        self.icr_calc     = -1.0  # Calculated input count rate from above expression
        self.cor_factor   = 1.0   # Calculated correction factor based on icr,ocr,lt,rt
                                  # data_corrected = data * cor_factor
                                  #
        # Calibration parameters
        self.offset       = 0.    # Offset
        self.slope        = 1.0   # Slope
        self.quad         = 0.    # Quadratic
        self.units        = 'keV' # Calibration units, a string
        self.two_theta    = 0.    # 2-theta of this Mca for energy-dispersive diffraction

        # Init data and params
        self.init_params(mca_params=kws)
        self.init_data(data=data,channels=channels) 

    ###############################################################################
    def init_params(self, mca_params={}):
        """
        set/reset parameters based on key word arguments
        """
        for key in mca_params.keys():
            setattr(self,key,mca_params[key])
        
        # Make sure correction updated
        self._calc_correction()

    ########################################################################
    def get_calib_params(self,):
        """
        return calibration data
        """
        ret = {'offset':self.offset,
               'slope':self.slope,
               'quad':self.quad,
               'units':self.units,
               'tth':self.two_theta}
        return ret

    ########################################################################
    def init_data(self, data=None, channels=None):
        """
        Inputs:
            data: A numpy array of data (counts).
            channels: Array of channel numbers
        """
        if data:
            self.data = Num.asarray(data,dtype=Num.int)
        elif self.data==[]:
            self.data = Num.zeros(self.nchans,dtype=Num.int)
            
        # Note if channels == None, assume the same 
        # length as data and channel[0] = 0
        if channels:
            self.channels = Num.asarray(channels,dtype=Num.int)
        else:
            self.channels = Num.arange(len(self.data),dtype=Num.int)
        
        # Check
        self.nchans = len(self.data)
        if len(self.channels) != self.nchans:
            raise "Data-channel length mismatch in MCA %s" % self.name
        return

    ########################################################################
    def update_correction(self,tau=None):
        """
        Update the deadtime correction 
        if tau == None just recompute,
        otherwise assign a new tau and recompute
        """
        if tau != None:
            self.tau = tau
        self._calc_correction()

    ########################################################################
    def _calc_correction(self):
        """
        if self.tau > 0 this will be used in the correction factor calculation 
        if self.tau = 0 then we assume ocr = icr in the correction factor calculation,
                      ie only lt correction
                     (note deadtime.calc_icr handles above two conditions)
        if self.tau < 0 (or None):
           if input_counts > 0  this will be used for icr in the factor calculation
           if input_counts <= 0 we assume ocr = icr in the correction factor calculation,
                                ie only lt correction
        """
        if (self.live_time <=0) or (self.real_time <=0):
            self.cor_factor  = 1.0
            return
        
        if self.total_counts > 0:
            ocr = self.total_counts / self.live_time
        else:
            ocr = None
        
        if self.tau >= 0:
            self.icr_calc = deadtime.calc_icr(ocr,self.tau)
            icr = self.icr_calc
        elif self.input_counts > 0:
            icr = self.input_counts / self.live_time
        else:
            icr = ocr = None
        self.cor_factor  = deadtime.correction_factor(self.real_time,
                                                      self.live_time,
                                                      icr = icr,
                                                      ocr = ocr)
        if self.cor_factor <= 0:
            print "Error computing data correction factor --> setting to 1"
            self.cor_factor = 1.0

    ########################################################################
    def get_data(self,correct=True):
        """
        Returns the data (counts) from the Mca
        Note if correct == True the corrected data is returned. However,
        this does not (re)compute the correction factor, therefore, make
        sure the correction factor is up to date before requesting corrected data...
        """
        if correct == True:
            d = self.cor_factor * self.data
            # note adding .5 rounds the data
            d = (d+0.5).astype(Num.int)
            return d
        else:
            return self.data

    ########################################################################
    def get_energy(self):
        """
        Returns a list containing the energy of each channel in the MCA spectrum.
        """
        energy = calib.channel_to_energy(self.channels,
                                         offset=self.offset,
                                         slope=self.slope,
                                         quad=self.quad)
        return energy
    
    ####################################################################
    """
    def count_totals(self):
        ""
        useful for doing a deadtime fit???
        ""
        rt = self.elapsed.real_time
        lt = self.elapsed.live_time
        OCR = self.elapsed.total_counts/mca.elapsed.live_time
        ICR = self.elapsed.input_counts/mca.elapsed.live_time
        ICR_CALC = self.elapsed.icr_calc
        COR = self.elapsed.cor_factor
        lst = {'rt':lt,'lt':lt,'OCR':OCR,'ICR':ICR,'ICR_CALC':ICR_CALC,'COR':COR}
        return lst
    """
    #########################################################################
            