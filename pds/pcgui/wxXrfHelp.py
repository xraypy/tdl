########################################################################
"""
Tom Trainor (fftpt@uaf.edu)
This should be run as a child window from wxXrf

Modifications:
--------------


"""
########################################################################

from PythonCard import model, dialog
import wx
import os
import types

from   wxUtil import wxUtil

########################################################################

XrfParamHelp = """ XRF Help

==== BgrParams ====
** bottom_width      = 4.0
     Bottom width

** bottom_width_flag = 0
     Bottom width flag
     0 = Optimize
     1 = Fixed

** top_width         = 0.0
     Top width

** top_width_flag    = 0
     Top width flag
     0 = Optimize
     1 = Fixed

** exponent          = 2
     Exponent

** tangent           = False
     Tangent flag

** compress          = 4
     Compress

==== PeakParams ====

** label          = ""
     Peak label

** ignore         = False
     Don't fit peak

** energy         = 0.
     Peak energy

** energy_flag    = 1
     Flag for fitting energy
       0 = Optimize energy
       1 = Fix energy 

** fwhm           = 0.
     Peak FWHM

** fwhm_flag      = 1
     Flag for fitting FWHM
       0 = Optimize FWHM
       1 = Fix FWHM to global curve
       2 = Fix FWHM to input value

** ampl           = 0.
     Peak amplitude

** ampl_factor    = 0.
     Fixed amplitude ratio to previous peak
       0.  = Optimize amplitude of this peak
       >0. = Fix amplitude to this value relative
             to amplitude of previous free peak
      -1.0 = Fix amplitude at 0.0

==== FitParams ==== 

** energy_offset         =  0.
     Energy calibration offset (keV)

** energy_slope          =  1.
     Energy calibration slope 

** energy_flag           =  0
     Energy flag
       0 = Optimize energy calibration coefficients
       1 = Fix energy calibration coefficients

** fwhm_offset           =  0.15
     FWHM model offset (keV)

** fwhm_slope            =  0.
     FWHM model slope

** fwhm_flag             =  0
     Fwhm flag
       0 = Optimize FWHM coefficients
       1 = Fix FWHM coefficients

** chi_exp               =  0.
     Exponent of data for weights
       w = 1/y**chi_exp
       0.0 = ones
       0.5 = sqrt (poisson stats)

==== Fit Scan ====
** Fits all the xrf objects in the scan list

** Use Previous :
     If selected then use index-1 as seed parameters
     for fitting index

** Initialize:
     The index of the first scan in the list to fit
     It will then be used as the seed value for fitting
     all the xrf objects restarting at index zero.

     If Initialize > -1 and use previous is checked then
     the initial fit will only be used as the seed for fitting
     the first index.
       
"""

###########################################################
class wxXrfHelp(model.Background, wxUtil):

    ###########################################################
    def on_initialize(self, event):
        # Initialization
        # including sizer setup, do it here
        # self.setupSizers()
        self.help = XrfParamHelp
        self.components.HelpText.text = self.help

    ###########################################################
    def on_menuFileExit_select(self,event): 
        self.close()
