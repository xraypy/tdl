# T. Trainor (2006)
# Wrapper functions for ScanData
# T. Trainor, 5-16-2007
#
# --------------
# Modifications
# --------------
#
#
# --------------
# Todo
# --------------
import Util
from Num import Num
#import SD as ScanData
import ScanData
import XRF
import deadtime
import types

####################################################
def get_scaler(sd,label=None,**kws):
    return sd.get_scaler(label=label)
    
def get_positioner(sd,label=None,**kws):
    return sd.get_positioner(label=label)
    
def get_spectrum(sd,pnt,**kws):
    return sd.get_spectrum(pnt)

def get_spectra_OCR(sd,**kws):
    return sd.get_spectra_OCR()

def fit_deadtime(sd,io=None,plot=True,tdl=None,**kws):
    ocr = sd.get_spectra_OCR()
    if io == None:
        io = sd.get_scaler()
    elif type(io) == types.StringType:
        io = sd.get_scaler(label=io)
    io = Num.array(io)
    #ocr = Num.array(ocr)
    tau = []
    for j in range(len(ocr)):
        (params,msg) = deadtime.fit_deadtime_curve(io,ocr[j])
        a = params[0]
        t = params[1]
        tau.append(t)
        # test
        if plot:
            ocr_calc = deadtime.calc_ocr(params,io)
            tdl.setVariable('xocr_c',val=ocr_calc)
            tdl.setVariable('xocr',val=ocr[j])
            tdl.setVariable('xio',val=io)
            cmd = "plot(xio,xocr,fmt='o')"
            tdl.eval(cmd)
            cmd = "plot(xio,xocr_c,fmt='r')"
            tdl.eval(cmd)
            
    return tau


#################################################
_groups_ = [('scan',True)]

_func_ = {"scan.get_scaler":get_scaler,
          "scan.get_positioner":get_positioner,
          "scan.get_spectrum":get_spectrum,
          "scan.get_ocr":get_spectra_OCR,
          "scan.fit_deadtime":fit_deadtime}

#_scripts_ = ['']
