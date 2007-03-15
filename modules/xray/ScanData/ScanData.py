# ScanData
# T. Trainor (fftpt@uaf.edu)
# Class to store and operate on generic 'scan data'
# Note spec.get_scan should ret one of these....
#################

from Num import Num

#################

"""
ScanData Class holds:
name = scan name
scan_dims -> would be [npts1] for a 1-d scan
            [npts1,npts2] for a 2-d scan
             where npts1 is the outer loop
             and npts1 is for the inner loop
mcas = []                    -> array has dims = scan_dims or is [] for No mca
scalers = {'io':[],'i1':[]}  -> each array has dims = scan_dims
positioners  = {'E':[]}      -> these may be single values, or same size as scalers
primary_axis = ['E']         -> if multi-dim, list of outer, inner
primary_det  = 'i1'          -> use for default plotting
"""

class ScanData:
    def __init__(self,name='',scan_dims=[],spectra=[],scalers={},
                 positioners={},primary_axis=[],primary_det=None):
        self.name         = name
        self.scan_dims    = scan_dims
        self.spectra      = spectra
        self.scalers      = scalers
        self.positioners  = positioners
        self.primary_axis = primary_axis
        self.primary_det  = primary_det

    def __repr__(self):
        lout = "Scan Name = %s\n" % self.name
        lout = lout + "Scan Dimensions = %s\n" % str(self.scan_dims)
        lout = lout + "Scalers:\n"
        for sc in self.scalers.keys():
            lout = lout + "    %s" % sc
        lout = lout + "\nPositioners:\n"
        for p in self.positioners.keys():
            lout = lout + "    %s" % p
        if self.spectra != []:
            lout = lout + "\nScan includes spectra data:\n"
        else:
            lout = lout + "\nScan does not include mca data:\n"
        lout = lout + "Primary scan axis = %s\n"  % str(self.primary_axis)
        lout = lout + "Primary detector = %s\n"  % self.primary_det
        return lout

    def get_scaler(self,label=None):
        if label == None:
            label = self.primary_det[0]
        return self.scalers.get(label)
        
    def get_positioner(self,label=None):
        if label == None:
            label = self.primary_axis[0]
        return self.positioners.get(label)

    def get_spectrum(self,pnt=0):
        pnt = int(pnt)
        if self.spectra == []:
            return None
        if pnt not in range(self.scan_dims[0]):
            return None
        return self.spectra[pnt]

    def get_spectra_OCR(self):
        ndet = self.spectra[0].ndet
        ocr = []
        for s in self.spectra:
            cts = s.get_count_totals()
            det_ocr = []
            for j in range(ndet):
                det_ocr.append(cts[j]['OCR'])
            ocr.append(det_ocr)
        return Num.transpose(ocr)

##############################################################################

