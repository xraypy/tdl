"""
Data analysis package

This package contains modules used for
processing/analyzing x-ray data. More specifically
the modules are geared towards processing surface
diffraction, reflectivity and reflection
standing wave (and other related) data.  
"""
from reader       import Reader
from reader       import spec_scan as spec
from scan_data    import ScanData 
from ctr_data     import CtrData
from rasd_data    import RasdData

