#######################################################################
"""
Tom Trainor (tptrainor@alaska.edu)

Menu function to handle interactive
ctr processing

Modifications:
--------------


"""
#######################################################################

import types
import copy
import numpy as num
from   matplotlib import pyplot

from   pds.shellutil import Menu, show_more 
from   pds.shellutil import get_tf, get_yn, get_int, get_flt
from   plotter import cursor
import ctr_data
import image_data
import data

########################################################################
IMG_HEADER = """
Number of points  = %s
Current point     = %s
"""

IMG_LABELS = ['display','imax','rotangle','setroi','plotsums',
              'selectroi','bgr','copyall','integrate','intall',
              'point','next','previous','flag','help','quit']
IMG_DESCR = ["Display image",
             "Set max image intensity value",
             "Set image rotation angle (deg ccw)",
             "Set roi from image zoom (Figure 1)",
             "Plot row/column sums (Figure 2)",
             "Select roi from sum plots (Figure 2)",
             "Set background parameters",
             "Apply current roi and background params to all images",
             "Integrate current image",
             "Integrate all images",
             "Select scan point",
             "Select next point ",
             "Select previous point", 
             "Flag as bad point",
             "Show options",
             "Quit / All Done"]

"""
NEED:
- append new stuff if passed and find idx range of new stuff.
- integrate:
   - a point (idx)
   - range of points (idx or new)
   - selection from plot (single or range)
   (note in integrate menu operate on range of points...)
- modify params I, Inorm, corr_params etc... (apply to point, all, new)
- plot correction (e.g. area)
"""

########################################################################
def ctr_menu(ctr,scans=[],I='I_c',Inorm='io',Ierr='Ierr_c',
             corr_params={},scan_type='image'):
    """
    Interactively inspect/integrate CtrData 
    """
    if not isinstance(ctr,ctr_data.CtrData):
        print "Error CtrData object required"
        return
    prompt   = 'Select option >'
    npts     = len(ctr.L)
    scan_pnt = 0
    ret      = ''

    # make menu
    m = Menu(labels=IMG_LABELS,descr=IMG_DESCR,sort=False,matchidx=True)
    
    # loop
    while ret != 'quit':
        header   = IMG_HEADER % (str(npts),str(scan_pnt))
        m.header = header
        ret      = m.prompt(prompt)

        if ret == 'display':
            ctr.plot()
        else:
            pass

################################################################################

