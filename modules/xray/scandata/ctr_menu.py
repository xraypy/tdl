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
Number of points       = %s
New/selected points    = %s
Current selected point = %s
-----
Scan Type = '%s'
I='%s', Inorm='%s',Ierr='%s'
-----
Geom='%s'
Beam slits = %s
Det slits = %s
Sample  = %s
Scale = %s
------
Add some HKL, I F etc...
"""

IMG_LABELS = ['plot','pointselect','zoomselect','labels','correction',
              'cplot','setint','integrate','intselected','intall',
              'point','next','previous','help','quit']
IMG_DESCR = ["Plot all structure factors",
             "Select scan point from plot (Fig 0)",
             "Select scan set from plot zoom (Fig 0)",
             "Set Intensity labels", # --> sub menu apply to current, set, all
             "Set Correction Params", # "" 
             "Plot correction",
             "Set integration parameters", # "", also set bgr and flag bad
             "Integrate current point",
             "Integrate selected points",
             "Integrate all points",
             "Select scan point(s)",
             "Select next point ",
             "Select previous point", 
             "Show options",
             "Quit / All Done"]

########################################################################
def ctr_menu(ctr,scans=None,I=None,Inorm=None,Ierr=None,
             corr_params=None,scan_type=None):
    """
    Interactively inspect/integrate CtrData 
    """
    #if not isinstance(ctr,ctr_data.CtrData):
    #    print "Error CtrData object required"
    #    return
    prompt   = 'Select option >'
    ret      = ''
    point    = 0
    set      = []
    npts     = len(ctr.L)

    # append new data
    if scans!=None:
        ctr.append_scans(self,scans,I=I,Inorm=Inorm,Ierr=Ierr,
                         corr_params=corr_params,scan_type=scan_type)
        set   = [npts,len(ctr.L)]
        point = npts
        npts  = len(ctr.L)
        
    # make menu
    m = Menu(labels=IMG_LABELS,descr=IMG_DESCR,sort=False,matchidx=True)
    
    # loop
    while ret != 'quit':
        header   = IMG_HEADER % (str(npts),str(set),str(point),
                                 str(ctr.scan_type[point]),
                                 str(ctr.labels['I'][point]),
                                 str(ctr.labels['Inorm'][point]),
                                 str(ctr.labels['Ierr'][point]),
                                 str(ctr.corr_params[point].get('geom')),
                                 str(ctr.corr_params[point].get('beam_slits')),
                                 str(ctr.corr_params[point].get('det_slits')),
                                 str(ctr.corr_params[point].get('sample')),
                                 str(ctr.corr_params[point].get('scale')),
                                 )
        m.header = header
        ret      = m.prompt(prompt)

        if ret == 'plot':
            ctr.plot(fig=0)
        elif ret == 'point':
            point = get_int(prompt='Enter point',
                            default=point,min=0,max = npts-1)
        elif ret == 'next':
            if point + 1 < npts: 
                point = point + 1
        elif ret == 'previous':
            if point - 1 > -1: 
                point = point - 1
        else:
            pass


################################################################################

