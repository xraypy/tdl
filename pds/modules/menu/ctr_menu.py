"""
Menu function to handle interactive ctr data processing.

Authors / Modifications:
------------------------
Tom Trainor (tptrainor@alaska.edu)


"""
#######################################################################

import types
import copy
import numpy as num
from   matplotlib import pyplot

from   pds.lib.shellutil import Menu, show_more 
from   pds.lib.shellutil import get_str, get_tf, get_yn 
from   pds.lib.shellutil import get_int, get_flt, get_flt_list 
import ctr_data
import image_data
import image_menu

########################################################################

CTR_HEADER = """
##########################################
Number of points       = %s
New/selected points    = %s
Current selected point = %s
-----
Scan Type = '%s'
Labels: I='%s', Inorm='%s', Ierr='%s',Ibgr='%s'
Geom='%s'
Beam slits = %s
Det slits = %s
Sample  = %s
Scale = %s
Bad flag = %s,
H=%3.2f, K=%3.2f, L=%6.5f
I=%6.5g, Ierr=%6.5g, Ibgr=%6.5g, ctot=%6.5f
F=%6.5g, Ferr=%6.5g
------
"""

CTR_LABELS = ['plot','iplot','labels',
              'correction','cplot','integrate','copyint',
              'pointselect','pplot','zoomselect','selpoint',
              'next','previous','flag','help','done']
CTR_DESCR = ["Plot structure factors",
             "Plot intensities",
             "Set Intensity labels", 
             "Set Correction Params", 
             "Plot correction",
             "Set integration parameters / integrate",
             "Copy integration parameters",
             "Set 'selected point' from last plot click (Fig 0)",
             "Set+plot 'selected point' from last plot click (Fig 0)",
             "Set 'selected points' from last plot zoom (Fig 0)",
             "Select specific point number",
             "Select next point ",
             "Select previous point",
             "Toggle bad flag",
             "Show options",
             "All Done"]

########################################################################
def ctr_menu(ctr,scans=None,I=None,Inorm=None,Ierr=None,
             corr_params=None,scan_type=None):
    """
    Interactively inspect/integrate CtrData

    Parameters:
    -----------
    * ctr is a ctr data object
    * scans is a [list] of scan data objects
    * I, Inorm, Ierr are string tags corresponding to
      the values from the scandata object to used in computing
      structure factors
    * corr_params is a (optional) dictionary of correction parameters
    * scan_type is the type of scan (e.g. 'image')
    
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
        ctr.append_scans(scans,I=I,Inorm=Inorm,Ierr=Ierr,
                         corr_params=corr_params,scan_type=scan_type)
        set   = range(npts,len(ctr.L))
        point = npts
        npts  = len(ctr.L)
        
    # make menu
    m = Menu(labels=CTR_LABELS,descr=CTR_DESCR,sort=False,matchidx=True)

    # make a plot
    def _ctr_plot(ctr,I=False):
        #pyplot.close(0)
        pyplot.figure(0)
        pyplot.clf()
        if I == False:
            ctr.plot(fig=0,spnt=point)
        else:
            ctr.plot_I(fig=0,spnt=point)
    _ctr_plot(ctr)
    
    # loop
    while ret != 'done':
        if point in ctr.bad:
            bad_flag = True
        else:
            bad_flag = False
        header   = CTR_HEADER % (str(npts),str(set),str(point),
                                 str(ctr.scan_type[point]),
                                 str(ctr.labels['I'][point]),
                                 str(ctr.labels['Inorm'][point]),
                                 str(ctr.labels['Ierr'][point]),
                                 str(ctr.labels['Ibgr'][point]),
                                 str(ctr.corr_params[point].get('geom')),
                                 str(ctr.corr_params[point].get('beam_slits')),
                                 str(ctr.corr_params[point].get('det_slits')),
                                 str(ctr.corr_params[point].get('sample')),
                                 str(ctr.corr_params[point].get('scale')),
                                 str(bad_flag),
                                 ctr.H[point],ctr.K[point],ctr.L[point],
                                 ctr.I[point],ctr.Ierr[point],ctr.Ibgr[point],
                                 ctr.ctot[point],ctr.F[point],ctr.Ferr[point])
        m.header = header
        ret      = m.prompt(prompt)
        #
        if ret == 'plot':
            _ctr_plot(ctr)
        elif ret == 'iplot':
            _ctr_plot(ctr,I=True)
        elif ret == 'labels':
            (scan,spnt) = ctr.get_scan(point)
            lbls = scan.scalers.keys()
            if scan.image:
                lbls.extend(scan.image.peaks.keys())
            print "Labels: \n", lbls 
            Ilbl = get_str(prompt='Enter Intensity lbl',
                           default=ctr.labels['I'][point],valid=lbls)
            Ierrlbl = get_str(prompt='Enter Intensity error lbl',
                              default=ctr.labels['Ierr'][point],valid=lbls)
            Inormlbl = get_str(prompt='Enter Intensity norm lbl',
                               default=ctr.labels['Inorm'][point],valid=lbls)
            Ibgrlbl = get_str(prompt='Enter Intensity background lbl',
                               default=ctr.labels['Ibgr'][point],valid=lbls)
            pp = "Apply to 'p' current point, 's' to set, "
            pp = pp + "'a' to all, 'n' for none"
            app = get_str(prompt=pp,default='p',
                          valid=['p','a','s','n'])
            if app == 'p':
                ctr.labels['I'][point]     = Ilbl
                ctr.labels['Ierr'][point]  = Ierrlbl
                ctr.labels['Inorm'][point] = Inormlbl
                ctr.labels['Ibgr'][point]  = Ibgrlbl
                print "Integrate point index %i" % point
                ctr.integrate_point(point)
            elif app == 's':
                for j in set:
                    ctr.labels['I'][j]     = Ilbl
                    ctr.labels['Ierr'][j]  = Ierrlbl
                    ctr.labels['Inorm'][j] = Inormlbl
                    ctr.labels['Ibgr'][j]  = Ibgrlbl
                    print "Integrate point index %i" % j
                    ctr.integrate_point(j)
            elif app == 'a':
                for j in range(npts):
                    ctr.labels['I'][j]     = Ilbl
                    ctr.labels['Ierr'][j]  = Ierrlbl
                    ctr.labels['Inorm'][j] = Inormlbl
                    ctr.labels['Ibgr'][j]  = Ibgrlbl
                    print "Integrate point index %i" % j
                    ctr.integrate_point(j)
            elif app == 'n':
                print "Ignoring updates"
            _ctr_plot(ctr)
        elif ret == 'correction':
            corr_params = {}
            #
            check = get_yn(prompt='Edit geometry',default='y')
            if check=='yes':
                geom = get_str(prompt='Enter geometry',
                               default=ctr.corr_params[point]['geom'],
                               valid=['psic'])
                corr_params['geom'] = geom
            #
            check = get_yn(prompt='Edit beam slits',default='y')
            if check=='yes':
                beam_slits = ctr.corr_params[point]['beam_slits']
                if beam_slits == None: beam_slits = {'horz':1.,'vert':1.}
                print "Horizontal = total slit width in lab-z,or horizontal scatt plane"
                beam_slits['horz'] = get_flt('Beam horizontal slit (mm)',
                                             default=beam_slits['horz'])
                print "Vertical = total slit width in lab-x,or vertical scatt plane"
                beam_slits['vert'] = get_flt('Beam vertical slit (mm)',default=beam_slits['vert'])
                corr_params['beam_slits'] = beam_slits
            #
            check = get_yn(prompt='Edit detector slits',default='y')
            if check=='yes':
                det_slits = ctr.corr_params[point]['det_slits']
                if det_slits == None: det_slits = {'horz':1.,'vert':1.}
                print "Horizontal = total slit width in lab-z,or horizontal scatt plane"
                det_slits['horz'] = get_flt('Det horizontal slit (mm)',
                                            default=det_slits['horz'])
                print "Vertical = total slit width in lab-x,or vertical scatt plane"
                det_slits['vert'] = get_flt('Det vertical slit (mm)',
                                            default=det_slits['vert'])
                corr_params['det_slits'] = det_slits
            #
            check = get_yn(prompt='Edit sample shape/size',default='y')
            if check=='yes':
                sample = ctr.corr_params[point]['sample']
                check2 = get_str(prompt="Sample described by diameter 'd' or polygon 'p'?",
                                 default='d',valid=['d','p'])
                if check2 == 'd':
                    if type(sample) == types.DictType: sample = 1.
                    sample = get_flt('Sample diameter (mm)',default=sample)
                elif check2 == 'p':
                    xx = "Sample is described by a polygon which is a collection \n"
                    xx = xx + "of [x,y,(z)] points (minimum of three points);\n"
                    xx = xx + "x is lab frame vertical, y is along ki, z is horizontal.\n"
                    xx = xx + "These points need to be specified along with the set of\n"
                    xx = xx + "gonio angles at the time of the polygon determination.\n"
                    xx = xx + "If phi and chi are at the 'flat' values, and the sample is\n"
                    xx = xx + "on the beam center, then z = 0 (and only [x,y] points are needed).\n"
                    print xx
                    sample = {'polygon':[],'angles':{'phi':0.0,'chi':0.0}}
                    while 1:
                        xx = "Enter a sample polygon point, '[]' to end"
                        sp = get_flt_list(prompt=xx,default=[])
                        if sp == []:
                            break
                        else:
                            sample['polygon'].append(sp)
                    xx = "Enter gonio angles"
                    sample['angles']['phi'] = get_flt('phi',default=sample['angles']['phi'])
                    sample['angles']['chi'] = get_flt('chi',default=sample['angles']['chi'])
                    print 'Sample: ', sample
                corr_params['sample'] = sample
            #
            check = get_yn(prompt='Edit scale',default='y')
            if check=='yes':
                scale = get_flt(prompt='Enter scale factor',default=ctr.corr_params[point]['scale'])
                corr_params['scale'] = scale
            #
            pp = "Apply to 'p' current point, 's' to set, 'a' to all, 'n' for none"
            app = get_str(prompt=pp, default='p',
                          valid=['p','a','s','n'])
            if app == 'p':
                ctr.corr_params[point].update(corr_params)
                print "Integrate point index %i" % point
                ctr.integrate_point(point)
            elif app == 's':
                for j in set:
                    ctr.corr_params[j].update(corr_params)
                    print "Integrate point index %i" % j
                    ctr.integrate_point(j)
            elif app == 'a':
                for j in range(npts):
                    ctr.corr_params[j].update(corr_params)
                    print "Integrate point index %i" % j
                    ctr.integrate_point(j)
            elif app == 'n':
                print "Ignoring updates"
            #
            _ctr_plot(ctr)
        elif ret == 'cplot':
            (scan,spnt) = ctr.get_scan(point)
            corr_params = ctr.corr_params[point]
            corr = ctr_data._get_corr(scan,spnt,corr_params)
            if ctr.scan_type[point] == 'image':
                corr.ctot_stationary(plot=True,fig=2)
        elif ret == 'integrate':
            if ctr.scan_type[point] == 'image':
                (sc_idx,scan_pnt) = ctr.scan_index[point]
                scan = ctr.scan[sc_idx]
                image_menu.image_menu(scan,scan_pnt)
                #
                pp = "Apply / copy integration parameters to\n"
                pp = pp + "(note image max is not copied)\n"
                pp = pp + "'p' current point only, 's' to set, 'a' to all"
                app = get_str(prompt=pp,
                              default='p',valid=['p','a','s'])
                if app == 'p':
                    print "Integrate point index %i" % point
                    ctr.integrate_point(point)
                elif app == 's':
                    roi      = scan.image.rois[scan_pnt]
                    rotangle = scan.image.rotangle[scan_pnt]
                    bgrpar   = scan.image.bgrpar[scan_pnt]
                    for j in set:
                        print "Integrate point index %i" % j
                        ctr.integrate_point(j,roi=roi,rotangle=rotangle,bgrpar=bgrpar)
                elif app == 'a':
                    roi      = scan.image.rois[scan_pnt]
                    rotangle = scan.image.rotangle[scan_pnt]
                    bgrpar   = scan.image.bgrpar[scan_pnt]
                    for j in range(npts):
                        print "Integrate point index %i" % j
                        ctr.integrate_point(j,roi=roi,rotangle=rotangle,bgrpar=bgrpar)
            _ctr_plot(ctr)
        elif ret == 'copyint':
            idx = get_int(prompt='Enter point index to copy FROM',
                          default=point,min=0,max = npts-1)
            if ctr.scan_type[idx] == 'image':
                (sc_idx,scan_pnt) = ctr.scan_index[idx]
                scan = ctr.scan[sc_idx]
                #
                roi      = scan.image.rois[scan_pnt]
                rotangle = scan.image.rotangle[scan_pnt]
                bgrpar   = scan.image.bgrpar[scan_pnt]
                #
                pp = "Apply / copy integration parameters to\n"
                pp = pp + "(note image max is not copied)\n"
                pp = pp + "'p' current point only, 's' to set, 'a' to all"
                app = get_str(prompt=pp,
                              default='p',valid=['p','a','s'])
                if app == 'p':
                    print "Integrate point index %i" % point
                    ctr.integrate_point(point,roi=roi,rotangle=rotangle,bgrpar=bgrpar)
                elif app == 's':
                    for j in set:
                        print "Integrate point index %i" % j
                        ctr.integrate_point(j,roi=roi,rotangle=rotangle,bgrpar=bgrpar)
                elif app == 'a':
                    for j in range(npts):
                        print "Integrate point index %i" % j
                        ctr.integrate_point(j,roi=roi,rotangle=rotangle,bgrpar=bgrpar)
            _ctr_plot(ctr)
        elif ret == 'pointselect':
            idx = ctr.get_idx()
            point = idx[0]
        elif ret == 'pplot':
            idx = ctr.get_idx()
            point = idx[0]
            ctr.plot_point(idx=point,fig=3,show_int=True)
        elif ret == 'zoomselect':
            set = ctr.get_points()
            point = set[0]
        elif ret == 'selpoint':
            point = get_int(prompt='Enter point',
                            default=point,min=0,max = npts-1)
        elif ret == 'next':
            if point + 1 < npts: 
                point = point + 1
        elif ret == 'previous':
            if point - 1 > -1: 
                point = point - 1
        elif ret == 'flag':
            if point in ctr.bad:
                ctr.bad.remove(point)
            else:
                ctr.bad.append(point)
        else:
            pass

################################################################################
