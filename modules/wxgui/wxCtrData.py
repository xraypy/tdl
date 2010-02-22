"""
CTR Gui

Authors/Modifications:
----------------------
Tom Trainor (tptrainor@alaska.edu)


"""
###############################################################################

from PythonCard import model, dialog
import wx
import os
import glob
import types
import copy
import time
import numpy as num
from matplotlib import pyplot

from wxUtil import wxUtil
import scandata
import ctr_data
import image_data

###############################################################################

CTR_HEADER = """Number of points          = %s
New/selected points (set) = %s
Current selected point    = %s
H=%3.2f, K=%3.2f, L=%6.5f
Scan Type = '%s'
Labels: I='%s', Inorm='%s', Ierr='%s', Ibgr='%s'
Geom='%s'
Beam slits = %s
Det slits  = %s
Sample  = %s
Scale   = %s
Bad point flag = %s,
I=%6.5g, Ierr=%6.5g, Ibgr=%6.5g, ctot=%6.5f
F=%6.5g, Ferr=%6.5g
"""

PARAM_DESCR = {
'I':"""Enter label for intensity. Options:\n%s""",

'Inorm':"""Enter label for intensity normalization. Options:\n%s""",

'Ierr':"""Enter label for intensity errors. Options:\n%s""",

'Ibgr':"""Enter label for intensity background. Options:\n%s""",

'image roi':"""Enter the image roi. The format is [x1,y1,x2,y2] where values
are in pixels corresponding to two corners of the box. Use the set button to
select the roi from the image plot zoom value (fig 2) -> raw scan data plot.""",

'image rotangle':"""Enter a rotation angle for the image (counter clockwise degrees)""",

'bgr flag':"""Enter flag for image background method:
    0 determine row and column backgrounds after summation
    1 determine 2D background using 'c'olumn direction 
    2 determine 2D background using 'r'ow direction
    3 determine 2D background from the average 'r'ow and 'c'olumn directions""",

'bgr col nbgr':"""Number of background points for linear part of column direction
(y-direction) bgr fit. If nbgr = 0, no linear fit is included""",

'bgr col width':"""Peak width for the column (y) direction bgr fit.
The background function should fit features that are in general broader
than the width value. Estimate cwidth using width from peak width in the
column (y) direction. Note: width = 0 corresponds to no polynomial bgr""",

'bgr col power':"""Power of polynomial used in column (y) direction bgr fit.
Use pow < 0.5 for polynomials flatter than a circle (pow = 0 are linear)
and pow  > 0.5 for steeper.""",

'bgr col tan':"""Flag for use of tangents in column (y) direction bgr determination.
If 'True' then the local slope is removed when fitting a polynomial to a point.
Helpful for data sitting on a broad sloping background""",

'bgr row nbgr':"""Number of background points for linear part of row direction
(x-direction) bgr fit. If nbgr = 0, no linear fit is included""",

'bgr row width':"""Peak width for the row (x) direction bgr fit.
The background function should fit features that are in general broader
than the width value. Estimate rwidth using width from peak width in the
row (x) direction. Note: width = 0 corresponds to no polynomial bgr""",

'bgr row power':"""Power of polynomial used in row (x) direction bgr fit.
Use pow < 0.5 for polynomials flatter than a circle (pow = 0 are linear)
and pow  > 0.5 for steeper.""",

'bgr row tan':"""Flag for use of tangents in row direction (x) bgr determination.
If 'True' then the local slope is removed when fitting a polynomial to a point.
Helpful for data sitting on a broad sloping background""",

'beam_slits':"""Enter the incident beam slit settings: beam_slits = {'horz':.6,'vert':.8}.
horz = beam horz width (at the sample) in mm (total width in lab-z / horizontal scattering plane)
vert = beam vert hieght (at the sample) in mm (total width in lab-x / vertical scattering plane).
If beam slits are 'None' or {} no area correction will be done""",

'det_slits':"""Enter the detector slit settings: det_slits = {'horz':1.,'vert':1.}.
horz = det horz width in mm (total width in lab-z / horizontal scattering plane)
vert = det vert hieght in mm (total width in lab-x / vertical scattering plane).
If detector slits are 'None' or {} only a spil-off correction will be computed""",

'geom':"""Enter goniometer geometry.  Options: 'psic'""",

'sample dia':"""Enter the diameter (in mm) of a round sample mounted on center
of the goniometer.  If this is 'None' then use the sample polygon or
no sample description. """,

'sample polygon':"""A list of lists that describe a general polygon sample shape, e.g.:    
    polygon =  [[1.,1.], [.5,1.5], [-1.,1.],[-1.,-1.],[0.,.5],[1.,-1.]]
Each entry is an [x,y] or [x,y,z] vector pointing to an apex of the sample
polygon. These vectors should be given in general lab frame coordinates
(x is lab frame vertical, y is positive along the direction of the beam,
the origin is the rotation center). The vectors need to be specified at a 
given set of angles (sample angles).

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
being the +y direction.  """,

'sample angles':"""Sample angles used for description of the sample polygon. Use
the format" angles = {'phi':123.5,'chi':flatchi:0.3,'eta':0.,'mu':0.}""",

'scale':"""Enter scale factor.  The scale factor multiplies by all the intensity
values. e.g.if Io ~ 1million cps then using 1e6 as the scale makes the normalized
intensity close to cps.  ie y = scale*y/norm"""
}

###############################################################################

class wxCtrData(model.Background, wxUtil):

    ###########################################################
    # Init and util methods
    ###########################################################
    
    def on_initialize(self, event):
        """
        Initialization

        including sizer setup, do it here
        self.setupSizers()
        """
        self.startup  = True
        self.dir      = '.'
        self.set = []
        self.slider_cnt = 0
        
        # set up shell
        self.shell = None
        self.init_shell()

        # Make sure Scan Data is loaded
        self.exec_line("import scandata")

        # init all the menus
        self.init_shell_items()
        self.init_gui()

    def init_shell_items(self,):
        self.init_CtrDataVar()
        self.init_AppendScanVar()

    def on_Update_mouseClick(self,event):
        self.init_shell_items()
    
    ###########################################################
    #   Menu
    ###########################################################

    def on_menuFileExit_select(self,event): 
        self.close()

    def on_menuFileReadCtr_select(self,event):        
        cdir = self.eval_line("pwd()")
        result = dialog.fileDialog(self, 'Open', cdir, '',"*")
        if result.accepted:
            path        = result.paths[0]
            dir,fname   = os.path.split(path)
            if os.path.abspath(dir) != os.path.abspath(cdir):
                dir = dir.replace("\\","\\\\")
                line = "cd('%s')" % dir
                self.eval_line(line)
        else:
            self.post_message("File selection cancelled.")
            return
        #print dir, fname
        self.save_file = path
        line = "restore %s" % path
        print line
        self.exec_line(line)
        print "Select ctr data variable"
        self.init_shell_items()

    def on_menuFileSaveCtr_select(self,event):
        cdir = self.eval_line("pwd()")
        result = dialog.fileDialog(self, 'Save', cdir, '',"*")
        if result.accepted:
            path        = result.paths[0]
            dir,fname   = os.path.split(path)
            if os.path.abspath(dir) != os.path.abspath(cdir):
                dir = dir.replace("\\","\\\\")
                line = "cd('%s')" % dir
                self.eval_line(line)
        else:
            self.post_message("File selection cancelled.")
            return
        #print dir, fname
        ctr = self.get_ctr_name()
        line = "save %s %s" % (ctr,fname) 
        print line
        self.exec_line(line)

    def on_menuFileWriteHKL_select(self,event):
        cdir = self.eval_line("pwd()")
        result = dialog.fileDialog(self, 'Save', cdir, '',"*")
        if result.accepted:
            path        = result.paths[0]
            dir,fname   = os.path.split(path)
            if os.path.abspath(dir) != os.path.abspath(cdir):
                dir = dir.replace("\\","\\\\")
                line = "cd('%s')" % dir
                self.eval_line(line)
        else:
            self.post_message("File selection cancelled.")
            return
        #print dir, fname
        ctr = self.get_ctr()
        if ctr != None:
            ctr.write_HKL(path)  

    def on_menuHelpDocumentation_select(self,event):
        self.exec_line("web 'http://cars9.uchicago.edu/ifeffit/tdl/Docs/Pds/CtrGui'")

    ###########################################################
    #  Ctr Data
    ###########################################################

    def init_CtrDataVar(self):
        """
        Initialize the list of variables. 
        """
        ctr_var = self.components.CtrDataVar.stringSelection
        tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        #tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False,instance=scandata.CtrData)
        tmp = tmp['var'] + tmp['ins']
        tmp.sort()
        self.components.CtrDataVar.items = tmp
        if ctr_var in tmp:
            self.components.CtrDataVar.stringSelection = ctr_var
        else:
            self.components.CtrDataVar.text = ''
        return

    def on_CtrDataVar_select(self,event):
        "select a Ctr data name and check it out"
        self.init_gui()
        if self.components.AutoPlotCtr.checked==True:
            self._plot_ctr()
        return

    def on_CtrDataVar_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init_gui()
            if self.components.AutoPlotCtr.checked==True:
                self._plot_ctr()
        else:
            event.skip()
        return

    def get_ctr_name(self,):
        if len(self.components.CtrDataVar.stringSelection) > 0:
            self.components.CtrDataVar.text = self.components.CtrDataVar.stringSelection
        ctr = self.components.CtrDataVar.text
        if len(ctr.strip()) == 0: return None
        name = "%s" % ctr.strip()
        return name

    def get_ctr(self):
        name = self.get_ctr_name()
        if name == None:
            print "Please provide a Ctr Data instance name"
            return None
        ctr = self.get_data(name)
        return ctr

    def set_ctr(self,ctr):
        name = self.get_ctr_name()
        return self.set_data(name,ctr)

    def check_ctr(self,):
        try:
            name = self.get_ctr_name()
            ctr  = self.get_data(name)
            #if isinstance(ctr,scandata.CtrData):
            if hasattr(ctr,'L'):
                self.post_message("Valid ctr object")
                return True
            else:
                self.post_message("Invalid ctr object")
                return False
        except:
            self.post_message("Invalid ctr object")
            return False

    ###########################################################
    #  Scan Data
    ###########################################################

    def init_AppendScanVar(self):
        """
        Initialize the list of variables. 
        """
        scan_var = self.components.AppendScanVar.stringSelection
        tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        #tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False,
        #                                                  instance=scandata.ScanData)
        tmp = tmp['var'] + tmp['ins']
        tmp.sort()
        tmp.insert(0,'')
        self.components.AppendScanVar.items = tmp
        if scan_var in tmp:
            self.components.AppendScanVar.stringSelection = scan_var
        else:
            self.components.AppendScanVar.stringSelection = ''
        return

    def on_AppendScanVar_select(self,event):
        """
        select a scan data name and check it out
        """
        s = self.get_scan()
        return

    def get_scan(self):
        scan = self.components.AppendScanVar.stringSelection
        if len(scan.strip()) == 0: return None
        name = "%s" % scan.strip()
        try:
            scan = self.get_data(name)
            #if isinstance(ctr,scandata.ScanData):
            if type(scan) != types.ListType:
                scan = [scan]
            if hasattr(scan[0],'get_positioner'):
                self.post_message("Valid scan object")
                return scan
            else:
                self.post_message("Invalid scan object")
                self.components.AppendScanVar.stringSelection = ''
                return None
        except:
            self.post_message("Invalid scan object")
            self.components.AppendScanVar.stringSelection = ''
            return None

    def on_AppendScan_mouseClick(self,event):
        """
        Select a scan data name and check it out
        """
        ctr = self.get_ctr()
        if ctr == None: return
        npts = len(ctr.L)
        scans = self.get_scan()
        if scans == None: return
        ctr.append_scans(scans)
        set  = range(npts,len(ctr.L))
        self.set = set
        self.update_gui_from_ctr()
        self.process_set()
        self.components.PointNum.stringSelection = str(npts)
        #
        self.components.AppendScanVar.stringSelection = ''
        return

    ###########################################################
    #  Point Nums
    ###########################################################
    
    def init_PointNum(self,pnts=None):
        """
        Initialize the point nums. 
        """
        if pnts == None:
            self.components.PointNum.items = ['0']
            self.components.PointNum.stringSelection = '0'
            self.set = []
        else:
            if len(self.components.PointNum.stringSelection) > 0:
                self.components.PointNum.text = self.components.PointNum.stringSelection
            try:
                point  = int(self.components.PointNum.text)
            except:
                print "Enter valid point"
                return
            if point not in pnts: point = 0
            self.components.PointNum.items = map(str,pnts)
            self.components.PointNum.stringSelection = str(point)
            self.components.PointNum.text = str(point)
    
    def init_AnchorPointNum(self,pnts=None):
        """
        Initialize the point nums. 
        """
        if pnts == None:
            self.components.AnchorPointNum.items = ['','0']
            self.components.AnchorPointNum.stringSelection = ''
        else:
            point= self.components.AnchorPointNum.stringSelection
            pnts = map(str,pnts)
            pnts.insert(0,'')
            if point not in pnts: point = ''
            self.components.AnchorPointNum.items = pnts
            self.components.AnchorPointNum.stringSelection = point

    def on_PointNum_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN or keyCode == 372:
            self.update_point()
        else:
            event.skip()
        return
    
    def on_PointNum_select(self,event):
        self.update_point()
        return

    def on_NextPoint_mouseClick(self,event):
        try:
            point  = int(self.components.PointNum.stringSelection)
        except:
            print "Enter valid point"
            return
        point = str(point + 1)
        if point in self.components.PointNum.items:
            self.components.PointNum.stringSelection = point
            self.update_point()
        return

    def on_PrevPoint_mouseClick(self,event):
        try:
            point  = int(self.components.PointNum.stringSelection)
        except:
            print "Enter valid point"
            return
        point = str(point - 1)
        if point in self.components.PointNum.items:
            self.components.PointNum.stringSelection = point
            self.update_point()
        return

    def on_PointNum_select(self,event):
        self.update_point()
        return

    def on_SetPointClick_mouseClick(self,event):
        ctr = self.get_ctr()
        if ctr == None: return
        point = ctr.get_idx()
        try:
            point = str(point[0])
            if point in self.components.PointNum.items:
                self.components.PointNum.stringSelection = point
                self.update_point()
        except:
            pass
        return

    def on_SetPointsClick_mouseClick(self,event):
        ctr = self.get_ctr()
        if ctr == None: return
        points = ctr.get_points()
        if points!=None:
            self.set = points
            self.update_gui_from_ctr()
            self.process_set()
        return

    def on_PointUpdate_mouseClick(self,event):
        self.update_point()

    def on_SetUpdate_mouseClick(self,event):
        self.process_set()

    def on_ToggleBad_mouseClick(self,event):
        ctr = self.get_ctr()
        if ctr == None: return
        point  = int(self.components.PointNum.stringSelection)
        if point in ctr.bad:
            ctr.bad.remove(point)
        else:
            ctr.bad.append(point)
        self.update_point()
    
    def update_point(self,update_gui=True):
        ctr = self.get_ctr()
        if ctr == None: return
        apnt = self.components.AnchorPointNum.stringSelection
        if len(apnt) > 0:
            apnt = int(apnt)
            (intpar,corrpar) = ctr_data.get_params(ctr,apnt)
            if self.components.SetCorrParams.checked == False:
                corrpar = {}
            if self.components.SetIntParams.checked == False:
                intpar = {}
            point = int(self.components.PointNum.stringSelection)
            ctr_data.set_params(ctr,point,intpar=intpar,corrpar=corrpar)
        if self.components.AutoIntegrate.checked == True:
            self.integrate_point(update_plots=update_gui)
        if update_gui:
            self.update_gui_from_ctr()
            if self.components.AutoPlotScan.checked == True:
                self._plot_scan()
            if self.components.AutoPlotCorr.checked == True:
                self._plot_corr()

    def process_set(self):
        for point in self.set:
            self.components.PointNum.stringSelection = str(point)
            print "Process point ", point
            self.update_point(update_gui=False)
        self.components.PointNum.stringSelection = str(self.set[0])
        self.update_gui_from_ctr()
        if self.components.AutoPlotScan.checked == True:
            self._plot_scan()

    ###########################################################
    #  Parameters
    ###########################################################
    
    def init_IntParamList(self,params=None):
        if params == None:
            self.components.IntParamList.items = []
        else:
            items = []
            for (key,val) in params.items():
                items.append([key,str(val)])
            items.sort()
            self.components.IntParamList.items = items
    
    def init_CorrParamList(self,params=None):
        if params == None:
            self.components.CorrParamList.items = []
        else:
            items = []
            for (key,val) in params.items():
                items.append([key,str(val)])
            items.sort()
            self.components.CorrParamList.items = items
        
    def on_IntParamList_select(self,event):
        self.components.ParamDescr.text = ''
        selected = self.components.IntParamList.getStringSelection()
        (name,val) = selected[0]
        self.components.ParamName.text = name
        self.components.ParamVal.text = val
        #
        ctr = self.get_ctr()
        point = int(self.components.PointNum.stringSelection)
        (s,spnt) = ctr.get_scan(point)
        lbls = s.scalers.keys()
        if hasattr(s,'image'):
            lbls.extend(s.image.peaks.keys())
        lbls.sort()
        #
        if name == 'I':
            t = PARAM_DESCR['I'] % lbls
            self.components.ParamDescr.text = t
        elif name == 'Inorm':
            t = PARAM_DESCR['Inorm'] % lbls
            self.components.ParamDescr.text = t
        elif name == 'Ierr':
            t = PARAM_DESCR['Ierr'] % lbls
            self.components.ParamDescr.text = t
        elif name == 'Ibgr':
            t = PARAM_DESCR['Ibgr'] % lbls
            self.components.ParamDescr.text = t
        elif name == 'image roi':
            t = PARAM_DESCR['image roi']
            self.components.ParamDescr.text = t
        elif name == 'image rotangle':
            t = PARAM_DESCR['image rotangle']
            self.components.ParamDescr.text = t
        elif name == 'bgr flag':
            t = PARAM_DESCR['bgr flag']
            self.components.ParamDescr.text = t
        elif name == 'bgr col nbgr':
            t = PARAM_DESCR['bgr col nbgr']
            self.components.ParamDescr.text = t
        elif name == 'bgr col width':
            t = PARAM_DESCR['bgr col width']
            self.components.ParamDescr.text = t
        elif name == 'bgr col power':
            t = PARAM_DESCR['bgr col power']
            self.components.ParamDescr.text = t
        elif name == 'bgr col tan':
            t = PARAM_DESCR['bgr col tan']
            self.components.ParamDescr.text = t
        elif name == 'bgr row nbgr':
            t = PARAM_DESCR['bgr row nbgr']
            self.components.ParamDescr.text = t
        elif name == 'bgr row width':
            t = PARAM_DESCR['bgr row width']
            self.components.ParamDescr.text = t
        elif name == 'bgr row power':
            t = PARAM_DESCR['bgr row power']
            self.components.ParamDescr.text = t
        elif name == 'bgr row tan':
            t = PARAM_DESCR['bgr row tan']
            self.components.ParamDescr.text = t
        return

    def on_CorrParamList_select(self,event):
        self.components.ParamDescr.text = ''
        selected = self.components.CorrParamList.getStringSelection()
        (name,val) = selected[0]
        self.components.ParamName.text = name
        self.components.ParamVal.text = val
        #
        ctr = self.get_ctr()
        point = int(self.components.PointNum.stringSelection)
        (s,spnt) = ctr.get_scan(point)
        #
        if name == 'beam_slits':
            t = PARAM_DESCR['beam_slits']
            self.components.ParamDescr.text = t
        elif name == 'det_slits':
            t = PARAM_DESCR['det_slits']
            self.components.ParamDescr.text = t
        elif name == 'geom':
            t = PARAM_DESCR['geom']
            self.components.ParamDescr.text = t
        elif name == 'sample dia':
            t = PARAM_DESCR['sample dia']
            self.components.ParamDescr.text = t
        elif name == 'sample polygon':
            t = PARAM_DESCR['sample polygon']
            self.components.ParamDescr.text = t
        elif name == 'sample angles':
            t = PARAM_DESCR['sample angles']
            self.components.ParamDescr.text = t
        elif name == 'scale':
            t = PARAM_DESCR['scale']
            self.components.ParamDescr.text = t
        return
    
    def on_ParamVal_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN or keyCode == 372:
            self.update_params()
            self.update_ctr_from_params()
            if self.components.AutoIntegrate.checked == True:
                self.integrate_point()
            self.update_gui_from_ctr()
        else:
            event.skip()
        return

    def on_ParamSet_mouseClick(self,event):
        name = self.components.ParamName.text
        if len(name) == 0: return
        if name == 'image roi':
            #print 'image roi'
            pyplot.figure(1)
            (x1,x2,y1,y2) = pyplot.axis()
            roi  = [int(x1),int(y1),int(x2),int(y2)]
            roi = image_data._sort_roi(roi)
            self.components.ParamVal.text = str(roi)
        #
        self.update_params()
        self.update_ctr_from_params()
        if self.components.AutoIntegrate.checked == True:
            self.integrate_point()
        self.update_gui_from_ctr()

    def on_CopyParamsToAll_mouseClick(self,event):
        ctr = self.get_ctr()
        if ctr == None: return
        npts = len(ctr.L)
        if npts == 0: return
        for j in range(npts):
            print "Setting params for point:", j
            self.update_ctr_from_params(point=j)
            self.integrate_point(point=j,update_plots=False)
        if self.components.AutoPlotCtr.checked==True:
            self._plot_ctr()
    
    def update_params(self):
        name = self.components.ParamName.text
        if len(name) == 0: return
        val = self.components.ParamVal.text
        val = val.strip()
        if val == '': return
        #
        int_items = self.components.IntParamList.items
        corr_items = self.components.CorrParamList.items
        idx1 = self._get_par_idx(name,int_items)
        idx2 = self._get_par_idx(name,corr_items)
        if idx1 > -1:
            int_items[idx1][1] = val
            self.components.IntParamList.items = int_items
        elif idx2 > -1:
            corr_items[idx2][1] = val
            self.components.CorrParamList.items = corr_items
        else:
            return

    def _get_par_idx(self,name,items):
        for j in range(len(items)):
            (n,v) = items[j]
            if n == name: return j
        return -1
    
    def update_ctr_from_params(self,point=None):
        ctr = self.get_ctr()
        if ctr == None: return
        intpar = {}
        for (n,v) in self.components.IntParamList.items:
            intpar[n] = str(v)
        corrpar = {}
        for (n,v) in self.components.CorrParamList.items:
            corrpar[n] = str(v)
        if point == None:
            point = int(self.components.PointNum.stringSelection)
        ctr_data.set_params(ctr,point,intpar=intpar,corrpar=corrpar)

    ######################################################
    # image max
    ######################################################

    def on_Imax_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN or keyCode == 372:
            im_max = self.components.Imax.text
            (m1,m2) = self.image_max()
            if im_max == '-1':
                self.components.ImaxSlider.value = int(m1)
                self.components.ImaxSlider.max = int(m1)
            elif int(im_max) > m1:
                self.components.ImaxSlider.max = int(im_max)
                self.components.ImaxSlider.value = int(im_max)
            #
            self.set_image_max()
            if self.components.AutoPlotScan.checked == True:
                self._plot_scan()
        else:
            event.skip()
        return

    def on_ImaxSlider_select(self,event):
        if self.slider_cnt < 2:
            self.slider_cnt  = self.slider_cnt  + 1
            return
        else:
            self.slider_cnt  = 1
        #
        self.components.Imax.text = str(self.components.ImaxSlider.value)
        self.set_image_max()
        if self.components.AutoPlotScan.checked == True:
            self._plot_scan()
        return

    def set_image_max(self):
        ctr = self.get_ctr()
        if ctr == None: return
        point = int(self.components.PointNum.stringSelection)
        (scan,spnt) = ctr.get_scan(point)
        immax = int(self.components.Imax.text)
        scan.image.im_max[spnt] = immax

    def image_max(self):
        ctr = self.get_ctr()
        if ctr == None: return
        point = int(self.components.PointNum.stringSelection)
        (scan,spnt) = ctr.get_scan(point)
        immax = num.max(scan.image.image[spnt])
        im_max = scan.image.im_max[spnt]
        return (immax,im_max)
    
    ######################################################
    # update stuff from ctr etc
    ######################################################
    
    def integrate_point(self,point=None,update_plots=True):
        ctr = self.get_ctr()
        if ctr== None: return
        if point == None:
            point = int(self.components.PointNum.text)
        if update_plots == True:
            if self.components.AutoPlotIntegration.checked == True:
                plot = True
                fig = 3
            else:
                plot = False
                fig = None
            ctr.integrate_point(point,plot=plot,fig=fig)
            if self.components.AutoPlotCtr.checked==True:
                self._plot_ctr()
        else:
            ctr.integrate_point(point)
    
    def init_gui(self):
        if self.startup:
            self.startup = False
            cmap = copy.copy(pyplot.cm.cmapnames)
            cmap.insert(0,'')
            self.components.ColorMap.items = cmap 
            self.components.ColorMap.stringSelection = ''
        check = self.check_ctr()
        if check == False:
            self.init_PointNum()
            self.init_AnchorPointNum()
            self.init_IntParamList()
            self.init_CorrParamList()
            self.components.ParamDescr.text = ''
            self.components.ParamVal.text = ''
            self.components.ParamName.text = ''
            self.components.Imax.text = ''
        else:
            self.update_gui_from_ctr()
        return

    def update_gui_from_ctr(self,):
        """
        update gui from a ctr data instance
        """
        ctr = self.get_ctr()
        if ctr == None: return
        npts  = len(ctr.L)
        pnts  = range(npts)
        self.init_PointNum(pnts=pnts)
        self.init_AnchorPointNum(pnts=pnts)
        #
        point = int(self.components.PointNum.text)
        set   = self.set
        if point in ctr.bad:
            bad_flag = True
        else:
            bad_flag = False
        header   = CTR_HEADER % (str(len(ctr.L)),str(set),str(point),
                                 ctr.H[point],ctr.K[point],ctr.L[point],
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
                                 ctr.I[point],ctr.Ierr[point],ctr.Ibgr[point],
                                 ctr.ctot[point],ctr.F[point],ctr.Ferr[point])
        self.components.PointData.text = header
        #
        (intpar,corrpar) = ctr_data.get_params(ctr,point)
        self.init_IntParamList(params=intpar)
        self.init_CorrParamList(params=corrpar)
        #
        if ctr.scan_type[point] == 'image':
            (immax,im_max) = self.image_max()
            self.components.ImaxSlider.min = 0
            self.components.ImaxSlider.max = int(immax)
            if im_max == -1:
                self.components.Imax.text = '-1'
                self.components.ImaxSlider.value = int(immax)
            else:
                self.components.Imax.text = str(im_max)
                self.components.ImaxSlider.value = int(im_max)
                
    ###########################################################
    #   Plot
    ###########################################################

    ######################################################
    def on_PlotCtr_mouseClick(self,event):
        self._plot_ctr()

    def _plot_ctr(self):
        ncol = self.components.PlotCol.text
        try:
            ncol = int(ncol)
        except:
            ncol = 2
        ctr = self.get_ctr()
        if ctr == None: return
        point = int(self.components.PointNum.text)
        if self.components.CtrPlotInt.checked == True:
            ctr.plot_I(fig=0,num_col=ncol,spnt=point)
        else:
            ctr.plot(fig=0,num_col=ncol,spnt=point)
    
    ######################################################
    def on_PlotScanData_mouseClick(self,event):
        self._plot_scan()
        
    def _plot_scan(self):
        ctr = self.get_ctr()
        if ctr == None: return
        pnt = int(self.components.PointNum.stringSelection)
        if ctr.scan_type[pnt] == 'image':
            if self.components.ColorMap.stringSelection.strip()!='':
                cmap = self.components.ColorMap.stringSelection.strip()
                cmap = str(cmap)
            else:
                cmap=None
            ctr.plot_point(idx=pnt,fig=1,cmap=cmap)
        else:
            ctr.plot_point(idx=pnt,fig=1)

    ######################################################
    def on_CorrPlot_mouseClick(self,event):
        self._plot_corr()
        
    def _plot_corr(self):        
        ctr = self.get_ctr()
        if ctr == None: return
        point = int(self.components.PointNum.stringSelection)
        (scan,spnt) = ctr.get_scan(point)
        corr_params = ctr.corr_params[point]
        corr = ctr_data._get_corr(scan,spnt,corr_params)
        if ctr.scan_type[point] == 'image':
            corr.ctot_stationary(plot=True,fig=2)

##################################################
if __name__ == '__main__':
    app = model.Application(wxCtrData)
    app.MainLoop()
