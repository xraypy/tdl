"""
Menu function to handle interactive rasd data analyses

Authors/Modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

Todo:
Create RasdList not from *.rsd files but directly from RasdData objects
Functions to set Emin, Emax for individual scans and to remove bad datapoints
Functions to append/remove scans to an existing RasdList object
Least squares fitting of R, theta, and DW.
anisotropic DW-Factors
"""
#######################################################################
import numpy as num
from matplotlib import pyplot
import os.path
from   pds.lib.shellutil import Menu 
from   pds.lib.shellutil import get_tf, get_yn, get_int, get_flt, get_str, get_flt_list
import rasd_ana
from atomic import f0data as database

########################################################################
RASD_HEADER = """
##########################################
Number of rasd scans   = %s
Current scan           = %s
Current scan AR        = %s
Current scan PR        = %s
Current scan E0shift   = %s
"""

RASD_LABELS = ['plot_norm','setE0shift','setPRstart','useInFourier','useInRefine','FourierParams',
              'runFourier','RefinementParams','runRefine','select','next','previous','done']
RASD_DESCR = ["Plot scan with (default) or without normalization",
             "Set E0 shift for this scan",
             "set a starting value for the resonant Phase (PR)",
             "use AR and PR of this scan in Fourier synthesis (True/False)",
             "use this scan in the structure refinement (True/False)",
             "Set Fourier synthesis parameters",
             "Run Fourier synthesis",
             "Set simmulated annealing parameters",
             "Run simmulated annealing structure refinement",
             "Select scan",
             "Select next scan ",
             "Select previous scan", 
             "All Done"]

########################################################################
def rasd_menu(rasddata, cell = None, bulk_file = None, sur_file = None, f1f2_file = None, E0 = 0.):
    """
    Interactively create RasdList object and inspect/analyze data in it 
    """
    if cell == None:
        ok = False
        while not ok:
            cell = get_flt_list(prompt='Enter unit cell [a,b,c,alpha,beta,gamma,delta1,delta2]')
            if len(cell) != 8:
                print 'invalid unit cell'
            else: ok = True
        
    if bulk_file == None:
        ok = False
        while not ok:
            bulk_file = get_str(prompt='Enter path and\or name of your bulk file (ROD .bul file without header and cell definition)')
            if not os.path.exists(bulk_file):
                print bulk_file+' not found'
            else:
                ok = True
    bulk = []
    f = file(bulk_file,'r')
    data = f.readlines()
    f.close()
    for i in range(len(data)):
        if '#' not in data[i]:
            tmp = str.rsplit(data[i])
            bulk.append([tmp[0],float(tmp[1]),float(tmp[2]),float(tmp[3]),float(tmp[4])])
                
    if sur_file == None:
        ok = False
        while not ok:
            sur_file = get_str(prompt='Enter path and\or name of your surface file (ROD .sur file without header and cell definition)')
            if not os.path.exists(sur_file):
                print sur_file+' not found'
            else:
                ok = True
    surface = []
    f = file(sur_file,'r')
    data = f.readlines()
    f.close()
    for i in range(len(data)):
        if '#' not in data[i]:
            tmp = str.rsplit(data[i])
            surface.append([tmp[0],float(tmp[1]),float(tmp[2]),float(tmp[3]),float(tmp[4]),float(tmp[5])])

    if f1f2_file == None:
        ok = False
        while not ok:
            f1f2_file = get_str(prompt='Enter path and\or name of your f1f2 file (HEPHAESTUS or experimental f1f2 file)')
            if not os.path.exists(f1f2_file):
                print f1f2_file+' not found'
            else:
                ok = True
    f = file(f1f2_file, 'r')
    data = f.readlines()
    f.close()
    f1f2 = num.ndarray((0,3),float)
    for i in range(len(data)):
        if '#' not in data[i]:
            tmp = str.rsplit(data[i])
            f1f2 = num.append(f1f2, [[int(float(tmp[0])),float(tmp[1]),float(tmp[2])]], axis = 0)
    f1f2 = f1f2.transpose()

    if E0 == 0:
        E0 = get_flt(prompt='Enter E0 in eV')

    allrasd = rasd_ana.read_RSD(cell, bulk, surface, database, rasddata, f1f2, E0)
        
    prompt   = 'Select option >'
    ret      = ''

    # make menu
    m = Menu(labels=RASD_LABELS,
                 descr=RASD_DESCR,
                 sort=False, matchidx=True)
    scan_pnt = 0
    norm = True
    
    # loop
    while ret != 'done':
        allrasd.list[scan_pnt] = rasd_ana.RASD_Fourier(allrasd,scan_pnt)
        header   = RASD_HEADER % (str(allrasd.dims),str(scan_pnt),
                                 str(allrasd.list[scan_pnt].AR),str(allrasd.list[scan_pnt].PR),
                                 str(allrasd.list[scan_pnt].e0shift))
        m.header = header
        allrasd.list[scan_pnt].plot(norm = norm, fig = 1)
        ret      = m.prompt(prompt)
        if ret == 'plot_norm':
            norm = get_tf(prompt = 'Plot normalized or not (True/False)', default = norm)
        elif ret == 'setE0shift':
            allrasd.list[scan_pnt].e0shift = get_flt(prompt = 'Enter e0 shift for this scan',
                                                       default = allrasd.list[scan_pnt].e0shift)
            allrasd.list[scan_pnt].E = allrasd.list[scan_pnt].Eorig + allrasd.list[scan_pnt].e0shift
            allrasd.list[scan_pnt].E0 = allrasd.E0 + allrasd.list[scan_pnt].e0shift
        elif ret == 'setPRstart':
            allrasd.list[scan_pnt].PR = get_flt(prompt = 'Enter PR start (0-1)', default = allrasd.list[scan_pnt].PR,
                                                min = 0, max = 1)
        elif ret == 'useInFourier':
            allrasd.list[scan_pnt].use_in_Fourier = get_tf(prompt= 'Use this scan in Fourier synthesis (True/False)?',
                                                           default = allrasd.list[scan_pnt].use_in_Fourier)
        elif ret == 'useInRefine':
            allrasd.list[scan_pnt].use_in_Refine = get_tf(prompt= 'Use this scan for Structure Refinement (True/False)?',
                                                           default = allrasd.list[scan_pnt].use_in_Refine)
        elif ret == 'FourierParams':
            allrasd.ZR = get_int(prompt= 'Atomic number of resonant element', default = allrasd.ZR)
            allrasd.xf = get_flt(prompt= 'number of unit cells along a', default = allrasd.xf)
            allrasd.yf = get_flt(prompt= 'number of unit cells along b', default = allrasd.yf)
            allrasd.zf = get_flt(prompt= 'number of unit cells along c', default = allrasd.zf)
            allrasd.an = get_int(prompt= 'number of datapoints along a', default = allrasd.an)
            allrasd.bn = get_int(prompt= 'number of datapoints along b', default = allrasd.bn)
            allrasd.cn = get_int(prompt= 'number of datapoints along c', default = allrasd.cn)
            allrasd.Plusminus = get_tf(prompt= 'calculate rho(e-) above and below surface', default = allrasd.Plusminus)
        elif ret == 'runFourier':
            allrasd.Fourier = []
            data = []
            for rasd in allrasd.list:
                if rasd.use_in_Fourier:
                    allrasd.Fourier.append([rasd.Q[0],rasd.Q[1],rasd.Q[2],rasd.AR,rasd.PR])
                    data.append(rasd.file)
            if allrasd.Fourier == []: print 'No scans specified for use in Fourier'
            else:
                rasd_ana.Fourier_synthesis(allrasd.Fourier, allrasd.cell, allrasd.ZR,allrasd.xf,allrasd.yf,allrasd.zf,
                                             allrasd.an,allrasd.bn,allrasd.cn,allrasd.Plusminus)
                print 'Scans used for Fourier synthesis: \n'+str(data) 
        elif ret == 'RefinementParams':
            allrasd.natoms = get_int(prompt = 'How many atoms in structure model?', default = allrasd.natoms)
            if allrasd.Rmin == None or len(allrasd.DWmin) != allrasd.natoms:
                allrasd.Rmin = num.zeros((allrasd.natoms,3),float)
                allrasd.Rmax = num.ndarray((0,3),float)
                for i in range(allrasd.natoms):
                    allrasd.Rmax = num.append(allrasd.Rmax, [[allrasd.cell[0],allrasd.cell[1],allrasd.cell[2]]],axis = 0)
                allrasd.thetamin = num.ones((allrasd.natoms),float) *0.05
                allrasd.thetamax = num.ones((allrasd.natoms),float) *0.999
                allrasd.DWmin = num.ones((allrasd.natoms),float) * 12
                allrasd.DWmax = num.ones((allrasd.natoms),float) * 90
            for i in range(allrasd.natoms):
                allrasd.Rmin[i] = get_flt_list(prompt = ('Rmin of atom '+str(i+1)+' [xmin,ymin,zmin] (Angstroem)'), default = allrasd.Rmin[i])
            for i in range(allrasd.natoms):
                allrasd.Rmax[i] = get_flt_list(prompt = ('Rmax of atom '+str(i+1)+' [xmax,ymax,zmax] (Angstroem)'), default = allrasd.Rmax[i])
            for i in range(allrasd.natoms):
                allrasd.thetamin[i] = get_flt(prompt = ('min. occupancy of atom '+str(i+1)), default = allrasd.thetamin[i])
            for i in range(allrasd.natoms):
                allrasd.thetamax[i] = get_flt(prompt = ('max. occupancy of atom '+str(i+1)), default = allrasd.thetamax[i])
            for i in range(allrasd.natoms):
                allrasd.DWmin[i] = get_flt(prompt = ('min. DW-factor for atom '+str(i+1)), default = allrasd.DWmin[i])
            for i in range(allrasd.natoms):
                allrasd.DWmax[i] = get_flt(prompt = ('max. DW-factor for atom '+str(i+1)), default = allrasd.DWmax[i])
            allrasd.Tstart = get_flt(prompt = 'Startimg Temperature for simmulated annealing', default = allrasd.Tstart)
            allrasd.Tend = get_flt(prompt = 'End Temperature for simmulated annealing', default = allrasd.Tend)
            allrasd.cool = get_flt(prompt = 'cooling factor for simmulated annealing (usually 0.7 - 0.95)', default = allrasd.cool, min = 0.1, max = 0.99)
            allrasd.maxrun = get_flt(prompt = 'max. repeats at one Temp.', default = allrasd.maxrun)
            allrasd.MC = get_flt(prompt = 'max fractional parameter change at T = 100', default = allrasd.MC)
            allrasd.RMS_count_max = get_flt(prompt = 'max. annealing (T) steps without improvement', default = allrasd.RMS_count_max)
            allrasd.factor = get_flt(prompt = 'some factor (Boltz. = exp(d(RMS)*factor/T)) in decision process (~1e4 - 1e7)', default = allrasd.factor)
        elif ret == 'runRefine':
            if allrasd.Rmin == None: print 'Please define simmulated annealing parameters first'
            else:
                allrasd.reflist = []
                for rasd in allrasd.list:
                    if rasd.use_in_Refine:
                        allrasd.reflist.append(rasd)
                if allrasd.reflist == []: print 'No scans specified for use in refinement'
                else: rasd_ana.simulated_annealing(allrasd)
        elif ret == 'select':
            scan_pnt = get_int(prompt='Enter scan number',
                               default=scan_pnt,min=0,max = allrasd.dims-1)
        elif ret == 'next':
            if scan_pnt + 1 < allrasd.dims: 
                scan_pnt = scan_pnt + 1
        elif ret == 'previous':
            if scan_pnt - 1 > -1: 
                scan_pnt = scan_pnt - 1
        else:
            pass

########################################################################
