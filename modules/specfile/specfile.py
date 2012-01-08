"""
Class to handle spec files

Authors/Modifications:
----------------------
* Matt Newville (newville@cars.uchicago.edu)
* Tom Trainor (tptrainor@alaska.edu)

Notes on the spec format:
-------------------------
An example scan would look like the following:

{{{
#S 16  hklscan  -0.000113467 -0.000113467  0.000280176 0.000280176  6.06 6.56  9 300
#D Mon Feb 23 07:16:59 2009
#T 300  (Seconds)
#G0 0 0 1 0.003991574459 0.0007565094145 1 0 0 0 0 0 0 50 0 0 1 4 4 5 4 0 0
#G1 8.094 4.988 6.071 90 90 90 0.7762769097 1.259660246 1.034950635 90 90 90 0 0 4 2 0 2.481 -0.0009 0.0008 -0.1244 175.2192 31.965 16.27375 15.0903 7.5478 11.9002 -96.0482 17.59425 0.35625 0.835801 0.835801
#G3 -0.2392618672 1.198330644 -0.002648198747 -0.7384726 -0.3882605276 -0.005057797345 -0.004221190899 0.001168841303 1.034934888
#G4 -0.0001134668114 0.0002801762336 5.940014158 0.835801 23.95543252 24.31406748 0.2847499978 48.2695 0.07510498876 0.1793176371 -0.0004019917579 -0.0001447829911 0.4831 108.0007 2 0 0 0 0 12 0 0 2.0803 123.1461 0 0 0 0 -180 -180 -180 -180 -180 -180 -180 -180 -180 0 0 0 0 0 0 0 0 0 0 0 0
#Q -0.000113467 0.000280176 5.94001
#P0 -0.0003 0.0002 -0.1344 178.1354 48.2695 24.4195 -0.0562 -0.1754
#P1 178.079 -0.3000524 -0.60010481 2.3800315 13.979 -5350 -800 -3079.1521
#P2 -4852.8589 0 0 0 0 -1868.5726 -270 -989.82892
#P3 -300 0.031894205 8000.0661 -0.0019447686 8000.0241 
#ATTEN 0 0 0 0
#ALP_BET 23.9554 24.3141
#ENERGY 14834.2
#N 15
#L H  K  L  del  eta  chi  phi  nu  mu  Alpha  Beta  Epoch  io  i1  Bicron  AmpTek_sc  ROI1  ROI2  ROI3  Seconds  IROI
-0.000113467 0.000280176 6.06     0.0007    -0.0001    -0.1338   178.0821    49.3115    24.9395  24.475   24.836  23051 9.4473e+06 51761 0 0 0 0 0 300.004 1.33933e+06
-0.000113467 0.000280176 6.11556     0.0009    -0.0004    -0.1334   178.0920    49.7983    25.1803  24.716   25.082  23352 9.4516e+06 51080 0 0 0 0 0 300.004 425333
-0.000113467 0.000280176 6.17111     0.0009    -0.0002    -0.1339   178.1012    50.2840    25.4237   24.96   25.324  23653 9.45093e+06 51891 0 0 0 0 0 300.004 236574
-0.000113467 0.000280176 6.22667     0.0012     0.0000    -0.1339   178.1102    50.7698    25.6640    25.2    25.57  23954 9.47686e+06 51705 0 0 0 0 0 300.004 168015
-0.000113467 0.000280176 6.28222     0.0011     0.0002    -0.1339   178.1192    51.2522    25.9070  25.443   25.809  24255 9.47265e+06 52003 0 0 0 0 0 300.004 135775
-0.000113467 0.000280176 6.33778     0.0012     0.0000    -0.1341   178.1282    51.7398    26.1525  25.688   26.051  24556 9.43717e+06 51906 0 0 0 0 0 300.004 115198
-0.000113467 0.000280176 6.39333     0.0011     0.0002    -0.1339   178.1372    52.2278    26.3938   25.93   26.298  24857 9.47247e+06 52130 0 0 0 0 0 300.004 105605
-0.000113467 0.000280176 6.44889     0.0016     0.0001    -0.1342   178.1463    52.7188    26.6382  26.174   26.545  25158 9.45268e+06 52698 0 0 0 0 0 300.004 95127
-0.000113467 0.000280176 6.50444     0.0009     0.0001    -0.1342   178.1543    53.2058    26.8848  26.421   26.785  25459 5.48945e+06 52787 0 0 0 0 0 300.004 51738
-0.000113467 0.000280176 6.56     0.0011     0.0004    -0.1344   178.1634    53.6970    27.1283  26.664   27.033  25760 56408 52766 0 0 0 0 0 300.004 29
}}}

* #S indicates a new scan and gives scan number, and scan command
* #D gives time and data
* #T gives count time
* #G0,#G1,#G3,#G4 give geometry data and are described below
* #Q gives the HKL value at the start of the scan
* #P0,#P1,#P2,#P3 give the motor settings at the start of the scan
* Other tags such as ATTEN, ENERGY etc may be added by the user
* #L gives labels of the variables recorded during the scan
* Then the scan data

The P arrays have corresponding motor labels.  Somewhere in the file
(at the top for sure, then later if they ever get redefined) there
should be somthing like the following in the file:

{{{
#O0 TwoTheta     theta       chi       phi        Nu       Psi     Omega     Kappa
#O1      Phi   SampleX   SampleY   SampleZ  Mono_Ang     S1_Hp     S1_Hw     S1_Vp
#O2    S1_Vw     S2_Hp     S2_Hw     S2_Vp     S2_Vw     S3_Hp     S3_Hw     S3_Vp
#O3    S3_Vw     S4_Hp     S4_Hw     S4_Vp     S4_Vw
}}}

Since there should always be a 1:1 correspondence between these motor names
and the P values, the below routine parses the names and values into a
dictionary.  

The G array returned by the below routing is a concatenation of
all the #G lines.  The contents and meaning of the G array 
likely depends on which geometry you are using.  At the bottom
of this module are some additional notes on parsing the G array for the
psic geometry

"""
#######################################################################

import numpy as num
import os
import types

#######################################################################
class SpecFile:
    """
    A spec file
    """
    def __init__(self, fname):
        """
        Initialize

        Parameters:
        -----------
        * fname is the specfile name (including full path)
        """
        self.path, self.fname = os.path.split(fname)
        self.max_scan = 0
        self.min_scan = 0
        self._mtime    = 0
        self._lines    = []
        self._summary  = []
        self._ok       = False
        self.read()

    def __repr__(self):
        """ display """
        self.read()
        lout = "Spec file: %s" % self.fname
        lout = "%s\nPath: %s" % (lout, os.path.join(self.path))
        lout = "%s\nFirst scan number: %i" % (lout,self.min_scan)
        lout = "%s\nLast scan number:  %i" % (lout,self.max_scan)
        lout = "%s\nLast scan: %s" % (lout, self._summary[self.max_scan-1]['date'])
        return lout

    def read(self):
        """
        Read the specfile

        This will re-read the file if its time stamp
        has changed since the last read
        """
        try:
            fname = os.path.join(self.path, self.fname)
            if os.path.getmtime(fname) != self._mtime:
                #print "Reading spec file %s" % fname
                #f  = Util.file_open(fname,default_path=self.path)
                f  = open(fname)
                self._mtime = os.path.getmtime(fname)
                self._lines = f.readlines()
                f.close()
                self._summarize()
                self._ok = True
        except IOError:
            print  '**Error reading file ', fname
            self._ok = False

    def _summarize(self):
        """
        summarize
        """
        lineno = 0
        (mnames,cmnd,date,xtime,Gvals,q,Pvals,atten,energy,lab,lStart,lStop,aborted) = (None,None,None,None,None,None,None,None,None,None,None,None,False)
        (index, ncols, n_sline, dataStart, dataStop) = (0,0,0,0,0)
        self._lines.append('\n')
        allLines = iter(self._lines)
        for i in allLines:
            lineno = lineno + 1
            i  = i[:-1]
            # get motor names: they should be at the top of the file
            # but they can be reset anywhere in the file
            if (i[0:2] == '#O'):
                if i[2] == '0': mnames = ''
                mnames = mnames + i[3:]
            # get scan number
            elif (i[0:3] == '#S '):
                v     = i[3:].split()
                index = int(v[0])
                cmnd  = i[4+len(v[0]):]
                n_sline= lineno
            elif (i[0:3] == '#D '):
                date = i[3:]
            elif (i[0:3] == '#T '):
                xtime = i[3:]
            elif (i[0:2] == '#G'):
                if i[2] == '0': Gvals = ''
                Gvals = Gvals + i[3:]
            elif (i[0:3] == '#Q '):
                q = i[3:]
            elif (i[0:2] == '#P'):
                if i[2] == '0': Pvals = ''
                Pvals = Pvals + i[3:]
            elif (i[0:3] == '#N '):
                ncols = int(i[3:])
            elif (i[0:3] == '#AT'):
                atten = i[6:]
            elif (i[0:3] == '#EN'):
                energy = i[8:]
            elif (i[0:3] == '#L '):
                lab = i[3:]
                dataStart = lineno + 1
                if cmnd.split()[0] == 'rodscan':
                    lineno = lineno + 1
                    buffOld = buffNew = allLines.next()
                    if buffNew[0:3] != '\n' and buffNew[0:3] != '#C ':
                        lStart = buffOld.split()[2]
                    else:
                        lStart = lStop = cmnd.split()[3]
                    while (buffNew[0:3] != '\n' and buffNew[0:3] != '#C '):
                        buffOld = buffNew
                        lineno = lineno + 1
                        buffNew = allLines.next()
                        """
                        try:
                            buffNew = allLines.next()
                        except:
                            print '\n\n'
                            print 'excepting on line', lineno
                            buffNew = '\n'
                        """
                    if buffOld[0:3] != '\n' and buffOld[0:3] != '#C ':
                        lStop = buffOld.split()[2]
                    dataStop = lineno
                    if buffNew[0:3] == '#C ':
                        if buffNew.split()[7] == 'aborted':
                            aborted = True
                else:
                    lineno = lineno + 1
                    buff = allLines.next()
                    while buff[0:3] != '\n' and buff[0:3] != '#C ':
                        lineno = lineno + 1
                        buff = allLines.next()
                    dataStop = lineno
                    if buff[0:3] == '#C ':
                        if buff.split()[7] == 'aborted':
                            aborted = True
                self._summary.append({'index':index,
                                     'nl_start':n_sline,
                                     'cmd':cmnd,
                                     'date':date,
                                     'time':xtime,
                                     'G':Gvals,
                                     'Q':q,
                                     'mot_names':mnames,
                                     'P':Pvals,
                                     'ncols':ncols,
                                     'labels':lab,
                                     'atten':atten,
                                     'energy':energy,
                                     'lineno':lineno,
                                     'lStart':lStart,
                                     'lStop':lStop,
                                     'aborted':aborted,
                                     'dataStart':dataStart,
                                     'dataStop':dataStop})
                (cmnd,date,xtime,Gvals,q,Pvals,atten,energy,lab,lStart,lStop,aborted) = (None,None,None,None,None,None,None,None,None,None,None,False)
                (index, ncols, n_sline, dataStart, dataStop) = (0,0,0,0,0)

        self.min_scan = self._summary[0]['index']
        self.max_scan = self._summary[0]['index']
        for i in self._summary:
            k = i['index']
            if (k > self.max_scan): self.max_scan = k
            if (k < self.min_scan): self.min_scan = k

    def scan_min(self):
        """
        get the minimum scan number
        """
        self.read()
        return self.min_scan

    def scan_max(self):
        """
        get the max scan number
        """
        self.read()
        return self.max_scan

    def nscans(self):
        """
        get the number of scans
        """
        self.read()
        return len(self._summary)

    def _check_range(self,i):
        """
        check if scan number is in range
        """
        self.read()
        j = True
        if ((i > self.max_scan) or (i < self.min_scan)): j = False
        return j

    #def scan_info(self):
    #    self.read()
    #    return self._summary

    #def get_summary(self,sc_num):
    def scan_info(self,sc_num):
        """
        return the scan info in a dictionary
        """
        self.read()
        for s in self._summary:
            if (sc_num == s['index']):
                return s
        return None
    
    def scan_data(self, sc_num):
        """
        return just the column data from the scan 
        """
        self.read()
        s = self.scan_info(sc_num)
        if (s == None): return None
        dat = []
        nl  = s['lineno']
        for i in (self._lines[nl:]):
            if (i[0:3] == '#S '):
                break
            elif (i[0:1] ==  '#'):
                x = 1
            elif (len(i)  > 3):
                q = i.split()
                dat.append(map(float,q))
        return dat

    def scan_dict(self, sc_num):
        """
        return scan information and data in a dictionary 
        """
        self.read()
        sc_dict = {'file':self.fname,
                   'index':sc_num,
                   'cmd':'',
                   'date':'',
                   'G':[],
                   'Q':[],
                   'P':{},
                   'ATTEN':[],
                   'ENERGY':[],
                   'labels':[],
                   'ncol':0,
                   'nrow':0,
                   'data':{}
                   }
        s = self.scan_info(sc_num)
        if (s == None): return sc_dict
        dat = self.scan_data(sc_num)
        
        # parse the various data into the dict
        sc_dict['cmd']  = s['cmd']
        sc_dict['date'] = s['date']
        if s['G'] != None: sc_dict['G']    = map(float,s['G'].split())
        if s['Q'] != None: sc_dict['Q']    = map(float,s['Q'].split())
        if s['atten'] != None: sc_dict['ATTEN'] = map(int,s['atten'].split())
        if s['energy'] != None: sc_dict['ENERGY'] = map(float,s['energy'].split())
        # get the motor positions
        p_dict = {}
        m_names = s['mot_names'].split()
        p_vals  = s['P'].split()
        if len(m_names) != len(p_vals):
            print "Mismatch in Motor Names and Motor Values"
        else:
            for j in range(len(m_names)):
                p_dict.update({m_names[j]:float(p_vals[j])})
        sc_dict['P'] = p_dict
        #
        lbls = s['labels'].split()
        sc_dict['labels'] = lbls
        ncol = len(lbls)
        nrow = len(dat)
        sc_dict['ncol']   = ncol
        sc_dict['nrow']   = nrow
        # data
        data_dict = {}
        for j in range(ncol):
            xx = []
            for k in range(nrow):
                xx.append(dat[k][j])
            data_dict.update({lbls[j]:xx})
        sc_dict['data'] = data_dict
        # all done
        return sc_dict

    def list_scans(self):
        """
        return a list of scans in the file 
        """
        self.read()
        sc_list = [] 
        for s in self._summary:
            line = "%s:%4.4i  (%s)\n   %s" % (self.fname,
                                              s['index'],
                                              s['date'],
                                              s['cmd'])
            sc_list.append(line)
        return sc_list


#######################################################################
#######################################################################
#######################################################################
"""
Some notes on spec G array
--------------------------

The following code snipet from spec standard.mac
shows how the #G lines in the data file are defined 

{{{
** in the _head macro:
 
  _head_par G 0
  _head_par U 1
  _head_par UB 3
  _head_par Q 4

** the definition of _head_par 
# $1 is parameter name.  $2 is #G suffix
def _head_par '
    if (whatis("$1[0]")) {
            local i
            printf("#G$2")
            for (i=0;;i++)
                    if (whatis(sprintf("$1[%d]", i)))
                            printf(" %.10g", $1[i])
                    else break
            printf("\n")
    }
'
}}}

This crazy bit of code correlates the spec internal arrays
G, U, UB and Q with the #G values in the file:
#G0 = G   --> Geometry information
#G1 = U   --> Orientation and lattice information
#G3 = UB  --> The UB matrix (9 parameters for 3x3 matrix)
#G4 = Q   --> More Geometry information

The contents of these arrays is outlined for psic geometry
in the psic_clean_start_idc.mac and described below.

Note preceding number below indicates the index in the
G array returned from read_spec.  The index jumps around
to account for padding that is in the spec arrays.
Remember this is specific to psic geometry...

########## G array ####################
# sector preference flags
0.  g_prefer=0       # G[0], Holds sector preference value. 1=vert, 2=horz, 
1.  g_sect=0         # G[1], Holds sector mode. 0-16, 0=use ranking method
                     #       if g_sect=g_prefer=0 then no angle transforms are done
2.  g_frz=1          # G[2], Nonzero when frozen mode is on. 
                     #       ie use the values in the F_... vars for angle constraints
                     #       if this is zero it just uses the current values
# reference vector
3.  g_haz=0          # G[3], H of reference vector
4.  g_kaz=0          # G[4], K of reference vector
5.  g_laz=1          # G[5], L of reference vector
# zone vectors
6.  g_zh0=0          # G[6], H of first zone mode vector
7.  g_zk0=0          # G[7], K of first zone mode vector
8.  g_zl0=0          # G[8], L of first zone mode vector
9.  g_zh1=0          # G[9], H of second zone mode vector
10. g_zk1=0          # G[10], K of second zone mode vector
11. g_zl1=0          # G[11], L of second zone mode vector
12. g_kappa = 50.031 # G[12], Kappa tilt angle
#
15. g_sigtau=0       # G[15], Nonzero when use SIGMA_AZ and TAU_AZ for the reference vector
                     #        orientation instead of using the g_haz, g_kaz, g_laz
# set geometry constraint mode
# these settings give a four-circle alpha-fixed mode if F_NU = F_MU = 0
16. g_mode1=2        # G[16], Nu-fixed
17. g_mode2=2        # G[17], Alpha-fixed
18. g_mode3=2        # G[18], Mu-fixed
19. g_mode4=0        # G[19], --
20. g_mode5=0        # G[20], --

############ U array  ##############
# lattice parameters
22. g_aa=1.54        # U[0], real-space a
23. g_bb=1.54        # U[1], real-space b
24. g_cc=1.54        # U[2], real-space c
25. g_al=90          # U[3], real-space alpha
26. g_be=90          # U[4], real-space beta
27. g_ga=90          # U[5], real-space gamma
28. g_aa_s=4.07999   # U[6], recip-space a* (x 2pi)
29. g_bb_s=4.07999   # U[7], recip-space b* (x 2pi)
30. g_cc_s=4.07999   # U[8], recip-space c* (x 2pi)
31. g_al_s=90        # U[9], recip-space alpha* 
32. g_be_s=90        # U[10], recip-space beta* 
33. g_ga_s=90        # U[11], recip-space gamma*
# orientation
# Note these settings define a unit orientation matrix.  
# Therefore, the azimuth vector (001) points along the 
# diffractometer Z-axis for all angles = 0. 
# reflections
# or0
34. g_h0=1           #  U[12], h
35. g_k0=0           #  U[13], k
36. g_l0=0           #  U[14], l
# or1
37. g_h1=1           #  U[15], h
38. g_k1=0           #  U[16], k
39. g_l1=1           #  U[17], l
# angle settings
# or0
40. g_u00=60         #  U[18], del
41. g_u01=30         #  U[19], eta
42. g_u02=0          #  U[20], chi
43. g_u03=0          #  U[21], phi
44. g_u04=0          #  U[22], nu
45. g_u05=0          #  U[23], mu
# or1
46. g_u10=90         #  U[24], del
47. g_u11=45         #  U[25], eta
48. g_u12=45         #  U[26], chi
49. g_u13=0          #  U[27], phi
50. g_u14=0          #  U[28], nu
51. g_u15=0          #  U[29], mu
# lambdas for ors
52. g_lambda0=1.54   #  U[30], lambda or0 
53. g_lambda1=1.54   #  U[31], lambda or1

############ Q array ######################
# These are recalc when do a calcHKL (or calc(2)).
63. H=0.0            # Q[0], H coord
64. K=0.0            # Q[1], K coord
65. L=0.0            # Q[2], L coord
66. LAMBDA=1.54      # Q[3], Wavelength (only recalculated if have a mono motor defined)
67. ALPHA=0          # Q[4], Angle between the reference vector and the xz-plane (incidence angle) 
68. BETA=0           # Q[5], Exit angle when the reference vector is the surface normal
69. OMEGA=0          # Q[6], Angle btwn Q and the plane perpendicular to the chi axis
70. TTH=0            # Q[7], Scattering angle
71. PSI=0            # Q[8], Azimuthal angle of the reference vector wrt Q and the scatt-plane
72. TAU=0            # Q[9], Longitudinal angle of the reference vector wrt Q and the scatt-plane
73. QAZ=0            # Q[10], Angle btwn Q and the yz-plane
74. NAZ=0            # Q[11], Angle btwn the reference vector and the yz-plane
75. SIGMA_AZ=0       # Q[12], Angle to specify reference vector, = -flat_chi
76. TAU_AZ=0         # Q[13], Angle to specify reference vector, = -flat_phi
# Values for the frozen angles, 
# These are only used if g_frz=1, otherwise the current values are used
77. F_ALPHA=0        # Q[14], Frozen value of ALPHA for alpha-fixed mode.
78. F_BETA=0         # Q[15], Frozen value of BETA for beta-fixed mode
79. F_OMEGA=0        # Q[16], Frozen value of OMEGA for omega-fixed mode
80. F_PSI=0          # Q[17], Frozen value of PSI for psi-fixed mode
81. F_NAZ=0          # Q[18], Frozen value of NAZ for naz-fixed mode
82. F_QAZ=0          # Q[19], Frozen value of QAZ for qaz-fixed mode
83. F_DEL=0          # Q[20], Frozen value of A[del] for delta-fixed mode
84. F_ETA=0          # Q[21], Frozen value of A[eta] for eta-fixed mode
85. F_CHI=0          # Q[22], Frozen value of A[chi] for chi-fixed mode
86. F_PHI=0          # Q[23], Frozen value of A[phi] for phi-fixed mode
87. F_NU=0           # Q[24], Frozen value of A[nu] for nu-fixed mode
88. F_MU=0           # Q[25], Frozen value of A[mu] for mu-fixed mode
89. F_CHI_Z=0        # Q[26], Value calculated for A[chi] in zone mode
90. F_PHI_Z=0        # Q[27], Value calculated for A[phi] in zone mode
# Cut points for the axis
91. CUT_DEL=-180     # Q[28], Cut point for del circle
92. CUT_ETA=-180     # Q[29], Cut point for eta circle
93. CUT_CHI=-180     # Q[30], Cut point for chi circle
94. CUT_PHI=-180     # Q[31], Cut point for phi circle
95. CUT_NU=-180      # Q[32], Cut point for mu circle
96. CUT_MU=-180      # Q[33], Cut point for nu circle
97. CUT_KETA=-180    # Q[34], Cut point for keta circle
98. CUT_KAP=-180     # Q[35], Cut point for kap circle
99. CUT_KPHI=-180    # Q[36], Cut point for kphi circle
"""
