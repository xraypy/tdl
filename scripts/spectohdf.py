'''
Spec to HDF5 converter
Author: Craig Biwer (cbiwer@uchicago.edu)
12/19/2011
'''

import h5py
import sys
import os
import numpy as num
from PIL import Image

import time

#Parses the specfile
def summarize(lines):

        #This list corresponds to the current set of parameters
        #written by spec to lines beginning #G. If that changes,
        #at least this list will need to be modified, if not more
        #code below.
        gLabs = ['g_prefer', 'g_sect', 'g_frz', 'g_haz', 'g_kaz', 'g_laz', 'g_zh0', 'g_zk0', 'g_z10',
                    'g_zh1', 'g_zk1', 'g_zl1', 'g_kappa', 'g_13', 'g_14', 'g_sigtau', 'g_mode1',
                    'g_mode2', 'g_mode3', 'g_mode4', 'g_mode5', 'g_21', 'g_aa', 'g_bb', 'g_cc',
                    'g_al', 'g_be', 'g_ga', 'g_aa_s', 'g_bb_s', 'g_cc_s', 'g_al_s', 'g_be_s', 'g_ga_s',
                    'g_h0', 'g_k0', 'g_l0', 'g_h1', 'g_k1', 'g_l1', 'g_u00', 'g_u01', 'g_u02', 'g_u03',
                    'g_u04', 'g_u05', 'g_u10', 'g_u11', 'g_u12', 'g_u13', 'g_u14', 'g_u15', 'g_lambda0',
                    'g_lambda1', 'g_54', 'g_55', 'g_56', 'g_57', 'g_58', 'g_59', 'g_60', 'g_61', 'g_62',
                    'g_H', 'g_K', 'g_L', 'g_LAMBDA', 'g_ALPHA', 'g_BETA', 'g_OMEGA', 'g_TTH', 'g_PSI',
                    'g_TAU', 'g_QAZ', 'g_NAZ', 'g_SIGMA_AZ', 'g_TAU_AZ', 'g_F_ALPHA', 'g_F_BETA',
                    'g_F_OMEGA', 'g_F_PSI', 'g_F_NAZ', 'g_F_QAZ', 'g_F_DEL', 'g_F_ETA', 'g_F_CHI',
                    'g_F_PHI', 'g_F_NU', 'g_F_MU', 'g_F_CHI_Z', 'g_F_PHI_Z', 'CUT_DEL', 'CUT_ETA',
                    'CUT_CHI', 'CUT_PHI', 'CUT_NU', 'CUT_MU', 'CUT_KETA', 'CUT_KAP', 'CUT_KPHI', 'g_100',
                    'g_101', 'g_102', 'g_103', 'g_104', 'g_105', 'g_106', 'g_107', 'g_108', 'g_109',
                    'g_110', 'g_111']

        summary = []
        lineno = 0
        (mnames,cmnd,sType,date,xtime,Gvals,q,Pvals,atten,energy,lab,lStart,lStop,aborted) = (None,None,None,None,None,None,None,None,None,None,None,None,None,False)
        pointData = []
        (index, ncols, n_sline, dataStart, dataStop) = (0,0,0,0,0)
        lines.append('\n')
        allLines = iter(lines)
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
                sType = cmnd.split()[0]
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
                lab = i[3:].split()
                dataStart = lineno + 1
                if cmnd.split()[0] == 'rodscan':
                    lineno = lineno + 1
                    buffOld = buffNew = allLines.next()
                    if buffNew[0:3] != '\n' and buffNew[0:3] != '#C ':
                        lStart = buffOld.split()[2]
                    while (buffNew[0:3] != '\n' and buffNew[0:3] != '#C '):
                        pointData.append(map(float, buffNew.split()))
                        buffOld = buffNew
                        lineno = lineno + 1
                        buffNew = allLines.next()
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
                        pointData.append(map(float, buff.split()))
                        lineno = lineno + 1
                        buff = allLines.next()
                    dataStop = lineno
                    if buff[0:3] == '#C ':
                        if buff.split()[7] == 'aborted':
                            aborted = True
                currentDict = {'index':index,
                                     'nl_start':n_sline,
                                     'cmd':cmnd,
                                     'sType':sType,
                                     'date':date,
                                     'time':xtime,
                                     'mnames':mnames.split(),
                                     'P':map(float, Pvals.split()),
                                     'gLabs':gLabs,
                                     'G':map(float, Gvals.split()),
                                     'Q':q,
                                     'ncols':ncols,
                                     'labels':lab,
                                     'atten':atten,
                                     'energy':energy,
                                     'lineno':lineno,
                                     'lStart':lStart,
                                     'lStop':lStop,
                                     'aborted':aborted,
                                     'dataStart':dataStart,
                                     'dataStop':dataStop,
                                     'pointData':pointData,}
                #for listSpot in range(len(mnames.split())):
                #    currentDict[mnames.split()[listSpot-1]] = float(Pvals.split()[listSpot-1])
                #for listSpot in range(len(gLabs)):
                #    currentDict[gLabs[listSpot-1]] = float(Gvals.split()[listSpot-1])
                summary.append(currentDict)
                (cmnd,sType,date,xtime,Gvals,q,Pvals,atten,energy,lab,lStart,lStop,aborted) = (None,None,None,None,None,None,None,None,None,None,None,None,False)
                pointData = []
                (index, ncols, n_sline, dataStart, dataStop) = (0,0,0,0,0)

        min_scan = summary[0]['index']
        max_scan = summary[0]['index']
        for i in summary:
            k = i['index']
            if (k > max_scan): max_scan = k
            if (k < min_scan): min_scan = k
            
        return summary
    
def readImage(file):
    try:
        im  = Image.open(file)
        arr = num.fromstring(im.tostring(), dtype='int32')
        arr.shape = (im.size[1],im.size[0])
        return arr
    except:
        print "Error reading file: %s" % file
        return None

#Main conversion function, takes user input and converts specified .spc file to .h5
def specToHDF(args):
    if len(args) == 0:
        print 'Current directory: ', os.getcwd()
        print 'Please enter the path to the specfile:'
        input = raw_input('> ')
        if not os.path.isfile(input):
            print 'Error: file not found'
            return
        print 'Please enter the name of the HDF file:'
        output = raw_input('> ')
        if output == '':
            output = 'default.h5'
        if output[-3:] != '.h5':
            output = output + '.h5'
        if os.path.isfile(output):
            choice = raw_input('File already exists: (A)ppend, (O)verwrite, or (C)ancel (A/O/C)? ').lower()
            if choice == 'c':
                return
            elif choice == 'a':
                print 'Sorry, not implemented yet.'
                return
    elif len(args) == 1:
        input = args[0]
        print 'Please enter the name of the HDF file:'
        output = raw_input('> ')
        if output == '':
            output = 'default.h5'
        if output[-3:] != '.h5':
            output = output + '.h5'
        if os.path.isfile(output):
            choice = raw_input('File already exists: (A)ppend, (O)verwrite, or (C)ancel (A/O/C)? ').lower()
            if choice == 'c':
                return
            elif choice == 'a':
                print 'Sorry, not implemented yet.'
                return
    elif len(args) == 2:
        input = args[0]
        if not os.path.isfile(input):
            print 'Error: file not found'
            return
        output = args[1]
        if output == '':
            output = 'default.h5'
        if output[-3:] != '.h5':
            output = output + '.h5'
        if os.path.isfile(output):
            choice = raw_input('File already exists: (A)ppend, (O)verwrite, or (C)ancel (A/O/C)? ').lower()
            if choice == 'c':
                return
            elif choice == 'a':
                print 'Sorry, not implemented yet.'
                return
    else:
        input = args[0]
        if not os.path.isfile(input):
            print 'Error: file not found'
            return
        output = args[1]
        if output == '':
            output = 'default.h5'
        if output[-3:] != '.h5':
            output = output + '.h5'
        if os.path.isfile(output):
            choice = raw_input('File already exists: (A)ppend, (O)verwrite, or (C)ancel (A/O/C)? ').lower()
            if choice == 'c':
                return
            elif choice == 'a':
                print 'Sorry, not implemented yet.'
                return
        options = args[2:]
    
    time1 = time.time()
    
    specDir = os.path.split(input)[0]
    specName = os.path.split(input)[-1]
    imageDir = specDir + '\\images\\%s\\' % specName[:-4]
    thisFile = open(input)
    lines = thisFile.readlines()
    thisFile.close()
    summary = summarize(lines)
    
    print imageDir
    #print summary[1].keys()
    
    #print summary[1]['labels']
    thatFile = h5py.File(output)
    masterGroup = thatFile.create_group('MasterCopy')
    specGroup = masterGroup.create_group(specName)
    #for key in summary[1].keys():
    #    if key not in ['labels', 'pointData', 'gLabs', 'mnames', 'G', 'P']:
    #        print key
    #        print type(summary[1][key])
    
    for scan in summary:
        scanGroup = specGroup.create_group(str(scan['index']))
        scanGroup.create_dataset('pointLabs', data=scan['labels'])
        scanGroup.create_dataset('pointData', data=scan['pointData'])
        scanGroup.create_dataset('paramLabs', data=scan['gLabs'] + scan['mnames'])
        scanGroup.create_dataset('paramData', data=scan['G'] + scan['P'])
        for key in scan.keys():
            if key not in ['labels', 'pointData', 'gLabs', 'mnames', 'G', 'P']:
                scanKey = scan[key]
                if scanKey == None:
                    scanKey = 'None'
                scanGroup.attrs[key] = scanKey
        thisDir = imageDir + 'S%03d\\' % scan['index']
        print thisDir
        if os.path.isdir(thisDir):
            imageData = []
            for imageFile in os.listdir(thisDir):
                try:
                    imageValue = None
                    imagePath = os.path.join(thisDir, imageFile)
                    if imagePath[-4:] == '.tif':
                        imageValue = readImage(imagePath)
                    if imageValue != None:
                        imageData.append(imageValue)
                except:
                    print 'Error reading image ' + imagePath
            if imageData != []:
                scanGroup.create_dataset('imageData', data=imageData, compression='szip')
    
    thatFile.close()
    time2 = time.time()
    print (time2-time1)/60

if __name__ == '__main__':    
    args = sys.argv[1:]
    specToHDF(args)
    
