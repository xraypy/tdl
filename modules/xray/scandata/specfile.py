#######################################################################
"""
Matt Newville, Tom Trainor (fftpt@uaf.edu)
Class to handle spec files

Modifications:
--------------


"""
#######################################################################

import numpy as num
import os
import types

#######################################################################
class SpecFile:
    def __init__(self, fname):
        self.path, self.fname = os.path.split(fname)
        self.max_scan = 0
        self.min_scan = 0
        self._mtime    = 0
        self._lines    = []
        self._summary  = []
        self._ok       = False
        self.read()

    def __repr__(self):
        self.read()
        lout = "Spec file: %s" % self.fname
        lout = "%s\nPath: %s" % (lout, os.path.join(self.path))
        lout = "%s\nFirst scan number: %i" % (lout,self.min_scan)
        lout = "%s\nLast scan number:  %i" % (lout,self.max_scan)
        lout = "%s\nLast scan: %s" % (lout, self._summary[self.max_scan-1]['date'])
        return lout

    def read(self):
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
        lineno = 0
        (mnames,cmnd,date,xtime,Gvals,q,Pvals,atten,lab) = (None,None,None,None,None,None,None,None,None)
        (index, ncols, n_sline) = (0,0,0)
        for i in self._lines:
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
                                     'lineno':lineno})
                (cmnd,date,xtime,Gvals,q,Pvals,atten,energy,lab) = (None,None,None,None,None,None,None,None,None)
                (index, ncols, n_sline) = (0,0,0)

        self.min_scan = self._summary[0]['index']
        self.max_scan = self._summary[0]['index']
        for i in self._summary:
            k = i['index']
            if (k > self.max_scan): self.max_scan = k
            if (k < self.min_scan): self.min_scan = k

    def scan_min(self):
        self.read()
        return self.min_scan

    def scan_max(self):
        self.read()
        return self.max_scan

    def nscans(self):
        self.read()
        return len(self._summary)

    def _check_range(self,i):
        self.read()
        j = True
        if ((i > self.max_scan) or (i < self.min_scan)): j = False
        return j

    #def scan_info(self):
    #    self.read()
    #    return self._summary

    #def get_summary(self,sc_num):
    def scan_info(self,sc_num):
        self.read()
        for s in self._summary:
            if (sc_num == s['index']):
                return s
        return None
    
    def scan_data(self, sc_num):
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
        sc_dict['G']    = map(float,s['G'].split())
        sc_dict['Q']    = map(float,s['Q'].split())
        sc_dict['ATTEN'] = map(int,s['atten'].split())
        sc_dict['ENERGY'] = map(float,s['energy'].split())
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
        self.read()
        sc_list = [] 
        for s in self._summary:
            line = "%s:%4.4i  (%s)\n   %s" % (self.fname,
                                              s['index'],
                                              s['date'],
                                              s['cmd'])
            sc_list.append(line)
        return sc_list
