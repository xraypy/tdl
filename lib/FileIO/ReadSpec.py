# Read spec files
# From Matt Newville
# Modifed by Tom Trainor

import string


class SpecFile:
    def __init__(self,fname=None):
        self.fname    = fname
        self.lines    = []
        self.summary  = []
        self.max_scan = 0
        self.min_scan = 0
        if (fname !=None):
            self.read_specfile(fname)       

    def read_specfile(self,fname):
        self.fname = fname
        try:
            f  = open(fname)
        except IOError:
            print  'error reading file ', fname
            return -1
        self.lines = f.readlines()
        f.close()
        self.summarize()
        return 0

    def summarize(self):
        lineno = 0
        (mnames,cmnd,date,xtime,Gvals,q,Pvals,lab) = (None,None,None,None,None,None,None,None)
        (index, ncols, n_sline) = (0,0,0)
        for i in self.lines:
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
            elif (i[0:3] == '#L '):
                lab = i[3:]
                self.summary.append({'index':index,
                                     'nl_start':n_sline,
                                     'cmnd':cmnd,
                                     'date':date,
                                     'time':xtime,
                                     'Gvals':Gvals,
                                     'q':q,
                                     'mot_names':mnames,
                                     'Pvals':Pvals,
                                     'ncols':ncols,
                                     'labels':lab,
                                     'nl_label':lineno})
                (cmnd,date,xtime,Gvals,q,Pvals,lab) = (None,None,None,None,None,None,None)
                (index, ncols, n_sline) = (0,0,0)

        self.min_scan = self.summary[0]['index']
        self.max_scan = self.summary[0]['index']
        for i in self.summary:
            k = i['index']
            if (k > self.max_scan): self.max_scan = k
            if (k < self.min_scan): self.min_scan = k

    def scan_min(self):
        return self.min_scan

    def scan_max(self):
        return self.max_scan

    def check_range(self,i):
        j = 1
        if ((i > self.max_scan) or (i < self.min_scan)): j = 0
        return j

    def nscans(self):
        return len(self.summary)

    def scan_info(self):
        return self.summary

    def get_summary(self,sc_num):
        for i in self.summary:
            if (sc_num == i['index']):
                return i
        return None
            
    def get_data(self, sc_num):
        s = self.get_summary(sc_num)
        if (s == None): return None
        dat = []
        nl  = s['nl_label']
        for i in (self.lines[nl:]):
            if (i[0:3] == '#S '):
                break
            elif (i[0:1] ==  '#'):
                x = 1
            elif (len(i)  > 3):
                q = i.split()
                dat.append(map(float,q))
        return dat

    def get_scan_dict(self, sc_num):
        sc_dict = {'file':self.fname,
                   'index':sc_num,
                   'cmd':'',
                   'date':'',
                   'G':[],
                   'Q':[],
                   'P':{},
                   'ATTEN':[],
                   'labels':[],
                   'ncol':0,
                   'nrow':0,
                   'data':{}
                   }
        s = self.get_summary(sc_num)
        if (s == None): return sc_dict
        dat = self.get_data(sc_num)
        
        # parse the various data into the dict
        sc_dict['cmd']  = s['cmnd']
        sc_dict['date'] = s['date']
        sc_dict['G']    = map(float,string.split(s['Gvals']))
        sc_dict['Q']    = map(float,string.split(s['q']))
        # get the motor positions
        p_dict = {}
        m_names = string.split(s['mot_names'])
        p_vals  = string.split(s['Pvals'])
        if len(m_names) != len(p_vals):
            print "Mismatch in Motor Names and Motor Values"
        else:
            for j in range(len(m_names)):
                p_dict.update({m_names[j]:float(p_vals[j])})
        sc_dict['P'] = p_dict
        #
        lbls = string.split(s['labels'])
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
            
            
 