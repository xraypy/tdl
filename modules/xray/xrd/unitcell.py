#!/usr/bin/python
##########################################################################
"""
Matt Newville
Unit cell calculations

Modifications:
--------------


"""
##########################################################################

import string
import Numeric

import Ifeffit
import atomic

##########################################################################
class UnitCell:
    """
    Crystallographic unit cell, using modified Atoms p1.inp file
    
    """
    def __init__(self):
        self.titles  = []
        self.atoms   = []
        self.hkl     = Numeric.zeros(3)
        self.abc     = Numeric.zeros(3)
        self.p1_data  = {'a':0.,'b':0.,'c':0.,
                        'alpha':90.,'beta':90.,'gamma':90.,
                        'space':'p 1', 'core': '',
                        'emin': -500., 'emax':1500., 
                        'edge':'K','rmax':0.}
        self.__file_read = 0
        self.structure_factor = 0.
        self.iff = Ifeffit.Ifeffit(use_numeric=1)


    def do_ifeffit(self,cmd):
        " execute ifefit command "
        # print cmd
        self.iff.ifeffit(cmd)

    def read_p1(self,file='p1.inp'):
        """read a p1.inp file, filling internal data.
        returns 0 on success, 1 on failure."""
        self.filename = file
        try:
            f = open(self.filename,'r')
            lines = f.readlines()
            f.close()
        except:
            print "couldn't read input file %s" % self.filename
            return 1
        
        atoms_mode = 0
        for line in lines:
            # remove end-of-line, any end-of-line comments,
            # and strip whitespace
            line = line[:-1]
            line = line.split('!')[0].strip()
            if (len(line)<1): continue

            if (line[:5] == 'title'):  # title lines
                t = line[5:]
                words = t.split('=')
                if (len(words)>1): t = words[1]
                self.titles.append(t.strip())
            elif (line[:5] == 'atoms'):
                atoms_mode = 1 
            elif (atoms_mode == 1):
                # read atoms list line
                dk_file = '';  bfac = 0.;   occ = 1.
                words = line.split()

                if (len(words) == 4): tag  = words[0]
                if (len(words) > 4):  tag  = words[4]
                if (len(words) > 5):  occ  = float(words[5])
                if (len(words) > 6):  bfac = float(words[6])
                if (len(words) > 7):  
                    dk_file = words[7]
                    print "diffkk file name = %s" %  dk_file
                self.atoms.append({'atom':words[0], 'x':float(words[1]),
                                   'y':float(words[2]), 'z':float(words[3]),
                                   'tag':tag, 'occupancy':occ,
                                   'dwf_b':bfac,'diffkk_file': dk_file })

            elif (line[:5] == 'space'):
                # space group doesn't quite fit with the rest of the
                # parameters below (multi-word, on a line by itself)
                self.p1_data['space'] = line[5:].strip()
            elif (line[:6] == 'miller'):
                # space group doesn't quite fit with the rest of the
                # parameters below (multi-word, on a line by itself)
                line = line.replace(',',' ')
                words = line[6:].split()
                if (len(words) < 3):
                    print "Only %i miller indices given." % (len(words))
                else:
                    self.hkl = Numeric.array(map(float,words[:3]))
            else:
                # read the free format keywords into self.p1_data
                # replace '=' and ',' by ' '.
                line = line.replace('=',' ')
                line = line.replace(',',' ')
                
                # split into words, and reverse the list so that
                # pop() can easily pull off the words
                words = line.split()
                words.reverse()

                while (len(words)>1):
                    # set data in self.p1_data dictionary
                    att = words.pop().strip()
                    if (att in self.p1_data.keys()):
                        try:
                            val = words.pop().strip()
                        except:
                            val = ''
                        self.p1_data[att] = val.strip()
        self.__file_read = 1

        # clean up, making sure numbers are held as floats.
        for p in ('a','b','c','emin','emax',
                  'alpha','beta','gamma', 'rmax'):
            self.p1_data[p] = float(self.p1_data[p])
            
        self.abc = Numeric.array([self.p1_data['a'],
                                  self.p1_data['b'],
                                  self.p1_data['c']])

        return 0


    def get_structure_factor(self,use_diffkk=1):
        """ calculate structure factor"""
        sfact = 0.
        if (self.__file_read == 0):
            print 'need to read a p1.inp file before structure factor calc.'
            return sfact
       
        # determine if we're including resonant corrections (is core atom known?)
        do_resonant = 0
        coretag = self.p1_data['core']
        coresym = ''
        for at in self.atoms:
            if (at['tag'] == coretag):  coresym = at['atom']
        if (coresym  != ''): do_resonant = 1

        if (do_resonant == 1):
            edge = self.p1_data['edge'].lower()
            fcn_edge = atomic.kedge
            if (edge == 'l1'): fcn_edge =atomic.l1edge
            if (edge == 'l2'): fcn_edge =atomic.l2edge
            if (edge == 'l3'): fcn_edge =atomic.l3edge
            fcn_width = atomic.kwidth
            if (edge == 'l1'): fcn_width =atomic.l1width
            if (edge == 'l2'): fcn_width =atomic.l2width
            if (edge == 'l3'): fcn_width =atomic.l3width
            
            e0   = fcn_edge(atomic.z(coresym))
            ewid = fcn_width(atomic.z(coresym))
            estep = ewid
            if (ewid > 5.0):                 estep = 5.0
            if (ewid > 1.0  and ewid < 5.0): estep = 1.0
            if (ewid > 0.5  and ewid < 1.0): estep = 0.5
            # create energy array in ifeffit
            self.do_ifeffit("set st.energy = range(%g,%g,%g)" % (e0+self.p1_data['emin'],
                                                                 e0+self.p1_data['emax'],
                                                                 estep) )
            energy = self.iff.get_array("st.energy")
            i_e0   = 0
            for i in range(len(energy)):
                if (abs(energy[i] - e0) < estep/10.): i_e0 = i
            sfact  = Numeric.zeros(len(energy))
                                
        # q:
        hkl_abc = self.hkl / self.abc
        d2inv   = Numeric.sum(hkl_abc*hkl_abc)     #      (1/d)^2 = Sum((h/a)^2 + (k/b)^2 + (l/c)^2)
        dspace  = Numeric.sqrt(1/d2inv)
        q       = 2 * Numeric.pi * hkl_abc
        qhkl    = 2 * Numeric.pi * self.hkl
        qnorm   = Numeric.sqrt(Numeric.sum(q*q))
        beta    = 0
        f0      = 0
        # sum over atoms in the unit cell
        for atom in self.atoms:
            sym   = atom['atom']
            tag   = atom['tag']
            f     = -1 * atomic.f0_cromer(sym,qnorm)
            z     = atomic.z(sym)
            qdotr = Numeric.dot(qhkl, Numeric.array( (atom['x'], atom['y'], atom['z'])))
            phase = Numeric.exp(1j * qdotr)
            # DWF, occupancy factors
            dwf   = atom['occupancy'] * Numeric.exp(-atom['dwf_b'] * d2inv / 4)

            # f0:   non-resonant structure factor
            f0    = f0   + dwf * f * phase
            # print sym, atom['tag'],  f, phase, q, qdotr

            # beta: sum over resonant sites of phase terms
            if (sym == coresym):   beta = beta + dwf * phase
                
            # resonant corrections to f
            if (do_resonant == 1):
                # first, see if f1&f2 are already calculated:
                f1  = self.iff.get_array("%s.f1" % sym)
                f2  = self.iff.get_array("%s.f2" % sym)
                # if not already calculated, f1 will be 0-length,
                # so calculate f1&f2 now
                group = tag.lower()
                if (len(f1) == 0):
                    if (use_diffkk==1 and len(atom['diffkk_file']) > 2):
                        self.do_ifeffit("""read_data(file=%s,group=tmp,label='en f1 f2') 
                        set %s.f1 = -1*interp(tmp.en,tmp.f1,st.energy)
                        set %s.f2 =    interp(tmp.en,tmp.f2,st.energy)
                        newplot st.energy, %s.f1, key='diffkk f1'
                           plot st.energy, %s.f2, key='diffkk f2'
                        """ % (atom['diffkk_file'],group,group,group,group) )
                        f1 = self.iff.get_array("%s.f1" % group)
                        f2 = self.iff.get_array("%s.f2" % group)
                    else:
                        self.do_ifeffit("f1f2(z=%i, energy=st.energy, group=%s,width=%g)" % (z,sym,ewid))
                        self.do_ifeffit("set %s.f1 = -1*%s.f1" % (sym,sym))
                        f1  = self.iff.get_array("%s.f1" % sym)
                        f2  = self.iff.get_array("%s.f2" % sym)
                        if (sym == coresym):
                            self.do_ifeffit(' newplot st.energy, %s.f1, key="CL f1"' % (sym))
                            self.do_ifeffit(' plot st.energy, %s.f2, key="CL f2"' % (sym))
                                
                f  = f + f1 + 1j*f2
                # f0 includes the resonant part of non-central atoms
                # (take it at e0)
                if (sym != coresym):  f0 = f0 + dwf * f[i_e0] * phase
                    
            # the structure factor sum
            sfact = sfact + dwf * f * phase
        # end structure factor sum
        self.do_ifeffit("show @arrays")
        
        # save results
        self.structure_factor = sfact
        self.beta = beta
        self.f0   = f0
        # save parameters for DAFS ITKK algorithm
        if (do_resonant ==1):
            self.iff.put_array('st.real', sfact.real)
            self.iff.put_array('st.imag', sfact.imag)

        # lorentz correction
        dd = Numeric.arcsin(atomic.hc / (2*dspace*energy[i_e0]))
        self.lorentz_normalization = ((atomic.hc / energy[i_e0])**3) / Numeric.sin(2*dd)

##########################################################################
##########################################################################
if (__name__ == '__main__'):
    import sys
    input_file = 'p1.inp'
    if (len(sys.argv)>1): input_file = sys.argv[1]

    cell = UnitCell()
    cell.read_p1(file=input_file)

    #for i in cell.p1_data.keys():
    #    print i,  ':\t', cell.p1_data[i]

    # for i in cell.atoms:
    #      print "atom  %s (%s) %8.5f %8.5f %8.5f" % (i['atom'], i['tag'],
    #                                                 i['x'], i['y'],i['z'])
        
    cell.get_structure_factor()
    # cell.do_ifeffit(' plot st.energy, st.real, key="real"')
    # cell.do_ifeffit(' plot st.energy, st.imag, key="imag"')
    
    energy = cell.iff.get_array("st.energy")
    st_re  = cell.iff.get_array("st.real")
    st_im  = cell.iff.get_array("st.imag")
    out_file = 'unit_cell.out'
    if 1==1:  #try:
        f = open(out_file,'w')
        f.write("# beta_real = %12.6f\n" %  cell.beta.real)
        f.write("# beta_imag = %12.6f\n" %  cell.beta.imag)
        f.write("# f0_real   = %12.6f\n" %  cell.f0.real)
        f.write("# f0_imag   = %12.6f\n" %  cell.f0.imag)
        f.write("# lorentz   = %12.6f\n" %  cell.lorentz_normalization)
        f.write("#--------------------------------------------\n")
        f.write("#  energy          real            imag\n")
        for i in range(len(energy)):
            f.write("    %12.6f     %12.6f    %12.6f\n"  % ( energy[i], st_re[i], st_im[i]))
        f.close()        
    else: # except: 
        print "Oops, can't write ", out_file

    print "beta_real = %12.6f" %  cell.beta.real
    print "beta_imag = %12.6f" %  cell.beta.imag
    print "f0_real   = %12.6f" %  cell.f0.real
    print "f0_imag   = %12.6f" %  cell.f0.imag
    print "lorentz   = %12.6f" %  cell.lorentz_normalization
    
