##########################################################################
"""
# physcon Physical constants
# Note: type physcon.help() (after import physcon)

# Copyright 2004,2007 Herman J.C. Berendsen, <www.hjcb.nl/python>;
# email <herman@hjcb.nl>. The data originate from the CODATA 2006 release,
# see <http://www.physics.nist.gov/PhysRefData/contents.html>.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version. See <http://www.gnu.org/licenses/gpl.txt>
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
"""
##########################################################################

from math import pi

# define light velocity, because we need it in calculations
cloc =299792458.

# dictionary of physical constants in SI units. Values CODATA 2006 http://www.physics.nist.gov/cuu/Constants/
# each item: [description (string), symbol (string), value (float), sd (float), relat. sd (float),
#             value(sd) unit (string), source (string)]
all={'lightvel':['velocity of light in vacuum','c',cloc,0., 0.,'299 792 458(ex) m/s', 'CODATA 2006'],
     'planck':["Planck's constant",'h',6.62606896e-34,3.3e-41,5.0e-8,'6.626 068 96(33) e-34 J s', 'CODATA 2006'],
     'dirac':["Dirac's constant = h/(2 pi)",'hbar',1.054571628e-34,5.3e-42,5.0e-8,'1.054 571 628(53) e-34 J s', 'CODATA 2006 '],
     'magn-const':['magnetic permeability of vacuum','mu_0',4e-7*pi,0.,0.,'1.256 637 061... e-6 N A^-2',''],
     'elec-const':['dielectric permittivity of vacuum','eps_0',1.e7/(4*pi*cloc*cloc),0.,0.,'8.854 187 817... e-12 F/m',''],
     'gravit':['Newton constant of gravitation','G',6.67428e-11, 6.7e-15,1.0e-4,'6.674 28(67) e-11 m^3 kg^-1 s^-2','CODATA 2006'],
     'charge-e':['elementary charge','e',1.60217653e-19,1.4e-26,8.5e-8,'1.602 176 53(14) e-19 C','CODATA 2006'],
     'mass-e':['electron mass','m_e',9.10938215e-31,4.5e-38,5.0e-8,'9.109 382 15(45) e-31 kg','CODATA 2006'],
     'mass-e/u':['electron mass in u','m_e_u',5.4857990943e-4,2.3e-13,4.2e-10,'5.485 799 0945(23) u','CODATA 2006'],
     'mass-p':['proton mass','m_p',1.672621637e-27,8.3e-35,5.0e-8,'1.672 621 637(83) e-27 kg','CODATA 2006'],
     'mass-p/u':['proton mass in u','m_p_u',1.00727646677,1.0e-10,1.0e-10,'1.007 276 466 77(10) u','CODATA 2006'],
     'mass-n':['neutron mass','m_n',1.674927211e-27,8.4e-35,5.0e-8,'1.674 927 211(84) e-27 kg','CODATA 2006'],
     'mass-n/u':['neutron mass in u','m_n_u',1.00866491597,4.3e-10,4.3e-10,'1.008 664 915 97(43) u','CODATA 2006'],
     'mass-d':['deuteron mass','m_d',3.34358320e-27,1.7e-34,5.0e-8,'3.343 583 20(17) e-27 kg','CODATA 2006'],
     'mass-d/u':['deuteron mass in u','m_d_u',2.013553212724,7.8e-11,3.9e-11,'2.013 553 212 724(78) u','CODATA 2006'],
     'mass-mu':['muon mass','m_m',1.88353130e-28,1.1e-35,5.6e-8,'1.883 531 30(11) e-28 kg','CODATA 2006'],
     'mass-mu/u':['muon mass in u','m_m_u',0.1134289256,2.9e-9,2.5e-8,'0.113 428 9256(29) u','CODATA 2006'],
     'ratio-me/mp':['electron/proton mass ratio','ratio_memp',5.4461702177e-4,2.4e-13,4.3e-10,'5.446 170 2177(24) e-4','CODATA 2006'],
     'ratio-mp/me':['proton/electron mass ratio','ratio_mpme',1836.15267247,8.0e-7,4.3e-10,'1836.152 672 47(80)','CODATA 2006'],
     'amu':['unified atomic mass unit = 1/12 m(12C)','u',1.660538782e-27,8.3e-35,5.0e-8,'1.660 538 782(83) e-27 kg','CODATA 2006'],
     'avogadro':['Avogadro constant','N_A',6.02214179e23,3.0e16,5.0e-8,'6.022 141 79(30) e23 mol^-1','CODATA 2006'],
     'boltzmann':['Boltzmann constant','k_B',1.3806504e-23,2.4e-29,1.7e-6,'1.380 6504(24) J/K','CODATA 2006'],
     'gas':['molar gas constant = N_A k_B','R',8.314472,1.5e-5,1.7e-6,'8.314 472(15) J mol^-1 K^-1','CODATA 2006'],
     'faraday':['Faraday constant = N_A e','F',96485.3399,2.4e-3,2.5e-8,'96 485.3399(24) C/mol','CODATA 2006'],
     'bohrradius':['Bohr radius = 4 pi eps_0 hbar^2/(m_e e^2)','a_0',5.2917720859e-11,3.6e-20,6.8e-10,'0.529 177 208 59(36) e-10 m','CODATA 2006'],
     'magflux-qu':['magnetic flux quantum = h/(2 e)','Phi_0',2.067833667e-15,5.2e-23,2.5e-8,'2.067 833 667(52) Wb','CODATA 2006'],
     'conduct-qu':['conductance quantum = 2 e^2/h','G_0',7.7480917004e-5,5.3e-14,6.8e-10,'7.748 091 7004(53) e-5 S','CODATA 2006'],
     'josephson':['Josephson constant = 2 e/h','K_J',4.83597891e14, 1.2e6,2.5e-8,'4.835 978 91(12) e14 Hz/V','CODATA 2006'],
     'bohrmagn':['Bohr magneton = e hbar/(2 m_e)','mu_B',9.27400915e-24,2.3e-31,2.5e-8,'9.274 009 15(23) e-24 J/T','CODATA 2006'],
     'nuclmagn':['nuclear magneton = e hbar/(2 m_p)','mu_N',5.05078324e-27,1.3e-34,2.5e-8,'5.050 783 24(13) e-27 J/T','CODATA 2006'],
     'magnmom-e':['electron magnetic moment','mu_e',-9.28476377e-24,2.3e-31,2.5e-8,'-9.284 763 77(23) e-24 J/T','CODATA 2006'],
     'magnmom-p':['proton magnetic moment','mu_p',1.41060662e-26,3.7e-33,2.5e-8,'1.410 606 662(37) e-26 J/T','CODATA 2006'],
     'gfactor-e':['electron g-factor','g_e',-2.0023193043622,1.5e-12,7.4e-13,'-2.002 319 304 3622(15)','CODATA 2006'],
     'gfactor-p':['proton g-factor','g_p',5.585694713, 4.6e-8,8.2e-9,'5.585 694 713(46)','CODATA 2006'],
     'alpha':['fine-structure constant = e^2/(4 pi eps_0 hbar c)','alpha',7.2973525376e-3,5.0e-12,6.8e-10,'7.297 352 5376(50) e-3','CODATA 2006'],
     'alpha-1':['inverse fine-structure constant = 4 pi eps_0 hbar c/e^2','',137.035999679,9.4e-8,6.8e-10,'137.035 999 679(94)','CODATA 2006'],
     'gyromagratio-p':['proton gyromagnetic ratio','gamma_p',2.675222099e8,7.0,2.6e-8,'2.675 222 099(70) e8 s^-1 T^-1','CODATA 2006'],
     'magres-p':['magnetic resonance frequency proton = gamma_p/(2*pi)','',4.25774821e7,1.1,2.6e-8,'42.577 4821(11) MHz/T','CODATA 2006'],
     'rydberg':['Rydberg constant = alpha^2 m_e c/(2 h)','R_infty',10973731.568527,7.3e-5,6.6e-12,'10 973 731.568 527(73) m^-1','CODATA 2006'],
     'stefan-boltzm':['Stefan-Boltzmann constant = pi^2 k^4/(60 hbar^3 c^2)','sigma',5.670400e-8,4.0e-13,7.0e-6,'5.670 400(40) e-8 W m^-2 K^-4','CODATA 2006']}             
     

# many common values are also available as global constants:
global alpha,a_0,c,e,eps_0,F,G,g_e,g_p,gamma_p,h,hbar,k_B
global m_d,m_e,m_n,m_p,mu_B,mu_e,mu_N,mu_p,mu_0,N_A,R,sigma,u
alpha = all['alpha'][2]
a_0 =  all['bohrradius'][2]
c = cloc
e =  all['charge-e'][2]
eps_0 =  all['elec-const'][2]
F =  all['faraday'][2]
G =  all['gravit'][2]
g_e =  all['gfactor-e'][2]
g_p =  all['gfactor-p'][2]
gamma_p =  all['gyromagratio-p'][2]
h = all['planck'][2]
hbar = all['dirac'][2]
k_B =  all['boltzmann'][2]
m_d =  all['mass-d'][2]
m_e =  all['mass-e'][2]
m_n =  all['mass-n'][2]
m_p =  all['mass-n'][2]
mu_B =  all['bohrmagn'][2]
mu_e =  all['magnmom-e'][2]
mu_N =  all['nuclmagn'][2]
mu_p =  all['magnmom-p'][2]
mu_0 =  all['magn-const'][2]
N_A =  all['avogadro'][2]
R =  all['gas'][2]
sigma =  all['stefan-boltzm'][2]
u =  all['amu'][2]


def help():
    print 'Available functions:'
    print '[note: key must be a string, within quotes!]' 
    print '  value(key) returns value (float)'
    print '  sd(key)    returns standard deviation (float)'
    print '  relsd(key) returns relative standard deviation (float)'
    print '  descr(key) prints description with units\n'
    print 'Available global variables:'
    print '  alpha, a_0, c, e, eps_0, F, G, g_e, g_p, gamma_p, h, hbar, k_B'
    print '  m_d, m_e, m_n, m_p, mu_B, mu_e, mu_N, mu_p, mu_0, N_A, R, sigma, u\n'
    allkeys=all.keys()
    allkeys.sort()
    print 'Available keys:'
    print allkeys

def value(key):
    return all[key][2]

def sd(key):
    return all[key][3]

def relsd(key):
    return all[key][4]

def descr(key):
    print 'Description of ',key,':'
    print '  Name:               ',all[key][0]
    print '  Symbol (if avail.): ',all[key][1]
    print '  Value:              ',all[key][2]
    print '  Standard deviation: ',all[key][3] 
    print '  Relative stdev:     ',all[key][4]
    print '  value(sd) unit:     ',all[key][5]
    print '  Source:             ',all[key][6],'\n'
    
