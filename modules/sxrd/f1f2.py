"""
Functions used to transform experimental f2 to experimental f1 using
differential kramers kroning transformation and cromer libermann
calculated f1 and f2. 

Authors/Modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""
##############################################################################

import numpy as num
from pylab import *

############################################ f1f2 transformation - Differential Kramers Kroning ####################################
def f1f2(datafile, expfile, e0, e0shift, output ='exp.f1f2', n=30):

    """
    Function to calculate Differential Kramers Kronig transformation from
    experimental f2 (normalized XANES) to experimental f1

    Literature: Ohta and Ishida (1988) ... integration Methods for Kramers Kronig Transformation, Applied Spectroscopy 42,6
                Cross et al. (1998) ... theoretical x-ray resonant scattering amplitudes ...,Physical Review B, 58, 17

    datafile: Hephaestus *.f1f2 file (1eV steps) --> Cromer Liberman (CL) f1 and f2
    expfile: Athena normalized XANES file (same eV grid as .f1f2 file) --> experimental f2
    output: filename experimental E f1 f2 wll be written to (for use in rasd_menu)
    n = 30: number of datapoints to match f2 CL with f2 exp.
    e0: theoretical edge energy
    e0shift: to match theoretical edge with experimental edge
    """
    

    #f1f2 holds
    #[0] - E
    #[1] - f1 theo
    #[2] - f2 theo
    #[3] - diff f1theo f1exp 
    #[4] - f1 exp
    #[5] - f2 exp
    
    e0 = e0 + e0shift
    #read Cromer-Liberman calculated f1f2 (1eV grid) from HEPHAESTUS
    f = file(datafile, 'r')
    data = f.readlines()
    f.close()
    f1f2 = num.ndarray((0,6),float)
    for i in range(len(data)):
        if '#' not in data[i]:
            tmp = str.rsplit(data[i])
            if float(tmp[0])< e0:
                f1f2 = num.append(f1f2, [[int(round(float(tmp[0]))),float(tmp[1]),float(tmp[2]),0.,0.,0.]], axis = 0)
            else:
                f1f2 = num.append(f1f2, [[int(round(float(tmp[0]))),float(tmp[1]),float(tmp[2]),0.,0.,1.]], axis = 0)
    f1f2 = f1f2.transpose()
    f1f2[0] = f1f2[0] + e0shift

    #read experimental f2 (normalized XANES)
    f = file(expfile, 'r')
    data = f.readlines()
    f.close()
    f2exp = num.ndarray((0,2),float)
    for i in range(len(data)):
        if '#' not in data[i]:
            tmp = str.rsplit(data[i])
            f2exp = num.append(f2exp, [[int(round(float(tmp[0]))),float(tmp[1])]], axis = 0)
    f2exp = f2exp.transpose()

    #associate experimental values to calculated values
    i=0
    for i in range(len(f1f2[0])):
        j=0
        for j in range(len(f2exp[0])):
            if f1f2[0][i]== f2exp[0][j]:
                f1f2[5][i] = f2exp[1][j]


    lower = 0
    upper = 0
    i = 0
    for i in range(n):
        lower = lower + (f1f2[2][i] - f1f2[5][i])
        upper = upper + (f1f2[2][len(f1f2[0])-i-1]- f1f2[5][len(f1f2[0])-i-1]+ 1)
    
    lower = lower /n
    upper = upper /n
    f1f2[5] = f1f2[5]*(upper-lower) + lower

    #calculate f1exp from f2exp (Kramers-Kronig) (see Ohta 1988/ Cross 1998)
    for i in range(num.shape(f1f2)[1]):
        sum = 0
        if divmod(float(i),2)[1] == 0:
            j = 1
            for j in range(1, len(f1f2[0]),2):
                sum = sum + (f1f2[5][j]-f1f2[2][j])/(f1f2[0][j] - f1f2[0][i])+(f1f2[5][j]-f1f2[2][j])/(f1f2[0][j] + f1f2[0][i])
        else:
            j = 0
            for j in range(0, len(f1f2[0]),2):
                sum = sum + (f1f2[5][j]-f1f2[2][j])/(f1f2[0][j] - f1f2[0][i])+(f1f2[5][j]-f1f2[2][j])/(f1f2[0][j] + f1f2[0][i])
        f1f2[3][i] = (sum * 2 / num.pi)


    f1f2[4] = f1f2[1] + f1f2[3]

    #write experimental values to .f1f2 file
    f = file(output,'w')
    f.write('# file: "'+output+'" containing experimental f1 f2 values \n')
    f.write('# calculated using Cromer Liberman f1 f2 from: "'+datafile+'"\n')
    f.write('# and experimental f2 from: "'+expfile+'"\n')
    f.write('# E0 = '+str(e0)+', e0shift = '+str(e0shift)+'\n')
    f.write('# Energy f1exp f2exp \n')
    i=0
    for i in range(len(f1f2[0])):
        f.write(str(f1f2[0][i])+'   '+str(f1f2[4][i])+'   '+str(f1f2[5][i])+' \n')
    f.close()

    #plot results
    figure(1)
    clf()
    plot(f1f2[0],f1f2[1],'b-')
    plot(f1f2[0],f1f2[2],'b-')
    plot(f1f2[0],f1f2[4],'r-')
    plot(f1f2[0],f1f2[5],'g-')
