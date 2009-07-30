###################################################################################################
"""
Frank Heberling (Frank.Heberling@ine.fzk.de)
(based partly on matlab code)
todo:
    test
"""
###################################################################################################

import numpy as Num
from xtal.xtal_calcs_aa import *


def active_area(tmp,d1,d2,phi_d1,phi_d2,Bwidth,Bheight, Debug ):
    """
    active_area(d1,d2,phi_d1,phi_d2,Bwidth,Bheight, Debug )
    
    calculates footprint and spilloff for surface diffraction measurements for a
    scandata_object 'tmp'
    assumption is that sample is parallelogram shaped and perfectly centered
    Sample orientation needs to be defined in polar coordinates
    by phi angles at which sample diagonals are vertically aligned (phi_d1 and phi_d2) at eta=chi=0
    and length of corresponding crystal diagonals d1, d2(mm)
    Bwidth and Bheight define Beam proportions (mm) 
    Debug flag, True - to print intermediate results
    """
    ##########################################################################################################################
    #azimuthal reference vector, e0 (hkl) from spec file
    e0 = Num.array(tmp.state['G'][3:6],float)
    #lattice params and alignment params from spec file
    cell = tmp.state['G'][22:28]
    h1 = tmp.state['G'][34:37]
    h2 = tmp.state['G'][37:40]
    angle1 = tmp.state['G'][40:46]
    angle2 = tmp.state['G'][46:52]
    lambda1 = tmp.state['G'][52]
    lambda2 = tmp.state['G'][53]
    #calculate orientation matrix
    U, UB = calc_U(h1, h2, angle1, angle2, lambda1, lambda2, cell)
    #M matrix transforms lab to surface frame coordinates
    M = calc_M(e0,UB,Num.eye(3))
    #Sample corner points, surface frame at Z = I:
    phi_d1 = Num.radians(phi_d1)
    phi_d2 = Num.radians(phi_d2)
    E0 = Num.array([cos(phi_d1),-sin(phi_d1),0])*d1/2 
    F0 = Num.array([cos(phi_d2),-sin(phi_d2),0])*d2/2
    #sample Cornerpoints in lab frame at Z = I
    E0 = Num.dot(Num.linalg.inv(M),E0)
    F0 = Num.dot(Num.linalg.inv(M),F0)
    H0 = -F0
    if Debug:
        print 'Sample corner point E at zero pos. in lab frame: '+str(E0)
        print 'Sample corner point F at zero pos. in lab frame: '+str(F0)
    #Sample area
    S_area = 2* area(E0,F0,H0,Num.eye(3))
    if Debug: print 'Sample area: '+str(S_area)
    #Beam corner points in plane perpendicular to beam (0,1,0)
    A0 = Num.array([ Bheight/2,  0,-Bwidth/2],float)
    B0 = Num.array([ Bheight/2,  0, Bwidth/2],float)
    D0 = -B0
    B_area = 2* area(A0,B0,D0,Num.eye(3))
    #initialize correction arrays to which results are appended    
    spilloff_array = Num.array([],float)
    footprint_array = Num.array([],float)
    #loop through scan named 'tmp'
    i=0
    for i in range(tmp.dims[0]):
        #read diffractometer angles and calc Goniometer rotation matrix Z 
        eta = tmp.scalers['eta'][i]
        mu  = tmp.scalers['mu'][i]
        chi = tmp.positioners['chi'][i]
        phi = tmp.positioners['phi'][i]
        angles = Num.array([0,eta,chi,phi,0,mu])
        Z = calc_Z(angles)
        #calc M matrix at this 
        M = calc_M(e0, UB, Z)
        #calc surface normal in lab frame
        n = Num.dot(Num.linalg.inv(M),Num.array([0,0,1]))
        n = n / vec_norm(n, Num.eye(3))
        #calc alpha
        alpha = Num.degrees(Num.arcsin(Num.dot(n,([0,-1,0]))))
        if alpha < 0: print '"Warning Alpha is negative!!!"'
        #Beam corner points in lab frame
        A0[1]=(-n[0]*A0[0]-n[2]*A0[2])/n[1]
        B0[1]=(-n[0]*B0[0]-n[2]*B0[2])/n[1]
        #Beam corner points in surface frame
        A = Num.dot(Num.linalg.inv(M),A0)
        B = Num.dot(Num.linalg.inv(M),B0)
        #make sure A and B are on the surface
        A[2] = 0.0
        B[2] = 0.0
        C = -A
        D = -B
        #Beam area
        B_area = 2* area(A,B,D,Num.eye(3))
        #Sample corner points in surface frame
        E = Num.dot(M,Num.dot(Z,E0))
        F = Num.dot(M,Num.dot(Z,F0))
        #make sure E and F are on the surface
        E[2] = 0.0
        F[2] = 0.0
        G = -E
        H = -F
        #intersect all the sample edges with all the beam edges
        efab, EFAB = intersect(E,F,A,B)
        efdc, EFDC = intersect(E,F,D,C)
        efad, EFAD = intersect(E,F,A,D)
        efbc, EFBC = intersect(E,F,B,C)

        ghab, GHAB = intersect(G,H,A,B)
        ghdc, GHDC = intersect(G,H,D,C)
        ghad, GHAD = intersect(G,H,A,D)
        ghbc, GHBC = intersect(G,H,B,C)

        fgad, FGAD = intersect(F,G,A,D)
        fgdc, FGDC = intersect(F,G,D,C)
        fgab, FGAB = intersect(F,G,A,B)
        fgbc, FGBC = intersect(F,G,B,C)

        ehab, EHAB = intersect(E,H,A,B)
        ehbc, EHBC = intersect(E,H,B,C)
        ehdc, EHDC = intersect(E,H,D,C)
        ehad, EHAD = intersect(E,H,A,D)
        #25 different possibilities how beam and sample can overlap (see footprint.pdf)
        if efab and efdc and ghab and ghdc:#1
            condition = 1
            footprint = 2*area(EFAB, EFDC, GHAB,Num.eye(3))
        elif efad and efbc and ghad and ghbc:#2
            condition = 2
            footprint = 2*area(EFAD,GHAD,EFBC,Num.eye(3))
        elif fgad and fgbc and ehad and ehbc:#3
            condition = 3
            footprint = 2*area(FGAD,FGBC,EHAD,Num.eye(3))
        elif fgab and fgdc and ehab and ehdc :#4
            condition = 4
            footprint = 2* area(FGAB,FGDC,EHAB,Num.eye(3))
        elif efab and fgdc and not efad:#5
            footprint = 2*(area(EFAB,F,FGDC,Num.eye(3))+area(EFAB,FGDC,EHAB,Num.eye(3)))
            condition = 5
        elif ehad and efad and not efdc:#6
            condition = 6
            footprint = 2* (area(EHAD,EFAD,F,Num.eye(3))+area(F,H,EHAD,Num.eye(3)))
        elif ghab and ehdc and not ehad:#7
            condition = 7
            footprint = 2*(area(GHAB,H,FGAB,Num.eye(3)) + area(FGAB,H,F,Num.eye(3)))
        elif fgad and ghad and not ghdc:#8
            condition = 8
            footprint = 2*(area(F,FGAD,H,Num.eye(3))+area(FGAD,GHAD,H,Num.eye(3)))
        elif fgab and ghdc and not fgad:#9
            condition = 9
            footprint = 2* (area(FGAB,G,GHDC,Num.eye(3)) + area(FGAB,GHDC,EFAB,Num.eye(3)))
        elif efad and fgad and not fgdc:#10
            condition = 10
            footprint = 2*(area(E,EFAD,FGAD,Num.eye(3))+area(FGAD,G,E,Num.eye(3)))
        elif efdc and fgdc and not fgbc:#11
            condition = 11
            footprint = 2*(area(EHAB,GHAB,E,Num.eye(3))+area(E,G,GHAB,Num.eye(3)))
        elif efbc and fgbc and not efdc:#12
            condition = 12
            footprint = 2*(area(E,EFBC,FGBC,Num.eye(3))+area(E,G,FGBC,Num.eye(3)))
        elif efab and efad and fgdc:#13
            condition = 13
            footprint = 2* (area(EFAB,EFAD,FGDC,Num.eye(3))+area(EFAD,FGAD,FGDC,Num.eye(3))+area(FGDC,EHAB,EFAB,Num.eye(3)))
        elif ehab and ehad and efad:#14
            condition = 14
            footprint = 2*(area(EHAD,EHAB,GHAB,Num.eye(3))+area(EHAB,FGBC,GHAB,Num.eye(3))+area(GHAB,GHBC,FGBC,Num.eye(3)))
        elif ghab and ghad and ehad:#15
            condition = 15
            footprint = 2*(area(EFBC,GHAD,GHAB,Num.eye(3))+area(EFBC,GHAB,FGAB,Num.eye(3))+area(EFBC,FGBC,FGAB,Num.eye(3)))
        elif fgab and fgad and ghad:#16
            condition = 16
            footprint = 2*(area(EHBC,EFBC,EFAB,Num.eye(3))+area(EHBC,EFAB,FGAB,Num.eye(3))+area(EHBC,FGAB,FGAD,Num.eye(3)))
        elif not efab and not efad and not fgad and not fgab:#17
            if S_area > B_area:
                condition = 17.1
                footprint = B_area
            else:
                condition = 17.2
                footprint = S_area
        elif fgad and fgdc and not efad:#18
            condition = 18
            footprint = 2* (area(A,FGAD,FGDC,Num.eye(3))+area(A,FGDC,EHAB,Num.eye(3)))
        elif efad and efdc and not fgdc:#19
            condition = 19
            footprint = 2*(area(EFAD,EFDC,C,Num.eye(3))+area(EFAD,A,C,Num.eye(3)))
        elif ehad and ehdc and not efdc:#20
            condition = 20
            footprint = 2*(area(EHDC,EHAD,A,Num.eye(3))+area(EHDC,A,C,Num.eye(3)))
        elif ghad and ghdc and not ehdc:#21
            condition = 21
            footprint = 2*(area(GHAD,GHDC,A,Num.eye(3))+area(GHDC,A,C,Num.eye(3)))
        elif efab and efad and not fgad:#22
            condition = 22
            footprint = 2*(area(EFAB,EFAD,B,Num.eye(3))+area(EFAD,B,D,Num.eye(3)))
        elif ehab and ehad and not efad:#23
            condition = 23
            footprint = 2*(area(EHAB,EHAD,B,Num.eye(3))+area(EHAD,B,D,Num.eye(3)))
        elif ghab and ghad and not ehad:#24
            condition = 24
            footprint = 2*(area(GHAB,GHAD,B,Num.eye(3))+area(GHAD,D,B,Num.eye(3)))
        elif fgab and fgad and not ghad:#25
            condition = 25
            footprint = 2*(area(FGAB,FGAD,B,Num.eye(3))+area(FGAD,D,B,Num.eye(3)))
        
        else: footprint = -1

        if footprint < B_area:
            spilloff = footprint/B_area
        else: spilloff = 1

        footprint_array = Num.append(footprint_array, footprint)
        spilloff_array = Num.append(spilloff_array, spilloff)

        if Debug:
            print '\n'
            print 'L: '+str(tmp.scalers['L'][i])
            print 'surface normal in lab frame: '+str(n)
            print 'Beam corner point A in surface frame: '+str(A)
            print 'Beam corner point B in surface frame: '+str(B)
            print 'Sample corner point E in surface frame: '+str(E)
            print 'Sample corner point F in surface frame: '+str(F)
            print 'Beam area: '+str(B_area)
            print 'intersection condition (cp. footprint.pdf): '+str(condition)
            print 'Incident angle, Alpha: '+str(alpha)
            print 'Alpha from spec file: '+str(tmp.scalers['Alpha'][i])
            print 'footprint: '+str(footprint)
            print 'spilloff: '+str(spilloff)
            
    
    return footprint_array, spilloff_array


