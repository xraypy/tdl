###################################################################################################
"""
Frank Heberling (Frank.Heberling@ine.fzk.de)
(based partly on matlab code)

crystallograhic calculation functions for the active area correction
mostly copied from matlab files
"""
######################################################################

import numpy as num
#from xtal.xtal_calcs_aa import *

######################################################################

def calc_Z(angles):
    #calculate goniometer rotation Matrix Z
    #angles is an array of diffrractometer angles as in spec file (del,eta,chi,phi,nu,mu)
    phi = num.radians(angles[3])
    chi = num.radians(angles[2])
    eta = num.radians(angles[1])
    mu  = num.radians(angles[5])
    Phi = num.array([[cos(phi),sin(phi),0],[-sin(phi),cos(phi),0],[0,0,1]],float)
    Chi = num.array([[cos(chi),0,sin(chi)],[0,1,0],[-sin(chi),0,cos(chi)]],float)
    Eta = num.array([[cos(eta),sin(eta),0],[-sin(eta),cos(eta),0],[0,0,1]],float)
    Mu  = num.array([[1,0,0],[0,cos(mu),-sin(mu)],[0,sin(mu),cos(mu)]],float)
    Z = num.dot(num.dot(num.dot(Mu,Eta),Chi),Phi)
    return Z

def vec_norm(a,G):
    #norm of  vector a in coordinatesystem with metric tensor G
    n = num.dot(a,num.dot(G,a))**0.5
    return n

def vec_angle(a,b,G):
    #calculate angle between vectors a and b in coordinate system with metric tensor G
    alpha = num.arccos((num.dot(a,num.dot(G,b)))/(vec_norm(a,G)*vec_norm(b,G)))
    return alpha


def calc_Q(angles, lam):
    #calculate k_in, k_out, and momentum transfer, Q, in lab frame, and two-theta
    #angles is an array of diffrractometer angles as in spec file (del,eta,chi,phi,nu,mu)
    k_in = 2* num.pi / lam * num.array([0,1,0])
    delta = num.radians(angles[0])
    nu = num.radians(angles[4])
    k_out = 2* num.pi / lam * num.array([sin(delta),cos(nu)*cos(delta),sin(nu)*cos(delta)])
    Q = k_out - k_in
    tth = vec_angle(k_in, k_out, num.eye(3))
    return k_in, k_out, Q, tth

def calc_G(cell):
    #calculate the metric Tensor, G, of the coordinate system (cell)
    #and the reziprocal metric Tensor, G_rez, and coordinate system (cell_rez)
    a = cell[0]
    b = cell[1]
    c = cell[2]
    alpha = num.radians(cell[3])
    beta  = num.radians(cell[4])
    gamma = num.radians(cell[5])
    G = num.zeros((3,3), float)
    G[0][0] = a**2
    G[0][1] = a*b*cos(gamma)
    G[0][2] = a*c*cos(beta)
    G[1][0] = G[0][1]
    G[1][1] = b**2
    G[1][2] = b*c*cos(alpha)
    G[2][0] = G[0][2]
    G[2][1] = G[1][2]
    G[2][2] = c**2
    G_rez = num.linalg.inv(G)
    a_rez = G_rez[0][0]**0.5
    b_rez = G_rez[1][1]**0.5
    c_rez = G_rez[2][2]**0.5
    alpha_rez = num.degrees(num.arccos(G_rez[1][2]/(b_rez * c_rez)))
    beta_rez  = num.degrees(num.arccos(G_rez[0][2]/(a_rez * c_rez)))
    gamma_rez = num.degrees(num.arccos(G_rez[0][1]/(a_rez * b_rez)))
    cell_rez = num.array([a_rez, b_rez, c_rez, alpha_rez, beta_rez, gamma_rez])
    return G, G_rez, cell_rez

def calc_B(cell):
    #calculate the B matrix
    #
    # the matrix B transforms the indicies of the 
    # reciprocal vector h to cartesian coordinates -->hc
    # where the cartesian system has been chosen so that:
    # x is parrallel to b1, 
    # y is in the plane of b1 and b2, 
    # z is perp to the plane of b1 and b2 
    # (see Busing and Levy)
    # Therefore:   
    #              h_c = B*h
    ##########################################################
    a3 = cell[2]
    A1 = num.radians(cell[3])
    G, G_rez, cell_rez = calc_G(cell)
    b1 = cell_rez[0]
    b2 = cell_rez[1]
    b3 = cell_rez[2]
    B1 = num.radians(cell_rez[3])
    B2 = num.radians(cell_rez[4])
    B3 = num.radians(cell_rez[5])
    B = num.array([[b1, b2*cos(B3), b3*cos(B2)],
                   [0, b2*sin(B3), -b3*sin(B2)*cos(A1)],
                   [0, 0, 1/a3]])
    return B

def calc_U(h1, h2, angles1, angles2, lambda1, lambda2, cell):
    #calculates the orientation matrix, U, from two reflectons h1 ,h2 and corresponding
    #diffractometer angles, lambdas and a given crystal coordinate system
    B = calc_B(cell) #calc B Matrix, carthesian crystal coords.
    #calc goniom. rot. Matrix and Momentum Transfer for both reflection
    Z1 = calc_Z(angles1)
    kin1, kout1, Q1, tth1 = calc_Q(angles1, lambda1)
    Z2 = calc_Z(angles2)
    kin2, kout2, Q2, tth2 = calc_Q(angles2, lambda2)
    #calc phi frame coords for diffraction vectors
    vphi_1 = num.dot(num.linalg.inv(Z1), (Q1/(2*num.pi)))
    vphi_2 = num.dot(num.linalg.inv(Z2), (Q2/(2*num.pi)))
    #calc carthesion coords of reflections
    hc1 = num.dot(B, h1)
    hc2 = num.dot(B, h2)
    #define tc vectors from hc vectors 
    tc1 = hc1 / vec_norm(hc1, num.eye(3))
    tc3 = num.cross(tc1, hc2) / vec_norm(num.cross(tc1, hc2), num.eye(3))
    tc2 = num.cross(tc3, tc1) / vec_norm(num.cross(tc3, tc1), num.eye(3))
    #define tphi vectors from vphi vectors
    tphi_1 = vphi_1 / vec_norm(vphi_1,num.eye(3))
    tphi_3 = num.cross(tphi_1,vphi_2) / vec_norm(num.cross(tphi_1,vphi_2), num.eye(3))
    tphi_2 = num.cross(tphi_3,tphi_1) / vec_norm(num.cross(tphi_3,tphi_1), num.eye(3))
    #define cooresponding matrices
    Tc = num.transpose(num.array([tc1,tc2,tc3]))
    Tphi = num.transpose(num.array([tphi_1,tphi_2,tphi_3]))
    #calc orientation matrix U and UB
    U = num.dot(Tphi, num.linalg.inv(Tc))
    UB = num.dot(U,B)
    return U, UB
                   
def area(A,B,C,G):
    #calculates the area of a triangle defined by cornerpoints A, B, and C
    #in coordinatesystem with metric tensor G
    a = (C-B)
    b = (C-A)
    tri = vec_norm(num.cross(a,b),G)/2
    return tri

def intersect(A,B,E,F):
    #intersects the two lines A+ab*(B-A) and E+ef*(F-E) in the surface plane (2D: (z=0))
    ef = ((E[0]-A[0])/(B[0]-A[0]) - (E[1]-A[1])/(B[1]-A[1])) / ((F[1]-E[1])/(B[1]-A[1]) - (F[0]-E[0])/(B[0]-A[0]))
    ab = (-A[0]+E[0]+ ef*(F[0]-E[0]))/(B[0]-A[0])
    if ef > 0 and ef < 1 and ab > 0 and ab < 1:
        intersection = True
        Point = A+ab*(B-A)
        return intersection, Point
    else:
        intersection = False
        Point = A+ab*(B-A)
        return intersection, Point

def calc_M(e0,UB,Z):
    #define a surface coord system as righthanded orthonormal system with zs along
    #azimuthal referenz vector and ys as projection of [0,-1,0] to the surface
    nm = num.dot(Z,num.dot(UB, e0)) #lab frame coords of the azimuthal ref. vector
    nm = nm / vec_norm(nm, num.eye(3))
    zs = nm
    v = num.array([0,-1,0])
    ys = v - (num.dot(v,nm) * nm) / vec_norm(nm, num.eye(3))**2 
    ys = ys / vec_norm(ys, num.eye(3))
    xs = num.cross(ys,zs) / vec_norm(num.cross(ys,zs), num.eye(3))
    #define the Matrix, M that transforms lab frame vectors to surface frame vectors
    M = num.array([xs,ys,zs])
    return M

######################################################################
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
    e0 = num.array(tmp.state['G'][3:6],float)
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
    M = calc_M(e0,UB,num.eye(3))
    #Sample corner points, surface frame at Z = I:
    phi_d1 = num.radians(phi_d1)
    phi_d2 = num.radians(phi_d2)
    E0 = num.array([cos(phi_d1),-sin(phi_d1),0])*d1/2 
    F0 = num.array([cos(phi_d2),-sin(phi_d2),0])*d2/2
    #sample Cornerpoints in lab frame at Z = I
    E0 = num.dot(num.linalg.inv(M),E0)
    F0 = num.dot(num.linalg.inv(M),F0)
    H0 = -F0
    if Debug:
        print 'Sample corner point E at zero pos. in lab frame: '+str(E0)
        print 'Sample corner point F at zero pos. in lab frame: '+str(F0)
    #Sample area
    S_area = 2* area(E0,F0,H0,num.eye(3))
    if Debug: print 'Sample area: '+str(S_area)
    #Beam corner points in plane perpendicular to beam (0,1,0)
    A0 = num.array([ Bheight/2,  0,-Bwidth/2],float)
    B0 = num.array([ Bheight/2,  0, Bwidth/2],float)
    D0 = -B0
    B_area = 2* area(A0,B0,D0,num.eye(3))
    #initialize correction arrays to which results are appended    
    spilloff_array = num.array([],float)
    footprint_array = num.array([],float)
    #loop through scan named 'tmp'
    i=0
    for i in range(tmp.dims[0]):
        #read diffractometer angles and calc Goniometer rotation matrix Z 
        eta = tmp.scalers['eta'][i]
        mu  = tmp.scalers['mu'][i]
        chi = tmp.positioners['chi'][i]
        phi = tmp.positioners['phi'][i]
        angles = num.array([0,eta,chi,phi,0,mu])
        Z = calc_Z(angles)
        #calc M matrix at this 
        M = calc_M(e0, UB, Z)
        #calc surface normal in lab frame
        n = num.dot(num.linalg.inv(M),num.array([0,0,1]))
        n = n / vec_norm(n, num.eye(3))
        #calc alpha
        alpha = num.degrees(num.arcsin(num.dot(n,([0,-1,0]))))
        if alpha < 0: print '"Warning Alpha is negative!!!"'
        #Beam corner points in lab frame
        A0[1]=(-n[0]*A0[0]-n[2]*A0[2])/n[1]
        B0[1]=(-n[0]*B0[0]-n[2]*B0[2])/n[1]
        #Beam corner points in surface frame
        A = num.dot(num.linalg.inv(M),A0)
        B = num.dot(num.linalg.inv(M),B0)
        #make sure A and B are on the surface
        A[2] = 0.0
        B[2] = 0.0
        C = -A
        D = -B
        #Beam area
        B_area = 2* area(A,B,D,num.eye(3))
        #Sample corner points in surface frame
        E = num.dot(M,num.dot(Z,E0))
        F = num.dot(M,num.dot(Z,F0))
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
            footprint = 2*area(EFAB, EFDC, GHAB,num.eye(3))
        elif efad and efbc and ghad and ghbc:#2
            condition = 2
            footprint = 2*area(EFAD,GHAD,EFBC,num.eye(3))
        elif fgad and fgbc and ehad and ehbc:#3
            condition = 3
            footprint = 2*area(FGAD,FGBC,EHAD,num.eye(3))
        elif fgab and fgdc and ehab and ehdc :#4
            condition = 4
            footprint = 2* area(FGAB,FGDC,EHAB,num.eye(3))
        elif efab and fgdc and not efad:#5
            footprint = 2*(area(EFAB,F,FGDC,num.eye(3))+area(EFAB,FGDC,EHAB,num.eye(3)))
            condition = 5
        elif ehad and efad and not efdc:#6
            condition = 6
            footprint = 2* (area(EHAD,EFAD,F,num.eye(3))+area(F,H,EHAD,num.eye(3)))
        elif ghab and ehdc and not ehad:#7
            condition = 7
            footprint = 2*(area(GHAB,H,FGAB,num.eye(3)) + area(FGAB,H,F,num.eye(3)))
        elif fgad and ghad and not ghdc:#8
            condition = 8
            footprint = 2*(area(F,FGAD,H,num.eye(3))+area(FGAD,GHAD,H,num.eye(3)))
        elif fgab and ghdc and not fgad:#9
            condition = 9
            footprint = 2* (area(FGAB,G,GHDC,num.eye(3)) + area(FGAB,GHDC,EFAB,num.eye(3)))
        elif efad and fgad and not fgdc:#10
            condition = 10
            footprint = 2*(area(E,EFAD,FGAD,num.eye(3))+area(FGAD,G,E,num.eye(3)))
        elif efdc and fgdc and not fgbc:#11
            condition = 11
            footprint = 2*(area(EHAB,GHAB,E,num.eye(3))+area(E,G,GHAB,num.eye(3)))
        elif efbc and fgbc and not efdc:#12
            condition = 12
            footprint = 2*(area(E,EFBC,FGBC,num.eye(3))+area(E,G,FGBC,num.eye(3)))
        elif efab and efad and fgdc:#13
            condition = 13
            footprint = 2* (area(EFAB,EFAD,FGDC,num.eye(3))+area(EFAD,FGAD,FGDC,num.eye(3))+area(FGDC,EHAB,EFAB,num.eye(3)))
        elif ehab and ehad and efad:#14
            condition = 14
            footprint = 2*(area(EHAD,EHAB,GHAB,num.eye(3))+area(EHAB,FGBC,GHAB,num.eye(3))+area(GHAB,GHBC,FGBC,num.eye(3)))
        elif ghab and ghad and ehad:#15
            condition = 15
            footprint = 2*(area(EFBC,GHAD,GHAB,num.eye(3))+area(EFBC,GHAB,FGAB,num.eye(3))+area(EFBC,FGBC,FGAB,num.eye(3)))
        elif fgab and fgad and ghad:#16
            condition = 16
            footprint = 2*(area(EHBC,EFBC,EFAB,num.eye(3))+area(EHBC,EFAB,FGAB,num.eye(3))+area(EHBC,FGAB,FGAD,num.eye(3)))
        elif not efab and not efad and not fgad and not fgab:#17
            if S_area > B_area:
                condition = 17.1
                footprint = B_area
            else:
                condition = 17.2
                footprint = S_area
        elif fgad and fgdc and not efad:#18
            condition = 18
            footprint = 2* (area(A,FGAD,FGDC,num.eye(3))+area(A,FGDC,EHAB,num.eye(3)))
        elif efad and efdc and not fgdc:#19
            condition = 19
            footprint = 2*(area(EFAD,EFDC,C,num.eye(3))+area(EFAD,A,C,num.eye(3)))
        elif ehad and ehdc and not efdc:#20
            condition = 20
            footprint = 2*(area(EHDC,EHAD,A,num.eye(3))+area(EHDC,A,C,num.eye(3)))
        elif ghad and ghdc and not ehdc:#21
            condition = 21
            footprint = 2*(area(GHAD,GHDC,A,num.eye(3))+area(GHDC,A,C,num.eye(3)))
        elif efab and efad and not fgad:#22
            condition = 22
            footprint = 2*(area(EFAB,EFAD,B,num.eye(3))+area(EFAD,B,D,num.eye(3)))
        elif ehab and ehad and not efad:#23
            condition = 23
            footprint = 2*(area(EHAB,EHAD,B,num.eye(3))+area(EHAD,B,D,num.eye(3)))
        elif ghab and ghad and not ehad:#24
            condition = 24
            footprint = 2*(area(GHAB,GHAD,B,num.eye(3))+area(GHAD,D,B,num.eye(3)))
        elif fgab and fgad and not ghad:#25
            condition = 25
            footprint = 2*(area(FGAB,FGAD,B,num.eye(3))+area(FGAD,D,B,num.eye(3)))
        
        else: footprint = -1

        if footprint < B_area:
            spilloff = footprint/B_area
        else: spilloff = 1

        footprint_array = num.append(footprint_array, footprint)
        spilloff_array = num.append(spilloff_array, spilloff)

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


