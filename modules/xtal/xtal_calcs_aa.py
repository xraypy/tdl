######################################################################
"""
Frank Heberling (Frank.Heberling@ine.fzk.de)

crystallograhic calculation functions for the active area correction
(see scandata.active_area.py)
mostly copied from matlab files
"""
######################################################################
import numpy as Num

def calc_Z(angles):
    #calculate goniometer rotation Matrix Z
    #angles is an array of diffrractometer angles as in spec file (del,eta,chi,phi,nu,mu)
    phi = Num.radians(angles[3])
    chi = Num.radians(angles[2])
    eta = Num.radians(angles[1])
    mu  = Num.radians(angles[5])
    Phi = Num.array([[cos(phi),sin(phi),0],[-sin(phi),cos(phi),0],[0,0,1]],float)
    Chi = Num.array([[cos(chi),0,sin(chi)],[0,1,0],[-sin(chi),0,cos(chi)]],float)
    Eta = Num.array([[cos(eta),sin(eta),0],[-sin(eta),cos(eta),0],[0,0,1]],float)
    Mu  = Num.array([[1,0,0],[0,cos(mu),-sin(mu)],[0,sin(mu),cos(mu)]],float)
    Z = Num.dot(Num.dot(Num.dot(Mu,Eta),Chi),Phi)
    return Z

def vec_norm(a,G):
    #norm of  vector a in coordinatesystem with metric tensor G
    n = Num.dot(a,Num.dot(G,a))**0.5
    return n

def vec_angle(a,b,G):
    #calculate angle between vectors a and b in coordinate system with metric tensor G
    alpha = Num.arccos((Num.dot(a,Num.dot(G,b)))/(vec_norm(a,G)*vec_norm(b,G)))
    return alpha


def calc_Q(angles, lam):
    #calculate k_in, k_out, and momentum transfer, Q, in lab frame, and two-theta
    #angles is an array of diffrractometer angles as in spec file (del,eta,chi,phi,nu,mu)
    k_in = 2* Num.pi / lam * Num.array([0,1,0])
    delta = Num.radians(angles[0])
    nu = Num.radians(angles[4])
    k_out = 2* Num.pi / lam * Num.array([sin(delta),cos(nu)*cos(delta),sin(nu)*cos(delta)])
    Q = k_out - k_in
    tth = vec_angle(k_in, k_out, Num.eye(3))
    return k_in, k_out, Q, tth

def calc_G(cell):
    #calculate the metric Tensor, G, of the coordinate system (cell)
    #and the reziprocal metric Tensor, G_rez, and coordinate system (cell_rez)
    a = cell[0]
    b = cell[1]
    c = cell[2]
    alpha = Num.radians(cell[3])
    beta  = Num.radians(cell[4])
    gamma = Num.radians(cell[5])
    G = Num.zeros((3,3), float)
    G[0][0] = a**2
    G[0][1] = a*b*cos(gamma)
    G[0][2] = a*c*cos(beta)
    G[1][0] = G[0][1]
    G[1][1] = b**2
    G[1][2] = b*c*cos(alpha)
    G[2][0] = G[0][2]
    G[2][1] = G[1][2]
    G[2][2] = c**2
    G_rez = Num.linalg.inv(G)
    a_rez = G_rez[0][0]**0.5
    b_rez = G_rez[1][1]**0.5
    c_rez = G_rez[2][2]**0.5
    alpha_rez = Num.degrees(Num.arccos(G_rez[1][2]/(b_rez * c_rez)))
    beta_rez  = Num.degrees(Num.arccos(G_rez[0][2]/(a_rez * c_rez)))
    gamma_rez = Num.degrees(Num.arccos(G_rez[0][1]/(a_rez * b_rez)))
    cell_rez = Num.array([a_rez, b_rez, c_rez, alpha_rez, beta_rez, gamma_rez])
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
    A1 = Num.radians(cell[3])
    G, G_rez, cell_rez = calc_G(cell)
    b1 = cell_rez[0]
    b2 = cell_rez[1]
    b3 = cell_rez[2]
    B1 = Num.radians(cell_rez[3])
    B2 = Num.radians(cell_rez[4])
    B3 = Num.radians(cell_rez[5])
    B = Num.array([[b1, b2*cos(B3), b3*cos(B2)],[0, b2*sin(B3), -b3*sin(B2)*cos(A1)],[0, 0, 1/a3]])
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
    vphi_1 = Num.dot(Num.linalg.inv(Z1), (Q1/(2*Num.pi)))
    vphi_2 = Num.dot(Num.linalg.inv(Z2), (Q2/(2*Num.pi)))
    #calc carthesion coords of reflections
    hc1 = Num.dot(B, h1)
    hc2 = Num.dot(B, h2)
    #define tc vectors from hc vectors 
    tc1 = hc1 / vec_norm(hc1, Num.eye(3))
    tc3 = Num.cross(tc1, hc2) / vec_norm(Num.cross(tc1, hc2), Num.eye(3))
    tc2 = Num.cross(tc3, tc1) / vec_norm(Num.cross(tc3, tc1), Num.eye(3))
    #define tphi vectors from vphi vectors
    tphi_1 = vphi_1 / vec_norm(vphi_1,Num.eye(3))
    tphi_3 = Num.cross(tphi_1,vphi_2) / vec_norm(Num.cross(tphi_1,vphi_2), Num.eye(3))
    tphi_2 = Num.cross(tphi_3,tphi_1) / vec_norm(Num.cross(tphi_3,tphi_1), Num.eye(3))
    #define cooresponding matrices
    Tc = Num.transpose(Num.array([tc1,tc2,tc3]))
    Tphi = Num.transpose(Num.array([tphi_1,tphi_2,tphi_3]))
    #calc orientation matrix U and UB
    U = Num.dot(Tphi, Num.linalg.inv(Tc))
    UB = Num.dot(U,B)
    return U, UB
                   
def area(A,B,C,G):
    #calculates the area of a triangle defined by cornerpoints A, B, and C
    #in coordinatesystem with metric tensor G
    a = (C-B)
    b = (C-A)
    tri = vec_norm(Num.cross(a,b),G)/2
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
    nm = Num.dot(Z,Num.dot(UB, e0)) #lab frame coords of the azimuthal ref. vector
    nm = nm / vec_norm(nm, Num.eye(3))
    zs = nm
    v = Num.array([0,-1,0])
    ys = v - (Num.dot(v,nm) * nm) / vec_norm(nm, Num.eye(3))**2 
    ys = ys / vec_norm(ys, Num.eye(3))
    xs = Num.cross(ys,zs) / vec_norm(Num.cross(ys,zs), Num.eye(3))
    #define the Matrix, M that transforms lab frame vectors to surface frame vectors
    M = Num.array([xs,ys,zs])
    return M
