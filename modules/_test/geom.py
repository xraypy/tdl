import numpy as Num

def sample_rotation(v,phi,chi,eta,mu):
    #use goniometer rotations on vector v:
    Phi = Num.array([[cos(phi),sin(phi),0],[-sin(phi),cos(phi),0],[0,0,1]],float)
    Chi = Num.array([[cos(chi),0,sin(chi)],[0,1,0],[-sin(chi),0,cos(chi)]],float)
    Eta = Num.array([[cos(eta),sin(eta),0],[-sin(eta),cos(eta),0],[0,0,1]],float)
    Mu  = Num.array([[1,0,0],[0,cos(mu),-sin(mu)],[0,sin(mu),cos(mu)]],float)
    v = Num.dot(Phi,v)
    v = Num.dot(Chi,v)
    v = Num.dot(Eta,v)
    v = Num.dot(Mu,v)   
    return v

def distance(A,B):
    n = ((A[0]-B[0])**2+(A[1]-B[1])**2+(A[2]-B[2])**2)**0.5
    return n

def area(A,B,C):
    #calculates the area of a triangle defined by cornerpoints A, B, and C
    a = distance(A,B)
    b = distance(B,C)
    c = distance(C,A)
    tri =(((a**2+b**2+c**2)**2-2*(a**4+b**4+c**4))/16)**0.5
    return tri

def intersect(A,B,E,F):
    #intersects the two lines A+ab*(B-A) and E+ef*(F-E)
    if B[0]-A[0]!=0 and B[1]-A[1]!=0:
        q0 = 0
        q1 = 1
    elif B[0]-A[0]!=0 and B[2]-A[2]!=0:
        q0=0
        q1=2
    elif B[1]-A[1]!=0 and B[2]-A[2]!=0:
        q0 = 1
        q1 = 2
        
    ef = ((E[q0]-A[q0])/(B[q0]-A[q0]) - (E[q1]-A[q1])/(B[q1]-A[q1])) / ((F[q1]-E[q1])/(B[q1]-A[q1]) - (F[q0]-E[q0])/(B[q0]-A[q0]))
    ab = (-A[q0]+E[q0]+ ef*(F[q0]-E[q0]))/(B[q0]-A[q0])
    if ef > 0 and ef < 1 and ab > 0 and ab < 1:
        intersection = True
        Point = A+ab*(B-A)
        return intersection, Point
    else:
        intersection = False
        Point = A+ab*(B-A)
        return intersection, Point

######################################################################################################################################
#geom.py
#geometrical corrections for surface diffraction measurements
#needed input : scandata_object 'tmp'
#assumption is that sample is parallelogram shaped and perfectly centered

#{These values are special for our measurements

#Sample corner points defined in polar coordinates
#by phi angle at which the corner is above the center and length of corresponding crystal diagonal (mm)
d1 = 11.183
phi_d1 = 56.4
d2 = 13.69
phi_d2 = -17
#Measured Beam proportions as vertical and horizontal width (mm)
Vw = 1.25
Hw = 1.25
#}
Debug = True #Debug flag to print intermediate results
#######################################################################################################################################
#zero position of surface normal, e0 (diffractometer coordinates (You'99)) from spec file
e0 = Num.array(tmp.state['G'][3:6],float)
# normalize e0
e0 = e0/distance(([0,0,0]),e0)
e0 = e0.transpose()
if Debug: print e0
#Sample corner points in surface zero position:
phi_d1 = Num.radians(phi_d1)
phi_d2 = Num.radians(phi_d2)
E0 = Num.array([cos(phi_d1),-sin(phi_d1),(-e0[0]*cos(phi_d1)+e0[1]*sin(phi_d1))/e0[2]])*d1/2
F0 = Num.array([cos(phi_d2),-sin(phi_d2),(-e0[0]*cos(phi_d2)+e0[1]*sin(phi_d2))/e0[2]])*d2/2
E0 = E0.transpose()
F0 = F0.transpose()
H0 = -F0
if Debug:
    print E0
    print F0
#Sample area
S_area = 2* area(E0,F0,H0)
#Beam corner points at plane perpendicular to beam (0,-1,0)
A = Num.array([ Vw/2,  0,-Hw/2],float)
B = Num.array([ Vw/2,  0, Hw/2],float)
#initialize correction arrays to whch results are appended    
spilloff_array = Num.array([],float)
footprint_array = Num.array([],float)
#loop through tmp
i=0
for i in range(tmp.dims[0]):
    #read diffractometer angles(!!!angle conventions for 13BM unknown!!! --> unsure about +-) 
    eta = -Num.radians(tmp.scalers['eta'][i])
    mu  = -Num.radians(tmp.scalers['mu'][i])
    chi = Num.radians(tmp.positioners['chi'][i])
    phi = -Num.radians(tmp.positioners['phi'][i])
    #Rotate surface normal
    e = sample_rotation(e0,phi,chi,eta,mu)
    #Rotate Sample corner points 
    E = sample_rotation(E0,phi,chi,eta,mu)
    F = sample_rotation(F0,phi,chi,eta,mu)
    E = E.transpose()
    F = F.transpose()
    G = -E
    H = -F
    #Beam corner points in sample surface plane
    A[1]=(-e[0]*A[0]-e[2]*A[2])/e[1]
    B[1]=(-e[0]*B[0]-e[2]*B[2])/e[1]
    C = -A
    D = -B
    #Beam area
    B_area = 2* area(A,B,D)
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
        footprint = 2*area(EFAB, EFDC, GHAB)
    elif efad and efbc and ghad and ghbc:#2
        condition = 2
        footprint = 2*area(EFAD,GHAD,EFBC)
    elif fgad and fgbc and ehad and ehbc:#3
        condition = 3
        footprint = 2*area(FGAD,FGBC,EHAD)
    elif fgab and fgdc and ehab and ehdc :#4
        condition = 4
        footprint = 2* area(FGAB,FGDC,EHAB)
    elif efab and fgdc and not efad:#5
        footprint = 2*(area(EFAB,F,FGDC)+area(EFAB,FGDC,EHAB))
        condition = 5
    elif ehad and efad and not efdc:#6
        condition = 6
        footprint = 2* (area(EHAD,EFAD,F)+area(F,H,EHAD))
    elif ghab and ehdc and not ehad:#7
        condition = 7
        footprint = 2*(area(GHAB,H,FGAB) + area(FGAB,H,F))
    elif fgad and ghad and not ghdc:#8
        condition = 8
        footprint = 2*(area(F,FGAD,H)+area(FGAD,GHAD,H))
    elif fgab and ghdc and not fgad:#9
        condition = 9
        footprint = 2* (area(FGAB,G,GHDC) + area(FGAB,GHDC,EFAB))
    elif efad and fgad and not fgdc:#10
        condition = 10
        footprint = 2*(area(E,EFAD,FGAD)+area(FGAD,G,E))
    elif efdc and fgdc and not fgbc:#11
        condition = 11
        footprint = 2*(area(EHAB,GHAB,E)+area(E,G,GHAB))
    elif efbc and fgbc and not efdc:#12
        condition = 12
        footprint = 2*(area(E,EFBC,FGBC)+area(E,G,FGBC))
    elif efab and efad and fgdc:#13
        condition = 13
        footprint = 2* (area(EFAB,EFAD,FGDC)+area(EFAD,FGAD,FGDC)+area(FGDC,EHAB,EFAB))
    elif ehab and ehad and efad:#14
        condition = 14
        footprint = 2*(area(EHAD,EHAB,GHAB)+area(EHAB,FGBC,GHAB)+area(GHAB,GHBC,FGBC))
    elif ghab and ghad and ehad:#15
        condition = 15
        footprint = 2*(area(EFBC,GHAD,GHAB)+area(EFBC,GHAB,FGAB)+area(EFBC,FGBC,FGAB))
    elif fgab and fgad and ghad:#16
        condition = 16
        footprint = 2*(area(EHBC,EFBC,EFAB)+area(EHBC,EFAB,FGAB)+area(EHBC,FGAB,FGAD))
    elif not efab and not efad and not fgad and not fgab:#17
        if S_area > B_area:
            condition = 17.1
            footprint = B_area
        else:
            condition = 17.2
            footprint = S_area
    elif fgad and fgdc and not efad:#18
        condition = 18
        footprint = 2* (area(A,FGAD,FGDC)+area(A,FGDC,EHAB))
    elif efad and efdc and not fgdc:#19
        condition = 19
        footprint = 2*(area(EFAD,EFDC,C)+area(EFAD,A,C))
    elif ehad and ehdc and not efdc:#20
        condition = 20
        footprint = 2*(area(EHDC,EHAD,A)+area(EHDC,A,C))
    elif ghad and ghdc and not ehdc:#21
        condition = 21
        footprint = 2*(area(GHAD,GHDC,A)+area(GHDC,A,C))
    elif efab and efad and not fgad:#22
        condition = 22
        footprint = 2*(area(EFAB,EFAD,B)+area(EFAD,B,D))
    elif ehab and ehad and not efad:#23
        condition = 23
        footprint = 2*(area(EHAB,EHAD,B)+area(EHAD,B,D))
    elif ghab and ghad and not ehad:#24
        condition = 24
        footprint = 2*(area(GHAB,GHAD,B)+area(GHAD,D,B))
    elif fgab and fgad and not ghad:#25
        condition = 25
        footprint = 2*(area(FGAB,FGAD,B)+area(FGAD,D,B))
        
    else: footprint = -1

    if footprint < B_area:
        spilloff = footprint/B_area
    else: spilloff = 1

    footprint_array = Num.append(footprint_array, footprint)
    spilloff_array = Num.append(spilloff_array, spilloff)

    if Debug:
        print e
        print condition
        alpha = Num.degrees(Num.arccos(Num.dot(e,([0,-1,0]))))-90
        naz = 90-Num.degrees(Num.arccos(Num.dot(e,([1,0,0]))))
        print naz
        print alpha #-tmp.scalers['Alpha'][i]
    
print footprint_array
print spilloff_array
