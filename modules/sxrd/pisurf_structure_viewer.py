"""
Functions and classes for 3D structure viewing in pisurf

Authors/modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""

try:
    import vtk
except ImportError:
    pass
    
import numpy as Num
from tdl.modules.sxrd.atom_styles import atom_styles
from tdl.modules.sxrd.ctrfitcalcs import RB_update, param_unfold

######################################################################################################

def calcM(cell):
    alpha = Num.radians(cell[3])
    beta = Num.radians(cell[4])
    gamma = Num.radians(cell[5])
    x = Num.cos(beta)*cell[2]
    y = Num.cos(alpha)*cell[2]*Num.cos(Num.pi/2-gamma)
    z = (cell[2]**2-x**2-y**2)**0.5
    M = Num.array([[cell[0],0,0],\
                   [Num.cos(gamma)*cell[1],Num.cos(Num.pi/2-gamma)*cell[1],0],\
                   [x,y,z]],float)
    return M

def aff2euk(M, aff_points):
    euk_points = []
    for point in aff_points:
        euk_points.append( Num.dot(point,M) )
    return euk_points

def Ellipsoid(center, color, U):
    # create the quadric function definition
    a0 = U[0][0]
    a1 = U[1][1]
    a2 = U[2][2]
    a3 = 2*U[0][1]
    a4 = 2*U[1][2]
    a5 = 2*U[0][2]
    quadric = vtk.vtkQuadric()
    quadric.SetCoefficients(a0,a1,a2,a3,a4,a5,0,0,0,-1.5282)
    # F(x,y,z) = a0*x^2 + a1*y^2 + a2*z^2 + a3*x*y + a4*y*z + a5*x*z + a6*x + a7*y + a8*z + a9
 
    # sample the quadric function
    sample = vtk.vtkSampleFunction()
    sample.SetSampleDimensions(36,36,36)
    sample.SetImplicitFunction(quadric)
    sample.SetModelBounds(-2, 2, -2, 2, -2, 2)
 
    #create the 0 isosurface
    contours = vtk.vtkContourFilter()
    contours.SetInput(sample.GetOutput())
    contours.GenerateValues(1,1,1)

    transform = vtk.vtkTransform()
    transform.Translate(center)
 
    # map the contours to graphical primitives
    contourMapper = vtk.vtkPolyDataMapper()
    contourMapper.SetInput(contours.GetOutput())
    contourMapper.ScalarVisibilityOff()
 
    # create an actor for the contours
    contourActor = vtk.vtkActor()
    contourActor.SetMapper(contourMapper)
    contourActor.GetProperty().SetColor(color) # (R,G,B)
    contourActor.SetUserTransform(transform)
    return contourActor

def createStructureRenderer(surface,cell,param,param_use,rigid_bodies,atom_styles):
    global_parms, surface = \
            param_unfold(param, param_use, surface, False, False)
    surface = RB_update(rigid_bodies, surface, param, cell)
    actors = []
    M = calcM(cell)

    supercell = []
    for i in surface:
        if i[10] > 0:
            uxy = i[7] * (i[4] * i[5])**0.5
            uxz = i[8] * (i[4] * i[6])**0.5
            uyz = i[9] * (i[5] * i[6])**0.5            
            U = Num.array([[i[4],uxy,uxz],\
                           [uxy,i[5],uyz],\
                           [uxz,uyz,i[6]]],float)
            Uinv = Num.linalg.inv(U)
            supercell.append([i[0],Num.array([i[1],i[2],i[3]]),Uinv])
            supercell.append([i[0],Num.array([i[1]+1,i[2],i[3]]),Uinv])
            supercell.append([i[0],Num.array([i[1],i[2]+1,i[3]]),Uinv])
            supercell.append([i[0],Num.array([i[1]+1,i[2]+1,i[3]]),Uinv])


    #create an ellipsoid for every atom in supercell
    for i in range(len(supercell)):
        # create source
        point = supercell[i][1]
        [point] = aff2euk(M,[point])
        Uinv = supercell[i][2]
        if supercell[i][0] not in atom_styles.keys():
            color = [1,1,1]
        else:
            color = atom_styles[supercell[i][0]]
        actor = Ellipsoid(point, color, Uinv)
        actors.append(actor)
 
    # Create cell corner points
    points = []
    for i in [0.,1.]:
        for j in [0.,1.]:
            for k in [0.,1.]:
                points.append([i,j,k])
    points = aff2euk(M, points)

    #create cell edge lines
    line_id1 = [0,0,0,1,1,2,2,3,4,4,5,6]
    line_id2 = [1,2,4,3,5,3,6,7,5,6,7,7]
    for i in range(12):
        line = vtk.vtkLineSource()
        line.SetPoint1(points[line_id1[i]])
        line.SetPoint2(points[line_id2[i]])
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(line.GetOutput())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(1,1,1) # (R,G,B)
        actors.append(actor)

    # create a renderer
    ren = vtk.vtkRenderer()
    # assign actors to the renderer
    for actor in actors:
        ren.AddActor(actor)

    ren.SetBackground(0.0,0.0,0.0)
    ren.GetActiveCamera().ParallelProjectionOn()
    ren.ResetCamera()
    ren.GetActiveCamera().Azimuth(180)
    ren.GetActiveCamera().Elevation(270)
    ren.ResetCameraClippingRange()
    return ren



def saveStructureScreenshot(renWin, filename):
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(renWin)
    w2if.Update()
 
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(filename)
    writer.SetInput(w2if.GetOutput())
    writer.Write()    
    

