"""
Functions used in pi-surf for resonant surface diffraction data refinement

Authors/modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""
###############################################################################

import numpy as Num
import random
import wx

from scipy.optimize import leastsq
from tdl.modules.sxrd.ctrfitcalcs import param_unfold
from tdl.modules.sxrd.simplex import insert, check_limits, calc_average, min_max, compression, calc_ftol

############################### methods used by simplex ############################################################################################
def contraction(Xmax, Xav, beta, used_params, allrasd, surface, parameter, parameter_usage, use_bulk_water, Refine_Data, use_lay_el):
    Xcon = beta*Xmax+(1-beta)*Xav
    parameter = insert(used_params, Xcon, parameter)
    allrasd = Rasd_difference(allrasd, surface, parameter, parameter_usage, use_bulk_water, Refine_Data, use_lay_el)
    Ycon = allrasd.RMS
    return Xcon,Ycon

def reflection(Xmax, Xav, alpha, used_params, allrasd, surface, parameter, parameter_usage, use_bulk_water, Refine_Data, use_lay_el):
    Xref = (1+alpha)*Xav - alpha*Xmax
    Xref = check_limits(used_params, Xref, parameter)
    parameter = insert(used_params, Xref, parameter)
    allrasd = Rasd_difference(allrasd, surface, parameter, parameter_usage, use_bulk_water, Refine_Data, use_lay_el)
    Yref = allrasd.RMS
    return Xref, Yref

def expansion(Xref, Xav, gamma, used_params, allrasd, surface, parameter, parameter_usage, use_bulk_water, Refine_Data, use_lay_el):
    Xexp = (1+gamma)*Xref - gamma*Xav
    Xexp = check_limits(used_params, Xexp, parameter)
    parameter = insert(used_params, Xexp, parameter)
    allrasd = Rasd_difference(allrasd, surface, parameter, parameter_usage, use_bulk_water, Refine_Data, use_lay_el)
    Yexp = allrasd.RMS
    return Xexp, Yexp
#################################Simplex main routine###################################################################################################
def res_simplex(StatusBar, parameter,param_usage, allrasd, surface, simplex_params, use_bulk_water, use_lay_el):

    alpha, beta, gamma, delta, ftol, maxiter = simplex_params
    Refine_Data = allrasd.Refine_Data
    StatusBar.SetStatusText('Preparing Simplex',0)
    used_params = []
    used_params_values = []
    for i in parameter.keys():
        if parameter[i][3]:
            used_params.append(i)
            used_params_values.append(parameter[i][0])        
    
    function_values = Num.ndarray((len(used_params)+1),float)
    points = Num.ndarray((len(used_params)+1,len(used_params)),float)
    
    for i in range(len(used_params)+1):
        if i == 0:
            points[i] = used_params_values
            parameter = insert(used_params, points[i], parameter)
            allrasd = Rasd_difference(allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el)
            function_values[i] = allrasd.RMS
        else:
            for j in range(len(points[i])):
                key = used_params[j]
                points[i][j] = used_params_values[j] + random.uniform(((parameter[key][1]-used_params_values[j])*delta), ((parameter[key][2]-used_params_values[j])*delta))
            parameter = insert(used_params, points[i], parameter)    
            allrasd = Rasd_difference(allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el)
            function_values[i] = allrasd.RMS
    not_converged = True
    z = 0
    mini, maxi = min_max(function_values)
    statusstring ='iteration '+str(z)+', best chi**2 = '+str(round(function_values[mini],4))
    StatusBar.SetStatusText(statusstring,0)
    while not_converged:
        StatusBar.SetStatusText(str(z),1)
        while wx.GetApp().Pending():
            wx.GetApp().Dispatch()
            wx.GetApp().Yield(True)
        old_minimum = function_values[mini]
        Xav = calc_average(points)
        Xref, Yref = reflection(points[maxi], Xav, alpha, used_params, allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el)
        if Yref < function_values[mini]:
            Xexp, Yexp = expansion(Xref, Xav, gamma, used_params, allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el)
            if Yexp < function_values[mini]:
                points[maxi] = Xexp
                function_values[maxi] = Yexp
            else:
                points[maxi] = Xref
                function_values[maxi] = Yref
        else:
            test = False
            for i in range(len(points)):
                if Yref < function_values[i]:
                    if i == maxi:
                        test = False
                    else:
                        test = True
            if test:
                points[maxi] = Xref
                function_values[maxi] = Yref
            else:
                if Yref < function_values[maxi]:
                    Xcon,Ycon = contraction(Xref, Xav, beta, used_params, allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el)
                else:
                    Xcon,Ycon =contraction(points[maxi], Xav, beta, used_params, allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el)

                if Ycon < function_values[maxi]:
                    points[maxi] = Xcon
                    function_values[maxi] = Ycon
                else:
                    points = compression(points, mini)
                    for i in range(len(points)):
                        parameter = insert(used_params, points[i], parameter)
                        allrasd = Rasd_difference(allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el)
                        function_values[i] = allrasd.RMS
        mini, maxi = min_max(function_values)
        act_ftol = calc_ftol(function_values)
        if function_values[mini]<old_minimum:
            statusstring ='iteration '+str(z)+', best chi**2 = '+str(round(function_values[mini],4))+' ftol = '+str(round(act_ftol,6))
            StatusBar.SetStatusText(statusstring,0)
        if act_ftol < ftol:
            not_converged = False
            print ' CONVERGENCE REACHED DUE TO FTOL \n\n'
        if z >= maxiter:
            not_converged = False
            print ' NO CONVERGENCE, STOP DUE TO MAXITER \n\n'
        z = z+1
    print 'best fit chi**2 = '+str(round(function_values[mini],7))+'\n'
    StatusBar.SetStatusText('End of Downhill Simplex, best chi**2: '+str(round(function_values[mini],4)),0)
    StatusBar.SetStatusText('',1)
    param_best = points[mini]
    parameter = insert(used_params, param_best, parameter)
    allrasd = Rasd_difference(allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el)
    return allrasd, parameter
################################################################################################################################
##################  Refinement of Atom coordinates, occupancies and DW- Factors  ###################################################################

def Wiggle2(vec, Rasd):
    a, b = vec
    delta_F = Rasd.F - ((a+b*(Rasd.E-Rasd.E0))*((Rasd.re_FNR + Rasd.re_FR)**2 + (Rasd.im_FNR + Rasd.im_FR)**2))
    return delta_F

def Jacobi2(vec, Rasd):
    J = Num.ndarray((2,len(Rasd.E)),float)
    J[0][0:] = -(Rasd.re_FNR + Rasd.re_FR)**2 - (Rasd.im_FNR + Rasd.im_FR)**2
    J[1][0:] = -(Rasd.E-Rasd.E0) * ((Rasd.re_FNR + Rasd.re_FR)**2 + (Rasd.im_FNR + Rasd.im_FR)**2)
    return J

### main function that calculates the difference between measured and calculated F**2s for a set of resonant atoms ###
def Rasd_difference(allrasd, surface_tmp, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el, fix_bckg = False):
    allrasd.RMS = 0
    allrasd.ndata = 0
    
    global_parms, surface = param_unfold(parameter,param_usage, surface_tmp, use_bulk_water, use_lay_el)
    natoms = len(surface)
    occ_el, K,sig_el,sig_el_bar,d_el,d0_el,sig_water, sig_water_bar, d_water, zwater, Scale, specScale, beta= global_parms
    
    for Rasd in allrasd.list:
        if Rasd.use_in_Refine:
            Q = Num.zeros((3),float)
            Q[0] = allrasd.g_inv[0][0]**0.5 * Rasd.Q[0]
            Q[1] = allrasd.g_inv[1][1]**0.5 * Rasd.Q[1]
            Q[2] = allrasd.g_inv[2][2]**0.5 * Rasd.Q[2]
        
            j = 0
            re_Fq = 0
            im_Fq = 0
            for n in range(natoms):
                R = Num.array([surface[n][1],surface[n][2],surface[n][3]],float)

                U = Num.ndarray((3,3),float)
                U[0][0] = surface[n][4]
                U[0][1] = U[1][0] = surface[n][7]*(surface[n][4]*surface[n][5])**0.5
                U[0][2] = U[2][0] = surface[n][8]*(surface[n][4]*surface[n][6])**0.5
                U[1][1] = surface[n][5]
                U[1][2] = U[2][1] = surface[n][9]*(surface[n][5]*surface[n][6])**0.5
                U[2][2] = surface[n][6]
                
                theta = surface[n][10]
                re_Fq = re_Fq + (theta * Num.cos(2*Num.pi*(Rasd.Q[0]*R[0]+Rasd.Q[1]*R[1]+Rasd.Q[2]*R[2])) * Num.exp(-2* Num.pi**2 * Num.dot(Num.dot(Q,U),Q)))
                im_Fq = im_Fq + (theta * Num.sin(2*Num.pi*(Rasd.Q[0]*R[0]+Rasd.Q[1]*R[1]+Rasd.Q[2]*R[2])) * Num.exp(-2* Num.pi**2 * Num.dot(Num.dot(Q,U),Q)))
            if use_lay_el and Q[0] == 0 and Q[1] == 0:
                re_lay, im_lay = calc_F_lay_el(Rasd.Q, occ_el, K, sig_el, sig_el_bar, d_el, d0_el, allrasd.g_inv)
                re_Fq = re_Fq + re_lay
                im_Fq = im_Fq + im_lay                                    
            Rasd.re_FR = Rasd.f1 * re_Fq - Rasd.f2 * im_Fq
            Rasd.im_FR = Rasd.f1 * im_Fq + Rasd.f2 * re_Fq

            Rasd.re_Fq = re_Fq
            Rasd.im_Fq = im_Fq

            if Refine_Data:
                if not fix_bckg:
                    vec = [Rasd.a, Rasd.b]
                    result = leastsq(Wiggle2, vec, args=(Rasd), Dfun = Jacobi2, col_deriv = 1)
                    vec = result[0]
                    Rasd.a, Rasd.b = vec

                Rasd.F_calc = (Rasd.a+Rasd.b*(Rasd.E - Rasd.E0))*((Rasd.re_FNR + Rasd.re_FR)**2 + (Rasd.im_FNR + Rasd.im_FR)**2) 
                Rasd.delta_F = ((Rasd.F - Rasd.F_calc)/Rasd.Ferr)**2
                Rasd.norm = (Rasd.a+ Rasd.b* (Rasd.E - Rasd.E0))*(Rasd.re_FNR**2 + Rasd.im_FNR**2)
                Rasd.RMS = Num.sum(Rasd.delta_F)/Rasd.ndata
        
                allrasd.RMS = allrasd.RMS + Num.sum(Rasd.delta_F)
                allrasd.ndata = allrasd.ndata + Rasd.ndata
            else:
                PR = Num.arctan(Rasd.im_Fq/Rasd.re_Fq)/(2*Num.pi)
                AR = Rasd.re_Fq /Num.cos(2*Num.pi*PR)
                if AR < 0 and PR < 0.5 and PR > 0.:
                    PR = PR + 0.5
                    AR = -AR
                elif AR < 0 and PR > 0.5:
                    PR = PR - 0.5
                    AR = -AR
                elif PR < 0 and AR > 0:
                    PR = PR + 1.
                elif AR < 0 and PR < 0 and PR > -0.5:
                    PR = PR + 0.5
                    AR = -AR
                if PR > 1: PR = PR -1
                Rasd.AR_refine = AR
                Rasd.PR_refine = PR

                Rasd.RMS = (Rasd.re_Fq / Num.cos(2*Num.pi*Rasd.PR) - Rasd.AR)**2 + (Rasd.im_Fq / Num.sin(2*Num.pi*Rasd.PR) - Rasd.AR)**2
                allrasd.RMS = allrasd.RMS + Rasd.RMS
                allrasd.ndata = allrasd.ndata + 1
        
    allrasd.RMS = allrasd.RMS/allrasd.ndata

    return allrasd

def calc_F_lay_el(hkl, occ, K, sig, sig_bar, d, d0, g_inv):
    q = hkl[2]* g_inv[2][2]**0.5
    f = Num.exp(-2 * Num.pi**2 * q**2 * sig)*occ
    x = Num.pi * q * d
    al = 2 * Num.pi**2 * q**2 * sig_bar + K * d
    a = Num.exp(al)*Num.cos(2*x)-1
    b = Num.exp(al)*Num.sin(-2*x)
    c = 4 * Num.cos(x)**2 * Num.sinh(al/2)**2 - 4 * Num.sin(x)**2 * Num.cosh(al/2)**2
    d = -2 * Num.sin(2*x) * Num.sinh(al)
    wert = 2*Num.pi*hkl[2]*d0
    rez = Num.cos(wert)
    imz = Num.sin(wert)
    wert = c**2 + d**2
    relayer = (a*c + b*d)/(wert)
    imlayer = (b*c - a*d)/(wert)
    re = f* (relayer * rez - imlayer * imz)
    im = f* (relayer * imz + imlayer * rez)
    return re, im
####################################################################################################################################################
def res_param_statistics(fpc, parameter,param_usage, allrasd, surface, use_bulk_water, use_lay_el):
    n = 0
    used_params = []
    for i in range(len(allrasd.list)):
        n = n + allrasd.list[i].ndata
    w = Num.zeros((n,n))
    n = 0
    for i in range(len(allrasd.list)):
        for j in range(allrasd.list[i].ndata):
            w[n][n] = (1/allrasd.list[i].Ferr[j])**2
            n = n+1
                   
    for i in parameter.keys():
        if parameter[i][3]:
            used_params.append(i)
    
    b = len(used_params)
    X = Num.zeros((b,n))
    Refine_Data = allrasd.Refine_Data

    def extract_values(data):
        y = Num.array([])
        for i in range(len(data.list)):
            for j in range(len(data.list[i].F)):
                y = Num.append(y, data.list[i].F_calc[j])
        return y
    
    for i in range(b):
        h = parameter[used_params[i]][0] * fpc
        if h == 0:
            h = fpc
        elif h < 0:
            h = -h
        if parameter[used_params[i]][0] -0.5*h >= parameter[used_params[i]][1] and parameter[used_params[i]][0] +0.5*h <= parameter[used_params[i]][2]:
            parameter[used_params[i]][0] = parameter[used_params[i]][0] -0.5*h
            allrasd = Rasd_difference(allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el)
            y1 = extract_values(allrasd)
            parameter[used_params[i]][0] = parameter[used_params[i]][0] +h
            allrasd = Rasd_difference(allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el, True)
            y2 = extract_values(allrasd)
            y1 = (y2 - y1)/h
            parameter[used_params[i]][0] = parameter[used_params[i]][0] -0.5*h
            
        elif parameter[used_params[i]][0] -0.5*h < parameter[used_params[i]][1]:
            allrasd = Rasd_difference(allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el)
            y1 = extract_values(allrasd)
            parameter[used_params[i]][0] = parameter[used_params[i]][0] +h
            allrasd = Rasd_difference(allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el, True)
            y2 = extract_values(allrasd)
            y1 = (y2 - y1)/h
            parameter[used_params[i]][0] = parameter[used_params[i]][0] -h
            
        elif parameter[used_params[i]][0] +0.5*h > parameter[used_params[i]][2]:
            parameter[used_params[i]][0] = parameter[used_params[i]][0] -h
            allrasd = Rasd_difference(allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el)
            y1 = extract_values(allrasd)
            parameter[used_params[i]][0] = parameter[used_params[i]][0] +h
            allrasd = Rasd_difference(allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el, True)
            y2 = extract_values(allrasd)
            y1 = (y2 - y1)/h
            
        for j in range(n):
            X[i][j] = y1[j] * parameter[used_params[i]][0] * Num.sqrt(w[j][j])

    allrasd = Rasd_difference(allrasd, surface, parameter, param_usage, use_bulk_water, Refine_Data, use_lay_el)
            
    V = Num.dot(X, Num.dot( w, Num.transpose(X)))
    C = Num.zeros((b,b))
    try:
        V = allrasd.RMS * Num.linalg.inv(V)
        for i in range(b):
            for j in range(b):
                if i == j:
                    C[i][j] = Num.sqrt(V[i][j])
                    parameter[used_params[i]][4] = C[i][j]
                else:
                    C[i][j] = V[i][j]/Num.sqrt(V[i][i]*V[j][j])
        keys = parameter.keys()
        for key in keys:
            if parameter[key][5] in keys:
                parameter[key][4] = parameter[parameter[key][5]][4]
    except:
        print "There's a hole in the matrix !!! Omit unused or insignificant parameters in the calculation.\n"
        
    return parameter, C, used_params
    

    
