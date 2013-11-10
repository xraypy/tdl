"""
Nelder Mead Simplex optimization routine used in pi-surf

Authors/modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""
###############################################################################

import numpy as Num
import random
import wx

from tdl.modules.sxrd.ctrfitcalcs import *
############################### methods used by simplex ############################################################################################
def insert(used_params, point, parameter):
    for i in range(len(used_params)):
        key = used_params[i]
        parameter[key][0] = point[i]
    return parameter

def check_limits(used_params, point, parameter):
    for i in range(len(used_params)):
        key = used_params[i]
        if point[i] < parameter[key][1]: 
            point[i] = parameter[key][1]
        elif point[i] > parameter[key][2]: 
            point[i] = parameter[key][2]
    return point
    
def calc_average(points):
    av_point = Num.zeros((len(points[0])),float)
    for i in points:
        av_point = av_point + i
    av_point = av_point/(len(points))
    return av_point

def min_max(function_values):
    Ymin = function_values.min()
    mini = int(Num.where(function_values == Ymin)[0][0])
    Ymax = function_values.max()
    maxi = int(Num.where(function_values == Ymax)[0][0])
    return mini, maxi

def contraction(Xmax, Xav, beta, used_params, parameter,param_usage, dat, cell, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el):
    #print 'contraction'
    Xcon = beta*Xmax+(1-beta)*Xav
    parameter = insert(used_params, Xcon, parameter)
    dat, Ycon = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
    return Xcon,Ycon

def reflection(Xmax, Xav, alpha, used_params, parameter,param_usage, dat, cell, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el):
    #print 'reflection'
    Xref = (1+alpha)*Xav - alpha*Xmax
    Xref = check_limits(used_params, Xref, parameter)
    parameter = insert(used_params, Xref, parameter)
    dat, Yref = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
    return Xref, Yref

def expansion(Xref, Xav, gamma, used_params, parameter,param_usage, dat, cell, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el):
    Xexp = (1+gamma)*Xref - gamma*Xav
    Xexp = check_limits(used_params, Xexp, parameter)
    parameter = insert(used_params, Xexp, parameter)
    dat, Yexp = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
    return Xexp, Yexp

def compression(X, mini):
    #print 'compression'
    Y = Num.ndarray((0,len(X[0])),float)
    for x in X:
        x = (x + X[mini])/2
        Y = Num.append(Y,[x],axis = 0)
    return Y

def parameter_plot(fig,used_params,parameter,points,mini):
    if fig == None:
        fig = figure(3,figsize=[15,5])
    fig.clear()
    fig.suptitle('Parameter Plot', fontsize = 20)
    plot = fig.add_subplot(111)
    plot.set_xticks(range(len(used_params)))
    plot.set_xticklabels(used_params,rotation = 90)
    low = []
    spread = []
    for i in used_params:
        low.append(parameter[i][1])
        spread.append(parameter[i][2]-parameter[i][1])
    for i in range(len(points)):
        if i == mini:
            pass
        else:
            plot.plot(range(len(used_params)),(points[i]-low)/spread,'bo')
        plot.plot(range(len(used_params)),(points[mini]-low)/spread,'ro')
    plot.set_xlim(-1,len(used_params)+1)
    plot.vlines(range(len(used_params)),ymin = 0,ymax = 1, color = 'k', linestyles = 'dashed')
    plot.set_ylim(0,1)
    plot.figure.canvas.draw()
    return fig

def calc_ftol(function_values):
    av = Num.sum(function_values)/len(function_values)
    sigma = 0
    for i in function_values:
        sigma = sigma + (i-av)**2
    sigma = Num.sqrt(sigma/len(function_values))
    return sigma
#################################Simplex main routine###################################################################################################
def simplex(StatusBar,parameter,param_usage, dat, cell, surface_tmp, NLayers, database, rigid_bodies, panel):
    Rod_weight = panel.Rod_weight
    g_inv = calc_g_inv(cell)
    use_bulk_water = panel.UBW_flag
    use_BVC = panel.use_BVC
    BVclusters = panel.BVclusters
    RMS_flag = panel.RMS_flag
    use_lay_el = panel.use_lay_el
    el = panel.el
    alpha, beta, gamma, delta, ftol, maxiter, random_pars = panel.simplex_params
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
        while wx.GetApp().Pending():
            wx.GetApp().Dispatch()
            wx.GetApp().Yield(True)
        if i == 0 and not random_pars:
            points[i] = used_params_values
            parameter = insert(used_params, points[i], parameter)
            dat, function_values[i] = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
        else:
            for j in range(len(points[i])):
                key = used_params[j]
                points[i][j] = used_params_values[j] + random.uniform(((parameter[key][1]-used_params_values[j])*delta), ((parameter[key][2]-used_params_values[j])*delta))
            parameter = insert(used_params, points[i], parameter)    
            dat, function_values[i] = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
            
    not_converged = True
    z = 0
    mini, maxi = min_max(function_values)
    statusstring ='iteration '+str(z)+', best chi**2 = '+str(round(function_values[mini],4))
    StatusBar.SetStatusText(statusstring,0)
    panel.Figure3 = parameter_plot(panel.Figure3,used_params,parameter,points,mini)
    old_mini = function_values[mini]
    while not_converged:
        StatusBar.SetStatusText(str(z),1)
        Xav = calc_average(points)
        while wx.GetApp().Pending():
            wx.GetApp().Dispatch()
            wx.GetApp().Yield(True)
        if panel.StopFit:
            StatusBar.SetStatusText('Fit stopped after '+str(z)+' iterations',0)
            not_converged = False
            print 'Fit aborted by user after '+str(z)+' iterations'
        Xref, Yref = reflection(points[maxi], Xav, alpha, used_params, parameter,param_usage, dat, cell, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
        if Yref < function_values[mini]:
            Xexp, Yexp = expansion(Xref, Xav, gamma, used_params, parameter,param_usage, dat, cell, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
            if Yexp < function_values[mini]:
                #print 'expansion'
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
                    Xcon,Ycon = contraction(Xref, Xav, beta, used_params, parameter,param_usage, dat, cell, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
                else:
                    Xcon,Ycon = contraction(points[maxi], Xav, beta, used_params, parameter,param_usage, dat, cell, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)

                if Ycon < function_values[maxi]:
                    points[maxi] = Xcon
                    function_values[maxi] = Ycon
                else:
                    points = compression(points, mini)
                    for i in range(len(points)):
                        parameter = insert(used_params, points[i], parameter) 
                        dat, function_values[i] = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
        mini, maxi = min_max(function_values)
        act_ftol = calc_ftol(function_values)
        if function_values[mini]<old_mini:
            statusstring ='iteration '+str(z)+', best chi**2 = '+str(round(function_values[mini],4))+' ftol = '+str(round(act_ftol,6))
            StatusBar.SetStatusText(statusstring,0)
            panel.Figure1 = plot_rods(panel.Figure1, StatusBar.nb.data, \
                                                             StatusBar.nb.MainControlPage.plotdims, StatusBar.nb.MainControlPage.doplotbulk,\
                                                             StatusBar.nb.MainControlPage.doplotsurf, StatusBar.nb.MainControlPage.doplotrough,\
                                                             StatusBar.nb.MainControlPage.doplotwater, function_values[mini])
            panel.Figure1.canvas.draw()
            panel.Figure3 = parameter_plot(panel.Figure3,used_params,parameter,points,mini)
            old_mini = function_values[mini]
            
        if act_ftol < ftol:
            not_converged = False
            print '\n CONVERGENCE REACHED DUE TO FTOL \n'
        if z >= maxiter:
            not_converged = False
            print '\n NO CONVERGENCE, STOP DUE TO MAXITER \n'
        z = z+1
    if not panel.StopFit: print ' Downhill Simplex stopped after '+str(z-1)+' iterations'
    print 'best fit chi**2 = '+str(round(function_values[mini],7))+'\n'
    StatusBar.SetStatusText('End of Downhill Simplex, best chi**2: '+str(round(function_values[mini],4)),0)
    StatusBar.SetStatusText('',1)
    param_best = points[mini]
    parameter = insert(used_params, param_best, parameter)
    data_best, RMS_best = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
    return data_best, parameter, RMS_best

#################################################################################################################################
def extract_values(data):
    y = Num.array([])
    for i in range(len(data)):
        for j in range(len(data[i].L)):
            y = Num.append(y, data[i].Fcalc[j])
    return y

def statistics(fpc, parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el):
    n = 0
    used_params = []
    for i in range(len(dat)):
        n = n + len(dat[i].L)
    w = Num.zeros((n,n))
    n = 0
    for i in range(len(dat)):
        for j in range(len(dat[i].L)):
            w[n][n] = (1/dat[i].Ferr[j])**2
            n = n+1
                   
    for i in parameter.keys():
        if parameter[i][3]:
            used_params.append(i)
    
    b = len(used_params)
    X = Num.zeros((b,n))
    
    for i in range(b):
        h = parameter[used_params[i]][0] * fpc
        if h == 0:
            h = fpc
        elif h < 0.:
            h = -h
        if parameter[used_params[i]][0] -0.5*h >= parameter[used_params[i]][1] and parameter[used_params[i]][0] +0.5*h <= parameter[used_params[i]][2]:
            parameter[used_params[i]][0] = parameter[used_params[i]][0] -0.5*h
            data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
            y1 = extract_values(data)
            parameter[used_params[i]][0] = parameter[used_params[i]][0] +h
            data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
            y2 = extract_values(data)
            y1 = (y2 - y1)/h
            parameter[used_params[i]][0] = parameter[used_params[i]][0] -0.5*h
            
        elif parameter[used_params[i]][0] -0.5*h < parameter[used_params[i]][1]:
            data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
            y1 = extract_values(data)
            parameter[used_params[i]][0] = parameter[used_params[i]][0] +h
            data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
            y2 = extract_values(data)
            y1 = (y2 - y1)/h
            parameter[used_params[i]][0] = parameter[used_params[i]][0] -h
            
        elif parameter[used_params[i]][0] +0.5*h > parameter[used_params[i]][2]:
            parameter[used_params[i]][0] = parameter[used_params[i]][0] -h
            data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
            y1 = extract_values(data)
            parameter[used_params[i]][0] = parameter[used_params[i]][0] +h
            data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
            y2 = extract_values(data)
            y1 = (y2 - y1)/h
            
        for j in range(n):
            X[i][j] = y1[j] * parameter[used_params[i]][0] * Num.sqrt(w[j][j])

    data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
            
    V = Num.dot(X, Num.dot( w, Num.transpose(X)))
    C = Num.zeros((b,b))
    try:
        V = R* Num.linalg.inv(V)
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
        
    return parameter, X, C, used_params, R

def single_param_sensitivities(param_label, fpc, parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el):
    n = 0
    for i in range(len(dat)):
        n = n + len(dat[i].L)
    w = Num.zeros((n,n))
    n = 0
    for i in range(len(dat)):
        for j in range(len(dat[i].L)):
            w[n][n] = (1/dat[i].Ferr[j])**2
            n = n+1

    X = Num.zeros((n))
    
    h = parameter[param_label][0] * fpc
    if h == 0.:
        h = fpc
    elif h < 0.:
        h = -h
    if parameter[param_label][0] -0.5*h >= parameter[param_label][1] and parameter[param_label][0] +0.5*h <= parameter[param_label][2]:
        parameter[param_label][0] = parameter[param_label][0] -0.5*h
        data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
        y1 = extract_values(data)
        parameter[param_label][0] = parameter[param_label][0] +h
        data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
        y2 = extract_values(data)
        y1 = (y2 - y1)/h
        parameter[param_label][0] = parameter[param_label][0] -0.5*h
            
    elif parameter[param_label][0] -0.5*h < parameter[param_label][1]:
        data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
        y1 = extract_values(data)
        parameter[param_label][0] = parameter[param_label][0] +h
        data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
        y2 = extract_values(data)
        y1 = (y2 - y1)/h
        parameter[param_label][0] = parameter[param_label][0] -h
            
    elif parameter[param_label][0] +0.5*h > parameter[param_label][2]:
        parameter[param_label][0] = parameter[param_label][0] -h
        data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
        y1 = extract_values(data)
        parameter[param_label][0] = parameter[param_label][0] +h
        data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
        y2 = extract_values(data)
        y1 = (y2 - y1)/h
            
    for j in range(n):
        X[j] = y1[j] * parameter[param_label][0] * Num.sqrt(w[j][j])

    data, R = calc_CTRs(parameter,param_usage, dat, cell,surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters, RMS_flag, use_lay_el, el)
    V = Num.dot(X,Num.dot(w,Num.transpose(X)))
    if V > 0:
        dp = Num.sqrt(R/V)
    else:
        dp = 0.

    return X, dp
                              
                              
                              
            
    
    
