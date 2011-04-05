import numpy as Num
from ctrFitcalcs import *
import random

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

def contraction(Xmax, Xav, beta, used_params, parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters):
    print 'contraction'
    Xcon = beta*Xmax+(1-beta)*Xav
    parameter = insert(used_params, Xcon, parameter)
    dat, Ycon = calc_CTRs(parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters)
    return Xcon,Ycon

def reflection(Xmax, Xav, alpha, used_params, parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters):
    print 'reflection'
    Xref = (1+alpha)*Xav - alpha*Xmax
    Xref = check_limits(used_params, Xref, parameter)
    parameter = insert(used_params, Xref, parameter)
    dat, Yref = calc_CTRs(parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters)
    return Xref, Yref

def expansion(Xref, Xav, gamma, used_params, parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters):
    Xexp = (1+gamma)*Xref - gamma*Xav
    Xexp = check_limits(used_params, Xexp, parameter)
    parameter = insert(used_params, Xexp, parameter)
    dat, Yexp = calc_CTRs(parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters)
    return Xexp, Yexp

def compression(X, mini):
    print 'compression'
    Y = Num.ndarray((0,len(X[0])),float)
    for x in X:
        x = (x + X[mini])/2
        Y = Num.append(Y,[x],axis = 0)
    return Y

def calc_xdist(Xmin, Xmax):
    xdist = 0
    n = len(Xmin)
    for i in range(n):
        xdist =  xdist + (Xmax[i]-Xmin[i])**2
    xdist = xdist**0.5
    return xdist
    

#################################Simplex main routine###################################################################################################
def simplex(parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water,\
            simplex_params, use_BVC, BVclusters):
    alpha, beta, gamma, delta, ftol, xtol, maxiter = simplex_params

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
            dat, function_values[i] = calc_CTRs(parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters)
        else:
            for j in range(len(points[i])):
                key = used_params[j]
                points[i][j] = used_params_values[j] + random.uniform(((parameter[key][1]-used_params_values[j])*delta), ((parameter[key][2]-used_params_values[j])*delta))
            parameter = insert(used_params, points[i], parameter)    
            dat, function_values[i] = calc_CTRs(parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters)
            
    not_converged = True
    z = 0
    mini, maxi = min_max(function_values)
    while not_converged:     
        Xav = calc_average(points)
        Xref, Yref = reflection(points[maxi], Xav, alpha, used_params, parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters)
        if Yref < function_values[mini]:
            Xexp, Yexp = expansion(Xref, Xav, gamma, used_params, parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters)
            if Yexp < function_values[mini]:
                print 'expansion'
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
                    Xcon,Ycon = contraction(Xref, Xav, beta, used_params, parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters)
                else:
                    Xcon,Ycon = contraction(points[maxi], Xav, beta, used_params, parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters)

                if Ycon < function_values[maxi]:
                    points[maxi] = Xcon
                    function_values[maxi] = Ycon
                else:
                    points = compression(points, mini)
                    for i in range(len(points)):
                        dat, function_values[i] = calc_CTRs(parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters)
        mini, maxi = min_max(function_values)
        print 'iteration '+str(z)+', best RMS = '+str(round(function_values[mini],7))+', worst RMS = '+str(round(function_values[maxi],7))
        if function_values[mini] >= function_values[maxi]-ftol:
            not_converged = False
            print ' CONVERGENCE REACHED DUE TO FTOL \n\n'
        if calc_xdist(points[mini], points[maxi]) <= xtol:
            not_converged = False
            print ' CONVERGENCE REACHED DUE TO XTOL \n\n'
        if z >= maxiter:
            not_converged = False
            print ' NO CONVERGENCE, STOP DUE TO MAXITER \n\n'
        z = z+1
    param_best = points[mini]
    parameter = insert(used_params, param_best, parameter)
    data_best, RMS_best = calc_CTRs(parameter,param_usage, dat, cell, bulk_tmp, surface_tmp, NLayers, database, g_inv, Rod_weight, rigid_bodies, use_bulk_water, use_BVC, BVclusters)
    return data_best, parameter, RMS_best
################################################################################################################################
