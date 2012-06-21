"""
Genetic algorithm optimization routine used in pi-surf

Authors/modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

"""
##############################################################################

import random
import wx
import numpy as Num

from tdl.modules.sxrd.ctrfitcalcs import *
from tdl.modules.sxrd.simplex import insert, parameter_plot, calc_xdist 

###############################################################################
def sort_population(fitness, individuals):
    ranking = Num.argsort(fitness)
    fitness_sorted = []
    individuals_sorted = []
    for i in range(len(fitness)):
        fitness_sorted.append(fitness[ranking[i]])
        individuals_sorted.append(individuals[ranking[i]])
    return Num.array(fitness_sorted), Num.array(individuals_sorted)

################################################################################
def genetic(panel, dat, cell, surface, NLayers, database,\
            rigid_bodies, parameter, param_usage):

    km, kc, popmult, maxgen, jump, no_improve, xtol, ftol, \
        random_pars = panel.genetic_params
    g_inv = calc_g_inv(cell)
    Rod_weight = panel.Rod_weight
    use_bulk_water = panel.UBW_flag
    use_BVC = panel.use_BVC
    use_lay_el = panel.use_lay_el
    el = panel.el
    BVclusters = panel.BVclusters
    RMS_flag = panel.RMS_flag
    fig1 = panel.Figure1
    fig3 = panel.Figure3
    plot_dims = panel.plotdims
    plot_bulk = panel.doplotbulk
    plot_surf = panel.doplotsurf
    plot_rough = panel.doplotrough
    plot_water = panel.doplotwater
    statusbar = panel.nb.frame.statusbar

    statusbar.SetStatusText('Preparing Population',0)
    used_params = []
    used_params_values = []
    for i in parameter.keys():
        if parameter[i][3]:
            used_params.append(i)
            used_params_values.append(parameter[i][0])

    n = int(len(used_params)*popmult)
    m = len(used_params)
    if n < 15: n = 15
    fitness = Num.ndarray((n),float)
    individuals = Num.ndarray((n,m), float)
    
    for i in range(n):
        while wx.GetApp().Pending():
            wx.GetApp().Dispatch()
            wx.GetApp().Yield(True)
        if i == 0 and not random_pars:
            individuals[i] = used_params_values
            parameter = insert(used_params, individuals[i], parameter)
            dat, fitness[i] = calc_CTRs(parameter,param_usage,\
                                                dat, cell,surface,\
                                                NLayers, database, g_inv,\
                                                Rod_weight, rigid_bodies,\
                                                use_bulk_water, use_BVC,\
                                                BVclusters, RMS_flag,\
                                                use_lay_el, el)
        else:
            for j in range(m):
                key = used_params[j]
                individuals[i][j] = random.uniform(parameter[key][1],parameter[key][2])
            parameter = insert(used_params, individuals[i], parameter)    
            dat, fitness[i] = calc_CTRs(parameter,param_usage,\
                                                dat, cell,surface,\
                                                NLayers, database, g_inv,\
                                                Rod_weight, rigid_bodies,\
                                                use_bulk_water, use_BVC,\
                                                BVclusters, RMS_flag,\
                                                use_lay_el, el)

    z = 0
    fitness, individuals = sort_population(fitness,individuals)
    fig3 = parameter_plot(fig3,used_params,parameter,individuals,0)
    best = fitness[0]
    best_count = 0
    not_converged = True
    while z < maxgen and not_converged:
        while wx.GetApp().Pending():
            wx.GetApp().Dispatch()
            wx.GetApp().Yield(True)
        statusstring ='generation '+str(z)+', best R = '\
                       +str(round(fitness[0],7))
        statusbar.SetStatusText(statusstring,0)
        # mutation
        for i in range(n):
            while wx.GetApp().Pending():
                wx.GetApp().Dispatch()
                wx.GetApp().Yield(True)
            if random.uniform(0,1) < km:
                individual = individuals[i][:]
                j = random.randrange(m)
                key = used_params[j]
                individual[j] = individuals[i][j]+ random.uniform(((parameter[key][1]-individuals[i][j])*jump),\
                                                                  ((parameter[key][2]-individuals[i][j])*jump))
                    
                parameter = insert(used_params, individual, parameter)
                dat, fit = calc_CTRs(parameter,param_usage,\
                                                dat, cell,surface,\
                                                NLayers, database, g_inv,\
                                                Rod_weight, rigid_bodies,\
                                                use_bulk_water, use_BVC,\
                                                BVclusters, RMS_flag,\
                                                use_lay_el, el)
                individuals = Num.append(individuals,[individual],  axis =0)
                fitness = Num.append(fitness, fit)
                statusbar.SetStatusText('mutation '+str(i),1)
        #sort before crossover
        fitness, individuals = sort_population(fitness, individuals)
        # crossover
        for i in range(int(n*kc/2)):
            while wx.GetApp().Pending():
                wx.GetApp().Dispatch()
                wx.GetApp().Yield(True)
            statusbar.SetStatusText('crossover '+str(i),1)
            choice1 = random.randrange(random.randrange(5,len(fitness)))
            choice2 = random.randrange(random.randrange(5,len(fitness)))
            while choice2 == choice1:
                choice2 = random.randrange(random.randrange(5,len(fitness)))
            crosspoint1 = random.randrange(1,m-2)
            crosspoint2 = random.randrange(crosspoint1,m-1)
            selection1 = individuals[choice1]
            selection2 = individuals[choice2]
            x1 = selection1[:crosspoint1]
            y1 = selection2[:crosspoint1]
            x2 = selection2[crosspoint1:crosspoint2]
            y2 = selection1[crosspoint1:crosspoint2]
            x3 = selection1[crosspoint2:]
            y3 = selection2[crosspoint2:]
            child1 = Num.append(x1,x2)
            child1 = Num.append(child1,x3)
            child2 = Num.append(y1,y2)
            child2 = Num.append(child2,y3)
            parameter = insert(used_params, child1, parameter)
            dat, fit_child1 = calc_CTRs(parameter,param_usage,\
                                                dat, cell,surface,\
                                                NLayers, database, g_inv,\
                                                Rod_weight, rigid_bodies,\
                                                use_bulk_water, use_BVC,\
                                                BVclusters, RMS_flag,\
                                                use_lay_el, el)
            parameter = insert(used_params, child2, parameter)
            dat, fit_child2 = calc_CTRs(parameter,param_usage,\
                                                dat, cell,surface,\
                                                NLayers, database, g_inv,\
                                                Rod_weight, rigid_bodies,\
                                                use_bulk_water, use_BVC,\
                                                BVclusters, RMS_flag,\
                                                use_lay_el, el)
            fitness = Num.append(fitness, fit_child1)
            fitness = Num.append(fitness, fit_child2)
            individuals = Num.append(individuals, Num.array([child1]), axis = 0)
            individuals = Num.append(individuals, Num.array([child2]), axis = 0)
        # sort population 
        fitness, individuals = sort_population(fitness, individuals)
        parameter = insert(used_params, individuals[0], parameter)
        dat, fittest = calc_CTRs(parameter,param_usage,\
                                                dat, cell,surface,\
                                                NLayers, database, g_inv,\
                                                Rod_weight, rigid_bodies,\
                                                use_bulk_water, use_BVC,\
                                                BVclusters, RMS_flag,\
                                                use_lay_el, el)
        #make sure the sorting is correct
        while fitness[0] != fittest:
            fitness, individuals = sort_population(fitness, individuals)
            parameter = insert(used_params, individuals[0], parameter)
            dat, fittest = calc_CTRs(parameter,param_usage,\
                                                dat, cell,surface,\
                                                NLayers, database, g_inv,\
                                                Rod_weight, rigid_bodies,\
                                                use_bulk_water, use_BVC,\
                                                BVclusters, RMS_flag,\
                                                use_lay_el, el)
        #and keep only the n fittest individuals
        fitness = fitness[:n]
        individuals = individuals[:n]
        #report intermediate results
        fig3 = parameter_plot(fig3,used_params,parameter,individuals,0)
        fig1 = plot_rods(fig1, dat, plot_dims, plot_bulk, plot_surf,\
                         plot_rough, plot_water, fittest)
        fig1.canvas.draw()
        statusstring ='generation '+str(z)+', best R = '\
                       +str(round(fittest,7))
        statusbar.SetStatusText(statusstring,0)
        if not fittest < best:
            best_count = best_count+1
        else:
            best = fittest
            best_count = 0
        if panel.StopFit:
            statusbar.SetStatusText('Fit stopped after '+str(z)+\
                                    ' generations',0)
            print 'Fit aborted by user after '+str(z)+' generations'
            not_converged = False
        if calc_xdist(individuals[0],individuals[n-1],used_params,parameter) < xtol:
            not_converged = False
            print '\n CONVERGENCE REACHED DUE TO XTOL\n'
        if fitness[n-1]-fitness[0] < ftol:
            not_converged = False
            print '\n CONVERGENCE REACHED DUE TO FTOL\n'
        if best_count > no_improve:
            not_converged = False
            print '\n EVOLUTION STAGNATED: NO IMPROVEMENT FOR '+str(no_improve)+' GENERATIONS\n'
        z = z+1
        if z == maxgen:
            print '\n STOP DUE TO MAXITER\n'
    print 'genetic algorithm finished after '+str(z)+' generations,\nfittest individual R = '+str(round(fitness[0],7))
    return dat, parameter, fittest
    
            
            
    
            
            
            
            
        
        


