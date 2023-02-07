#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: D.Brueckner
"""

import numpy as np


def corr_msd(X,N_t_max,delta_t,frac_part=1): #X[j,t]
    sum_tot = np.zeros(N_t_max) #only calculate up to N_t_max to save time
    
    Ntrajectories = X.shape[0]
    if(len(X.shape)==2):
        N_dim = 1
    elif(len(X.shape)==3):
        N_dim = X.shape[2]
    
    MSD = np.zeros(N_t_max)
    MSD_err = np.zeros(N_t_max)
    MSD_all = np.zeros((Ntrajectories,N_t_max))
    N_tot_delta = np.zeros(N_t_max)
    for j in range(0,Ntrajectories):
        if(N_dim>1):
            tracklength = (~np.isnan(X[j,:,0])).cumsum(0).argmax(0)+1
        elif(N_dim==1):
            tracklength = (~np.isnan(X[j,:])).cumsum(0).argmax(0)+1
        
        if(tracklength<N_t_max):
            delta_max = tracklength
        else:
            delta_max = N_t_max
        
        MSD_j = np.zeros(N_t_max)
        MSD_sum_j = np.zeros(N_t_max)
        for delta in range(0,delta_max):
            N_tot_j = 0
            N = tracklength - delta
            sum_MSDj = 0
            for t in range(0,N):
                
                if(N_dim>1):
                    if(~np.isnan(X[j,t+delta,0]) and ~np.isnan(X[j,t,0])):
                        N_tot_delta[delta] += 1
                        N_tot_j += 1
                        for d in range(0,N_dim):
                            sum_MSDj += (X[j,t+delta,d]-X[j,t,d])**2
                elif(N_dim==1):
                    if(~np.isnan(X[j,t+delta]) and ~np.isnan(X[j,t])):
                        N_tot_delta[delta] += 1
                        N_tot_j += 1
                        sum_MSDj += (X[j,t+delta]-X[j,t])**2
            if(N_tot_j>0):    
                MSD_j[delta] = sum_MSDj/N_tot_j
            else:
                MSD_j[delta] = np.nan
            MSD_sum_j[delta] = sum_MSDj
        
        sum_tot += MSD_sum_j
        MSD_all[j,:] = MSD_j
    
    for delta in range(0,N_t_max):
        if(N_tot_delta[delta]>10):
            MSD[delta] = sum_tot[delta]/N_tot_delta[delta] 
            
            N_trajs_delta = sum(~np.isnan(MSD_all[:,delta]))
            MSD_err[delta] = np.sqrt(np.nanvar(MSD_all[:,delta])/N_trajs_delta)
        else:
            MSD[delta] = np.nan
            MSD_err[delta] = np.nan
    time = np.linspace(0,N_t_max-1,N_t_max)*delta_t
    
    return MSD,MSD_all,MSD_err,time


def fit_2loc_msd_rouse(data,mode_error='z'):
    
    array_data = [] #<list of (T, d) numpy arrays>
    for j in range(0,data.N_part):
        array_data.append(data.Xij[j,data.t_start[j]:data.t_end[j],:])
    
    import tracklib as tl
    data_tracklib = tl.TaggedSet()
    
    #transform data into tracklib format
    for traj in array_data:
        data_tracklib.add(tl.Trajectory.fromArray(traj))
     
    # Set up the Rouse fit
    fit = tl.analysis.msdfit.lib.TwoLocusRouseFit(data_tracklib)
    
    # determine which localisation errors should be included in the fit
    if mode_error == 'z':
        fit.fix_values += [(0, -np.inf), (3, -np.inf)]
    elif mode_error == 'none':
        fit.fix_values += [(0, -np.inf), (3, -np.inf), (6, -np.inf)]
    
    # Run
    res = fit.run()
    ## res is a dict with fields ‘logL’ and ‘params’; note that ‘params’ is in log-space (see docstring of TwoLocusRouseFit)
    params = res['params']
    
    sigma = np.sqrt(np.exp(np.array([params[0],params[3],params[6]]))/2) #divide by 2 to get error on single spot!
    gamma = 3*np.exp(params[1])*data.delta_t**(-0.5)
    J = 3*np.exp(params[2])
    tau = (1/np.pi)*(J/gamma)**2

    return sigma,gamma,J,tau

