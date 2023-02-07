#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: D.Brueckner
"""

import numpy as np

def transform_data(data_csv,delta_t,s):
    
    N_entries = len(data_csv.cell_id)
    N_tracks_total = int(max(data_csv.cell_id))+1
    N_t_max = int(max(data_csv.time_point))+1

    mode_all = ['blue','green','dist_vec','red','state']

    for mode in mode_all:
        if(mode == 'blue' or mode == 'green' or mode == 'dist_vec'):
            N_dim = 3
        else:
            N_dim = 1
            
        X = np.empty((N_tracks_total,N_t_max,N_dim))
        X[:] = np.nan

        for it in range(0,N_entries):
            j = int(data_csv.cell_id[it])
            t = int(data_csv.time_point[it])
            
            if(mode == 'blue'):
                X[j,t,0] = data_csv.x_blue[it]
                X[j,t,1] = data_csv.y_blue[it]
                X[j,t,2] = data_csv.z_blue[it]
            elif(mode == 'green'):
                X[j,t,0] = data_csv.x_green[it]
                X[j,t,1] = data_csv.y_green[it]
                X[j,t,2] = data_csv.z_green[it]
            elif(mode == 'dist_vec'):
                X[j,t,0] = data_csv.x_Rij[it]
                X[j,t,1] = data_csv.y_Rij[it]
                X[j,t,2] = data_csv.z_Rij[it]
            elif(mode == 'red'):
                X[j,t,0] = data_csv.red[it]
            elif(mode == 'state'):
                X[j,t,0] = data_csv.state[it]

        if(mode == 'blue'):
            X_blue = X
        elif(mode == 'green'):
            X_green = X
        elif(mode == 'dist_vec'):
            X_dist_vec = X
        elif(mode == 'red'):
            X_red = X
        elif(mode == 'state'):
            X_state = np.round(X)
    
    from fns_data_wrapper import StochasticTrajectoryData
    data = StochasticTrajectoryData(X_blue,X_green,X_dist_vec,X_red,X_state,delta_t,s)
    
    return data
    
def transform_data_nostate(data_csv,delta_t,s):
    
    N_entries = len(data_csv.cell_id)
    N_tracks_total = int(max(data_csv.cell_id))+1
    N_t_max = int(max(data_csv.time_point))+1

    mode_all = ['blue','green','dist_vec','red']

    for mode in mode_all:
        if(mode == 'blue' or mode == 'green' or mode == 'dist_vec'):
            N_dim = 3
        else:
            N_dim = 1
            
        X = np.empty((N_tracks_total,N_t_max,N_dim))
        X[:] = np.nan

        for it in range(0,N_entries):
            j = int(data_csv.cell_id[it])
            t = int(data_csv.time_point[it])
            
            if(mode == 'blue'):
                X[j,t,0] = data_csv.x_blue[it]
                X[j,t,1] = data_csv.y_blue[it]
                X[j,t,2] = data_csv.z_blue[it]
            elif(mode == 'green'):
                X[j,t,0] = data_csv.x_green[it]
                X[j,t,1] = data_csv.y_green[it]
                X[j,t,2] = data_csv.z_green[it]
            elif(mode == 'dist_vec'):
                X[j,t,0] = data_csv.x_Rij[it]
                X[j,t,1] = data_csv.y_Rij[it]
                X[j,t,2] = data_csv.z_Rij[it]
            elif(mode == 'red'):
                X[j,t,0] = data_csv.red[it]

        if(mode == 'blue'):
            X_blue = X
        elif(mode == 'green'):
            X_green = X
        elif(mode == 'dist_vec'):
            X_dist_vec = X
        elif(mode == 'red'):
            X_red = X
    
    X_state = None
    from fns_data_wrapper import StochasticTrajectoryData
    data = StochasticTrajectoryData(X_blue,X_green,X_dist_vec,X_red,X_state,delta_t,s)
    
    return data



def transform_data_hansen(data_csv,delta_t,s):
    
    N_entries = len(data_csv.id)
    N_tracks_total = int(max(data_csv.id))+1
    N_t_max = int(max(data_csv.t))+1

    mode_all = ['blue','green','dist_vec']

    for mode in mode_all:
        if(mode == 'blue' or mode == 'green' or mode == 'dist_vec'):
            N_dim = 3
        else:
            N_dim = 1
            
        X = np.empty((N_tracks_total,N_t_max,N_dim))
        X[:] = np.nan

        for it in range(0,N_entries):
            j = int(data_csv.id[it])
            t = int(data_csv.t[it])
            
            if(mode == 'blue'):
                X[j,t,0] = data_csv.x[it]
                X[j,t,1] = data_csv.y[it]
                X[j,t,2] = data_csv.z[it]
            elif(mode == 'green'):
                X[j,t,0] = data_csv.x2[it]
                X[j,t,1] = data_csv.y2[it]
                X[j,t,2] = data_csv.z2[it]
            elif(mode == 'dist_vec'):
                X[j,t,0] = data_csv.x[it] - data_csv.x2[it]
                X[j,t,1] = data_csv.y[it] - data_csv.y2[it]
                X[j,t,2] = data_csv.z[it] - data_csv.z2[it]


        if(mode == 'blue'):
            X_blue = X*1e3
        elif(mode == 'green'):
            X_green = X*1e3
        elif(mode == 'dist_vec'):
            X_dist_vec = X*1e3
            
    X_red = None
    X_state = None
    from fns_data_wrapper import StochasticTrajectoryData
    data = StochasticTrajectoryData(X_blue,X_green,X_dist_vec,X_red,X_state,delta_t,s)
    
    return data

