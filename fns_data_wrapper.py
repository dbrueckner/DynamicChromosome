#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 15:52:27 2020

@author: D.Brueckner
"""

import numpy as np

class StochasticTrajectoryData(object):
 
    def __init__(self,X_blue,X_green,X_dist_vec,X_red,X_state,delta_t,s):
        
        self.delta_t = delta_t
        self.s = s
        
        self.X_blue = X_blue
        self.X_green = X_green
        
        self.X_blue[np.isnan(self.X_green)] = np.nan
        self.X_green[np.isnan(self.X_blue)] = np.nan
        
        indices_nonan_particles = ~np.all(np.isnan(self.X_blue[:,:,0]), axis=1) #particles that have non-nan entries
        
        self.X_blue = self.X_blue[indices_nonan_particles,:,:]
        self.X_green = self.X_green[indices_nonan_particles,:,:]
        
        self.Xij = X_dist_vec[indices_nonan_particles,:,:]
        self.rij = np.sqrt(self.Xij[:,:,0]**2+self.Xij[:,:,1]**2+self.Xij[:,:,2]**2)
        
        if X_red is not None:
            self.X_red = X_red[indices_nonan_particles,:,:]
        else:
            self.X_red = None
            
        if X_state is not None:    
            self.X_state = X_state[indices_nonan_particles,:,:]
        else:
            self.X_state = None
        
        self.r_blue = np.sqrt(self.X_blue[:,:,0]**2+self.X_blue[:,:,1]**2+self.X_blue[:,:,2]**2)
        self.r_green = np.sqrt(self.X_green[:,:,0]**2+self.X_green[:,:,1]**2+self.X_green[:,:,2]**2)
        self.X_com = (self.X_blue+self.X_green)/2

        self.N_part = self.X_blue.shape[0]
        self.N_t_max = self.X_blue.shape[1]
        self.N_t = []
        for j in range(0,self.N_part):
            length = (~np.isnan(self.X_blue[j,:,0])).cumsum(0).argmax(0)
            if(length==0):
                self.N_t.append(length)
            else:
                self.N_t.append(length+1)
        
        self.t_start = []
        self.t_end = []
        self.No_t = []
        for j in range(0,self.N_part):
            non_nan_indices = [i for i, x in enumerate(~np.isnan(self.X_blue[j,:,0])) if x]
            if len(non_nan_indices) > 0:
                self.t_start.append(non_nan_indices[0])
                self.t_end.append(non_nan_indices[-1])
                self.No_t.append(non_nan_indices[-1]-non_nan_indices[0]+1)
            else:
                self.t_start.append(np.nan)
                self.t_end.append(np.nan)
                self.No_t.append(np.nan)
 