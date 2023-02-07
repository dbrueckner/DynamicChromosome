#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: D.Brueckner
"""

import numpy as np


def MSD_func(time,gamma,J,sigma=np.zeros(3)): 
    import scipy.special
    tau = (1/np.pi)*(J/gamma)**2
    return 4*np.linalg.norm(sigma)**2 + 2*gamma*np.sqrt(time)*(1-np.exp(-tau/time)) + 2*J*scipy.special.erfc(np.sqrt(tau/time))
