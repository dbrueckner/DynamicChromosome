import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import analysis
import data_read
import plotting

directory_data = "data"    
    
delta_t = 28 #in seconds
it = 3
s = 142

plt.close('all')
fig_size = [4,3]     
params = {
          'figure.figsize': fig_size,
          }
plt.rcParams.update(params)

data_csv = pd.read_csv(directory_data + '/data_line' + str(it) + '.csv')
data = data_read.transform_data_nostate(data_csv,delta_t,s)

N_t_max = 50


### distance distribution ###
plt.figure()
plt.hist(data.rij.ravel(),40)

plt.xlabel('distance')
plt.ylabel('probability density')
plt.tight_layout()


### single-locus MSD ###
MSD,MSD_all,MSD_err,time = analysis.corr_msd(data.X_blue,N_t_max,data.delta_t)

plt.figure()
plt.plot(time,MSD,'o')
plt.loglog()

plt.xlabel('time')
plt.ylabel('single-locus MSD')
plt.tight_layout()


### two-locus MSD ###

MSD,MSD_all,MSD_err,time = analysis.corr_msd(data.Xij,N_t_max,data.delta_t)

### two-locus MSD fitting ###
sigma,gamma,J,tau = analysis.fit_2loc_msd_rouse(data,mode_error='z')

plt.figure()
plt.plot(time,MSD,'o')

time = np.linspace(20,1400,100)
plt.plot(time,plotting.MSD_func(time,gamma,J,sigma),'-',color='r')

plt.loglog()

plt.xlabel('time')
plt.ylabel('two-locus MSD')
plt.tight_layout()
    
    
    
    
    
    
    
    
    
    
