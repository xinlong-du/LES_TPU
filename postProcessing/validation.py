# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 15:21:41 2024

@author: xinlo
"""

import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

big_fig_size = (6,3);
plt_line_width = 0.5; 
fig_font_size = 8;

#%% load data
#load TPU wind tunnel test data
tpuData=scipy.io.loadmat('./TPU_data/Cp_ts_g12042290.mat');
tpuCp=tpuData['Wind_pressure_coefficients'];

#load LES data
lesData=np.loadtxt('./LES_data/p');
lesP=lesData[1:,1:];
time=lesData[1:,0];

U=7.4; #reference wind speed (m/s)
lesCp=lesP/(0.5*U*U); #The pressure is kinematic pressure pk=ps/rho (m^2/s^2)

#%% compare TPU and LES data
for i in range(0,9):
    fig=plt.figure(figsize=big_fig_size)
    ax = fig.add_axes([0, 0, 1, 1])
    lesPlot,=ax.plot(time,lesCp[0:,i], linewidth=plt_line_width)
    tpuPlot,=ax.plot(time,tpuCp[0:len(time),i], linewidth=plt_line_width)
    plt.rc('xtick', labelsize=fig_font_size)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=fig_font_size)    # fontsize of the tick labels
    ax.tick_params(direction="in")
    ax.set_xlabel('Time (s)',fontsize=fig_font_size)
    ax.set_ylabel('Cp',fontsize=fig_font_size)
    ax.legend([lesPlot,tpuPlot],['LES','TPU'],prop={'size': fig_font_size})