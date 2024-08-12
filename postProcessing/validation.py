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
sml_fig_size = (3,3);
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

U=7.14822; #reference wind speed (m/s)
lesCp=lesP/(0.5*U*U); #The pressure is kinematic pressure pk=ps/rho (m^2/s^2)

#%% compare TPU and LES data: time series
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
    
#%% compare TPU and LES data: mean Cp
meanTpuCp=np.mean(tpuCp[0:len(time),:],axis=0);
meanLesCp=np.mean(lesCp,axis=0);
tapID=np.arange(200);

fig=plt.figure(figsize=big_fig_size)
ax = fig.add_axes([0, 0, 1, 1])
lesPlot,=ax.plot(tapID,meanLesCp,linewidth=plt_line_width)
tpuPlot,=ax.plot(tapID,meanTpuCp,linewidth=plt_line_width)
plt.rc('xtick', labelsize=fig_font_size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fig_font_size)    # fontsize of the tick labels
ax.tick_params(direction="in")
ax.set_xlabel('Tap ID',fontsize=fig_font_size)
ax.set_ylabel('mean Cp',fontsize=fig_font_size)
ax.legend([lesPlot,tpuPlot],['LES','TPU'],prop={'size': fig_font_size})

fig=plt.figure(figsize=sml_fig_size)
ax = fig.add_axes([0, 0, 1, 1])
ax.scatter(meanLesCp,meanTpuCp,marker='.',s=25)
ax.plot([-1.25,1],[-1.25,1],'r-');
plt.rc('xtick', labelsize=fig_font_size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fig_font_size)    # fontsize of the tick labels
ax.tick_params(direction="in")
ax.set_xlabel('Mean Cp from LES',fontsize=fig_font_size)
ax.set_ylabel('Mean Cp from TPU',fontsize=fig_font_size)
plt.xlim(-1.25,1);
plt.ylim(-1.25,1);
plt.axis('equal')

#%% compare TPU and LES data: histogram
for i in range(0,9):
    fig=plt.figure(figsize=big_fig_size)
    plt.hist(lesCp[0:,i], density=True, bins=30, alpha = 0.5, range=(-2.5,0.5), edgecolor='black')  # density=False would make counts
    plt.hist(tpuCp[0:len(time),i], density=True, bins=30, alpha = 0.5, range=(-2.5,0.5), edgecolor='black') 
    plt.ylabel('PDF')
    plt.xlabel('Data')
    plt.legend(['LES', 'TPU'])