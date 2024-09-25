# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 12:14:24 2024

@author: xinlo
"""

import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

U=7.14822; #reference wind speed (m/s)

#%% hip
#dirs=['0','15','30','45','60','75','90'];
dirs=['15']
for dir1 in dirs:
    lesData=np.loadtxt('./LES_data/Hip/p'+dir1);
    lesP=lesData[501:,1:];
    time=lesData[1:9001,0];
    lesCp=lesP/(0.5*U*U); #The pressure is kinematic pressure pk=ps/rho (m^2/s^2)
    np.savetxt("./LES_data/Cp/lesCpH"+dir1+".csv",lesCp,delimiter=",");

#%% compare TPU and LES data: peak
tpuData=scipy.io.loadmat('./TPU_data/Cp_ts_h12064515.mat');  #hip
tpuCp=tpuData['Wind_pressure_coefficients'];

tapID=np.arange(np.size(tpuCp,1));
maxTpuCp=np.max(tpuCp[0:len(time),:],axis=0);
maxLesCp=np.max(lesCp,axis=0);
minTpuCp=np.min(tpuCp[0:len(time),:],axis=0);
minLesCp=np.min(lesCp,axis=0);

big_fig_size = (6,3);
sml_fig_size = (3,3);
plt_line_width = 0.5; 
fig_font_size = 8;

fig=plt.figure(figsize=big_fig_size)
ax = fig.add_axes([0, 0, 1, 1])
maxLesPlot,=ax.plot(tapID,maxLesCp,linewidth=plt_line_width,color="k")
maxTpuPlot,=ax.plot(tapID,maxTpuCp,linewidth=plt_line_width,color="r")
minLesPlot,=ax.plot(tapID,minLesCp,linewidth=plt_line_width,color="k")
minTpuPlot,=ax.plot(tapID,minTpuCp,linewidth=plt_line_width,color="r")
plt.rc('xtick', labelsize=fig_font_size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fig_font_size)    # fontsize of the tick labels
ax.tick_params(direction="in")
ax.set_xlabel('Tap ID',fontsize=fig_font_size)
ax.set_ylabel('Cp',fontsize=fig_font_size)
ax.legend([maxLesPlot,maxTpuPlot],['LES max (min)','TPU max (min)'],prop={'size': fig_font_size})

fig=plt.figure(figsize=sml_fig_size)
ax = fig.add_axes([0, 0, 1, 1])
ax.scatter(maxLesCp,maxTpuCp,marker='.',s=25)
ax.scatter(minLesCp,minTpuCp,marker='.',s=25)
ax.plot([-10,5],[-10,5],'r-');
plt.rc('xtick', labelsize=fig_font_size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fig_font_size)    # fontsize of the tick labels
ax.tick_params(direction="in")
ax.set_xlabel('Max/min of Cp from LES',fontsize=fig_font_size)
ax.set_ylabel('Max/min of Cp from TPU',fontsize=fig_font_size)
plt.xlim(-10,5);
plt.ylim(-10,5);
plt.axis('equal')
ax.legend(['Max','Min'],prop={'size': fig_font_size})
