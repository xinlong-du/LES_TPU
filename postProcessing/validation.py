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
#tpuData=scipy.io.loadmat('./TPU_data/Cp_ts_g12042290.mat'); #gable
#tpuData=scipy.io.loadmat('./TPU_data/Cp_ts_g12060000.mat'); #flat
tpuData=scipy.io.loadmat('./TPU_data/Cp_ts_h12064515.mat');  #hip

tpuCp=tpuData['Wind_pressure_coefficients'];

#load LES data
#gable
# lesData1=np.loadtxt('./LES_data/Gable/0/p');
# lesData2=np.loadtxt('./LES_data/Gable/7.6/p');
# lesData3=np.loadtxt('./LES_data/Gable/17.1/p');
# lesData=np.vstack((lesData1[0:3800,:],lesData2[0:4750,:],lesData3));

#flat
# lesData1=np.loadtxt('./LES_data/Flat/0/p');
# lesData2=np.loadtxt('./LES_data/Flat/11.4/p');
# lesData=np.vstack((lesData1[0:5700,:],lesData2));

#hip
lesData=np.loadtxt('./LES_data/Hip/p');

lesP=lesData[501:,1:];
time=lesData[1:9001,0];

U=7.14822; #reference wind speed (m/s)
lesCp=lesP/(0.5*U*U); #The pressure is kinematic pressure pk=ps/rho (m^2/s^2)

#%% plot coordinates of pressure taps
tapCoord=tpuData['Location_of_measured_points'];
tapX=tapCoord[0,:];
tapY=tapCoord[1,:];

fig=plt.figure(figsize=sml_fig_size)
ax = fig.add_axes([0, 0, 1, 1])
ax.scatter(tapX,tapY,marker='.',s=25)
plt.rc('xtick', labelsize=fig_font_size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fig_font_size)    # fontsize of the tick labels
ax.tick_params(direction="in")
ax.set_xlabel('X (m)',fontsize=fig_font_size)
ax.set_ylabel('Y (m)',fontsize=fig_font_size)
plt.axis('equal')

#%% compare TPU and LES data: time series
for i in range(0,5):
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
tapID=np.arange(np.size(tpuCp,1));

fig=plt.figure(figsize=big_fig_size)
ax = fig.add_axes([0, 0, 1, 1])
lesPlot,=ax.plot(tapID,meanLesCp,linewidth=plt_line_width,color="k")
tpuPlot,=ax.plot(tapID,meanTpuCp,linewidth=plt_line_width,color="r")
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

#contour plot of mean Cp
fig,ax=plt.subplots(1, 1) 
[X,Y]=np.meshgrid(tapX,tapY)
plt.plot(X,Y, marker='.', color='k', linestyle='none')

fig,ax=plt.subplots(1, 1) 
t=ax.tricontourf(tapX, tapY, meanTpuCp) 
ax.set_title('TPU mean Cp') 
plt.colorbar(t)
plt.show()

fig,ax=plt.subplots(1, 1) 
t=ax.tricontourf(tapX, tapY, meanLesCp) 
ax.set_title('LES mean Cp')  
plt.colorbar(t)
plt.show()

#%% compare TPU and LES data: std
stdTpuCp=np.std(tpuCp[0:len(time),:],axis=0);
stdLesCp=np.std(lesCp,axis=0);

fig=plt.figure(figsize=big_fig_size)
ax = fig.add_axes([0, 0, 1, 1])
lesPlot,=ax.plot(tapID,stdLesCp,linewidth=plt_line_width,color="k")
tpuPlot,=ax.plot(tapID,stdTpuCp,linewidth=plt_line_width,color="r")
plt.rc('xtick', labelsize=fig_font_size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fig_font_size)    # fontsize of the tick labels
ax.tick_params(direction="in")
ax.set_xlabel('Tap ID',fontsize=fig_font_size)
ax.set_ylabel('Std of Cp',fontsize=fig_font_size)
ax.legend([lesPlot,tpuPlot],['LES','TPU'],prop={'size': fig_font_size})

fig=plt.figure(figsize=sml_fig_size)
ax = fig.add_axes([0, 0, 1, 1])
ax.scatter(stdLesCp,stdTpuCp,marker='.',s=25)
ax.plot([0.1,0.6],[0.1,0.6],'r-');
plt.rc('xtick', labelsize=fig_font_size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fig_font_size)    # fontsize of the tick labels
ax.tick_params(direction="in")
ax.set_xlabel('Std of Cp from LES',fontsize=fig_font_size)
ax.set_ylabel('Std of Cp from TPU',fontsize=fig_font_size)
plt.xlim(0.1,0.6);
plt.ylim(0.1,0.6);
plt.axis('equal')

#%% compare TPU and LES data: peak
maxTpuCp=np.max(tpuCp[0:len(time),:],axis=0);
maxLesCp=np.max(lesCp,axis=0);
minTpuCp=np.min(tpuCp[0:len(time),:],axis=0);
minLesCp=np.min(lesCp,axis=0);

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

#contour plot of peak Cp
fig,ax=plt.subplots(1, 1) 
t=ax.tricontourf(tapX, tapY, maxTpuCp, levels=np.linspace(-0.5,5,12)) 
ax.set_title('TPU max Cp') 
plt.colorbar(t)
plt.show()

fig,ax=plt.subplots(1, 1) 
t=ax.tricontourf(tapX, tapY, maxLesCp, levels=np.linspace(-0.5,5,12)) 
ax.set_title('LES max Cp')  
plt.colorbar(t)
plt.show()

fig,ax=plt.subplots(1, 1) 
t=ax.tricontourf(tapX, tapY, minTpuCp, levels=np.linspace(-7.0,0.0,8)) 
ax.set_title('TPU min Cp') 
plt.colorbar(t)
plt.show()

fig,ax=plt.subplots(1, 1) 
t=ax.tricontourf(tapX, tapY, minLesCp, levels=np.linspace(-7.0,0.0,8)) 
ax.set_title('LES min Cp')  
plt.colorbar(t)
plt.show()

#%% compare TPU and LES data: histogram
for i in range(0,5):
    fig=plt.figure(figsize=big_fig_size)
    plt.hist(lesCp[0:,i], density=True, bins=30, alpha = 0.5, range=(-2.5,0.5), edgecolor='black')  # density=False would make counts
    plt.hist(tpuCp[0:len(time),i], density=True, bins=30, alpha = 0.5, range=(-2.5,0.5), edgecolor='black') 
    plt.ylabel('PDF')
    plt.xlabel('Data')
    plt.legend(['LES', 'TPU'])