# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 15:57:29 2024

@author: xinlo
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io
from maxminest import maxminest

sml_fig_size = (3,3); 
fig_font_size = 8;

#%% input
dur_ratio=6.0; #ratio between the duration of the wanted peak and the duration of the data
plot_on=0;     #plots in the peak estimation process
dirs=['0','15','30','45','60','75','90'];

#flat
for dir1 in dirs:
    CpLES_DF=pd.read_csv('../postProcessing/LES_data/Cp/lesCpF'+dir1+'.csv',header=None);
    CpLES=CpLES_DF.to_numpy().transpose();
    maxCpLES, minCpLES, max_std, min_std=maxminest(CpLES,dur_ratio,plot_on);
    peakCpLES=np.hstack((maxCpLES,minCpLES));
    np.savetxt("../postProcessing/LES_data/Cp/lesCpF"+dir1+"peak.csv",peakCpLES,fmt='%2.6f',delimiter=",");
    
    if dir1=='0':
        CpTPU_DT=scipy.io.loadmat('../postProcessing/TPU_data/Cp_ts_g120600'+dir1+'0.mat'); 
    else:
        CpTPU_DT=scipy.io.loadmat('../postProcessing/TPU_data/Cp_ts_g120600'+dir1+'.mat');
    CpTPU=CpTPU_DT['Wind_pressure_coefficients'].transpose();
    maxCpTPU, minCpTPU, max_std, min_std=maxminest(CpTPU,dur_ratio,plot_on);
    peakCpTPU=np.hstack((maxCpTPU,minCpTPU));
    np.savetxt("../postProcessing/LES_data/Cp/tpuCpF"+dir1+"peak.csv",peakCpTPU,fmt='%2.6f',delimiter=",");
    
    fig=plt.figure(figsize=sml_fig_size)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.scatter(maxCpLES,maxCpTPU,marker='.',s=10,label='Max')
    ax.scatter(minCpLES,minCpTPU,marker='.',s=10,label='Min')
    ax.plot([-6,6],[-6,6],'k-');
    ax.plot([-6,6],[-6/1.1,6/1.1],'b-',label='+/-10%',linewidth=0.5);
    ax.plot([-6,6],[-6/0.9,6/0.9],'b-',linewidth=0.5);
    ax.plot([-6,6],[-6/1.2,6/1.2],'r-',label='+/-20%',linewidth=0.5);
    ax.plot([-6,6],[-6/0.8,6/0.8],'r-',linewidth=0.5);
    ax.plot([-6,6],[-6/1.3,6/1.3],'m-',label='+/-30%',linewidth=0.5);
    ax.plot([-6,6],[-6/0.7,6/0.7],'m-',linewidth=0.5);
    plt.rc('xtick', labelsize=fig_font_size)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=fig_font_size)    # fontsize of the tick labels
    ax.tick_params(direction="in")
    ax.set_xlabel('Max/min of Cp from LES',fontsize=fig_font_size)
    ax.set_ylabel('Max/min of Cp from TPU',fontsize=fig_font_size)
    # plt.xlim(-10,6);
    # plt.ylim(-10,6);
    plt.axis('equal')
    plt.legend(prop={'size': fig_font_size})
    figName='../postProcessing/LES_data/Cp/CompPeakCpF'+dir1+'.tif'
    plt.savefig(figName, transparent=False, bbox_inches='tight', dpi=100)
    plt.show()
    
#gable
for dir1 in dirs:
    CpLES_DF=pd.read_csv('../postProcessing/LES_data/Cp/lesCpG'+dir1+'.csv',header=None);
    CpLES=CpLES_DF.to_numpy().transpose();
    maxCpLES, minCpLES, max_std, min_std=maxminest(CpLES,dur_ratio,plot_on);
    peakCpLES=np.hstack((maxCpLES,minCpLES));
    np.savetxt("../postProcessing/LES_data/Cp/lesCpG"+dir1+"peak.csv",peakCpLES,fmt='%2.6f',delimiter=",");
    
    if dir1=='0':
        CpTPU_DT=scipy.io.loadmat('../postProcessing/TPU_data/Cp_ts_g120422'+dir1+'0.mat'); 
    else:
        CpTPU_DT=scipy.io.loadmat('../postProcessing/TPU_data/Cp_ts_g120422'+dir1+'.mat');
    CpTPU=CpTPU_DT['Wind_pressure_coefficients'].transpose();
    maxCpTPU, minCpTPU, max_std, min_std=maxminest(CpTPU,dur_ratio,plot_on);
    peakCpTPU=np.hstack((maxCpTPU,minCpTPU));
    np.savetxt("../postProcessing/LES_data/Cp/tpuCpG"+dir1+"peak.csv",peakCpTPU,fmt='%2.6f',delimiter=",");
    
    fig=plt.figure(figsize=sml_fig_size)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.scatter(maxCpLES,maxCpTPU,marker='.',s=10,label='Max')
    ax.scatter(minCpLES,minCpTPU,marker='.',s=10,label='Min')
    ax.plot([-6,6],[-6,6],'k-');
    ax.plot([-6,6],[-6/1.1,6/1.1],'b-',label='+/-10%',linewidth=0.5);
    ax.plot([-6,6],[-6/0.9,6/0.9],'b-',linewidth=0.5);
    ax.plot([-6,6],[-6/1.2,6/1.2],'r-',label='+/-20%',linewidth=0.5);
    ax.plot([-6,6],[-6/0.8,6/0.8],'r-',linewidth=0.5);
    ax.plot([-6,6],[-6/1.3,6/1.3],'m-',label='+/-30%',linewidth=0.5);
    ax.plot([-6,6],[-6/0.7,6/0.7],'m-',linewidth=0.5);
    plt.rc('xtick', labelsize=fig_font_size)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=fig_font_size)    # fontsize of the tick labels
    ax.tick_params(direction="in")
    ax.set_xlabel('Max/min of Cp from LES',fontsize=fig_font_size)
    ax.set_ylabel('Max/min of Cp from TPU',fontsize=fig_font_size)
    # plt.xlim(-10,6);
    # plt.ylim(-10,6);
    plt.axis('equal')
    plt.legend(prop={'size': fig_font_size})
    figName='../postProcessing/LES_data/Cp/CompPeakCpG'+dir1+'.tif'
    plt.savefig(figName, transparent=False, bbox_inches='tight', dpi=100)
    plt.show()
    
#hip
for dir1 in dirs:
    CpLES_DF=pd.read_csv('../postProcessing/LES_data/Cp/lesCpH'+dir1+'.csv',header=None);
    CpLES=CpLES_DF.to_numpy().transpose();
    maxCpLES, minCpLES, max_std, min_std=maxminest(CpLES,dur_ratio,plot_on);
    peakCpLES=np.hstack((maxCpLES,minCpLES));
    np.savetxt("../postProcessing/LES_data/Cp/lesCpH"+dir1+"peak.csv",peakCpLES,fmt='%2.6f',delimiter=",");
    
    if dir1=='0':
        CpTPU_DT=scipy.io.loadmat('../postProcessing/TPU_data/Cp_ts_h120645'+dir1+'0.mat'); 
    else:
        CpTPU_DT=scipy.io.loadmat('../postProcessing/TPU_data/Cp_ts_h120645'+dir1+'.mat');
    CpTPU=CpTPU_DT['Wind_pressure_coefficients'].transpose();
    maxCpTPU, minCpTPU, max_std, min_std=maxminest(CpTPU,dur_ratio,plot_on);
    peakCpTPU=np.hstack((maxCpTPU,minCpTPU));
    np.savetxt("../postProcessing/LES_data/Cp/tpuCpH"+dir1+"peak.csv",peakCpTPU,fmt='%2.6f',delimiter=",");
    
    fig=plt.figure(figsize=sml_fig_size)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.scatter(maxCpLES,maxCpTPU,marker='.',s=10,label='Max')
    ax.scatter(minCpLES,minCpTPU,marker='.',s=10,label='Min')
    ax.plot([-6,6],[-6,6],'k-');
    ax.plot([-6,6],[-6/1.1,6/1.1],'b-',label='+/-10%',linewidth=0.5);
    ax.plot([-6,6],[-6/0.9,6/0.9],'b-',linewidth=0.5);
    ax.plot([-6,6],[-6/1.2,6/1.2],'r-',label='+/-20%',linewidth=0.5);
    ax.plot([-6,6],[-6/0.8,6/0.8],'r-',linewidth=0.5);
    ax.plot([-6,6],[-6/1.3,6/1.3],'m-',label='+/-30%',linewidth=0.5);
    ax.plot([-6,6],[-6/0.7,6/0.7],'m-',linewidth=0.5);
    plt.rc('xtick', labelsize=fig_font_size)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=fig_font_size)    # fontsize of the tick labels
    ax.tick_params(direction="in")
    ax.set_xlabel('Max/min of Cp from LES',fontsize=fig_font_size)
    ax.set_ylabel('Max/min of Cp from TPU',fontsize=fig_font_size)
    # plt.xlim(-10,6);
    # plt.ylim(-10,6);
    plt.axis('equal')
    plt.legend(prop={'size': fig_font_size})
    figName='../postProcessing/LES_data/Cp/CompPeakCpH'+dir1+'.tif'
    plt.savefig(figName, transparent=False, bbox_inches='tight', dpi=100)
    plt.show()