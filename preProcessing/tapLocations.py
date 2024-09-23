# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 14:30:45 2024

@author: xinlo
"""

import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math

big_fig_size = (6,3);
sml_fig_size = (3,3);
plt_line_width = 0.5; 
fig_font_size = 8;

#load TPU wind tunnel test data
windDir=30; #deg
#tpuData=scipy.io.loadmat('../postprocessing/TPU_data/Cp_ts_g12060000.mat');
tpuData=scipy.io.loadmat('../postprocessing/TPU_data/Cp_ts_h12064500.mat');
tapCoord=tpuData['Location_of_measured_points'];
tapX=tapCoord[0,:];
tapY=tapCoord[1,:];

#plot coordinates of pressure taps
fig=plt.figure(figsize=sml_fig_size)
ax = fig.add_axes([0, 0, 1, 1])
ax.scatter(tapX,tapY,marker='.',s=25)
plt.rc('xtick', labelsize=fig_font_size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fig_font_size)    # fontsize of the tick labels
ax.tick_params(direction="in")
ax.set_xlabel('X (m)',fontsize=fig_font_size)
ax.set_ylabel('Y (m)',fontsize=fig_font_size)
plt.axis('equal')

#%% calculate 3D coordinates for pressure taps
"""
tapXYZface5=np.stack((tapX[0:96],tapY[0:96],np.full((96,),12.01)),axis=1);

tapXface1=tapX[96:126];
tapYface1=tapY[96:126];
tapXYZface1=np.stack((np.full((30,),-12.01),tapYface1,tapXface1+34),axis=1);

tapXface2=tapX[126:168];
tapYface2=tapY[126:168];
tapXYZface2=np.stack((tapXface2,np.full((42,),-8.01),tapYface2+24),axis=1);

tapXface3=tapX[168:198];
tapYface3=tapY[168:198];
tapXYZface3=np.stack((np.full((30,),12.01),tapYface3,34-tapXface3),axis=1);

tapXface4=tapX[198:];
tapYface4=tapY[198:];
tapXYZface4=np.stack((tapXface4,np.full((42,),8.01),24-tapYface4),axis=1);

tapXYZ=np.vstack((tapXYZface5,tapXYZface1,tapXYZface2,tapXYZface3,tapXYZface4));
tapXYZmodel=tapXYZ/100.0;

np.savetxt("flat_g12060000.csv",tapXYZmodel,delimiter=",");
"""
tapXYZface8=np.stack((tapX[0:28],tapY[0:28],20.01-tapY[0:28]),axis=1);

tapXface6=tapX[28:56];
tapYface6=tapY[28:56];
tapXYZface6=np.stack((tapXface6,tapYface6,20.01+tapYface6),axis=1);

tapXface5=tapX[56:76];
tapYface5=tapY[56:76];
tapXYZface5=np.stack((tapXface5,tapYface5,24.01+tapXface5),axis=1);

tapXface7=tapX[76:96];
tapYface7=tapY[76:96];
tapXYZface7=np.stack((tapXface7,tapYface7,24.01-tapXface7),axis=1);

tapXface1=tapX[96:126];
tapYface1=tapY[96:126];
tapXYZface1=np.stack((np.full((30,),-12.01),tapYface1,34+tapXface1),axis=1);

tapXface2=tapX[126:168];
tapYface2=tapY[126:168];
tapXYZface2=np.stack((tapXface2,np.full((42,),-8.01),tapYface2+24),axis=1);

tapXface3=tapX[168:198];
tapYface3=tapY[168:198];
tapXYZface3=np.stack((np.full((30,),12.01),tapYface3,34-tapXface3),axis=1);

tapXface4=tapX[198:];
tapYface4=tapY[198:];
tapXYZface4=np.stack((tapXface4,np.full((42,),8.01),24-tapYface4),axis=1);

tapXYZ=np.vstack((tapXYZface8,tapXYZface6,tapXYZface5,tapXYZface7,tapXYZface1,tapXYZface2,tapXYZface3,tapXYZface4));
tapXYZmodel=tapXYZ/100.0;

tapXrot= tapXYZmodel[:,0]*math.cos(windDir/180*math.pi)+tapXYZmodel[:,1]*math.sin(windDir/180*math.pi);
tapYrot=-tapXYZmodel[:,0]*math.sin(windDir/180*math.pi)+tapXYZmodel[:,1]*math.cos(windDir/180*math.pi);
tapXYZmodelRot=np.stack((tapXrot,tapYrot,tapXYZmodel[:,2]),axis=1);

np.savetxt("hip_h12064530.csv",tapXYZmodelRot,delimiter=",");