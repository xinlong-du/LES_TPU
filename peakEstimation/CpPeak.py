# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 15:57:29 2024

@author: xinlo
"""
import pandas as pd
from maxminest import maxminest

#%% input
recordDF=pd.read_csv('../postProcessing/LES_data/Cp/lesCpF0.csv',header=None);
record=recordDF.to_numpy().transpose()
dur_ratio=1.0;
plot_on=0;

#%% calculate peaks
max_est, min_est, max_std, min_std=maxminest(record,dur_ratio,plot_on);