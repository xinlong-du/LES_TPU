# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 15:00:09 2024

@author: xinlo
"""
import numpy as np
import pandas as pd
import scipy as sp

#%% input
recordDF=pd.read_csv('../postProcessing/LES_data/Cp/lesCpF0.csv',header=None);
record=recordDF.to_numpy().transpose()
dur_ratio=1.0;
plot_on=1;

#%% setup
n_cdf_pk = 1000;     # number of points in CDF for peaks
cdf_pk_min = 0.0005; # minimum value of CDF for peaks
cdf_pk_max = 0.9995; # maximum value of CDF for peaks
cdf_pk = np.linspace(cdf_pk_min,cdf_pk_max, n_cdf_pk); # linearly spaced CDF values for peaks

rsize = np.size(record,0);

max_est = np.zeros([rsize,1]);
min_est = np.zeros([rsize,1]);

for i in range(0,rsize):
    x = record[i,:];
    n = len(x);

    mean_x = np.mean(x);
    std_x = np.std(x); # biased estimate
    skew_x = sum((x - mean_x)**3)/ (n*std_x**3); # biased estimate (Matlab default for 'skewness')
    
    X = x*np.sign(skew_x);  # Change sign (if necessary) to have positive skewness
    sort_X = np.sort(X); # sort in ascending order:
    
    mean_X = mean_x*np.sign(skew_x);
    std_X = std_x;
   
    CDF_X = np.linspace(1,n,n)/(n+1); # Empirical Cumulative Distribution Function
    
    # resample CDF more coarsely for more efficient parameter estimation:
    n_coarse = min(n,1000);
    CDF_coarse = np.linspace(1/(n_coarse+1),n_coarse/(n_coarse+1),n_coarse);
    temp = sp.interpolate.interp1d(CDF_X, sort_X);
    X_coarse = temp(CDF_coarse);
    
    # Estimate shape parameter of gamma distribution from coarse CDF:
    mean_X_coarse = np.mean(X_coarse);
    std_X_coarse = np.std(X_coarse);

    gamma_min = 1;
    gamma_max = 125;
    n_gamma = 19; # number of trial values of gamma
    n_start = 8; # start checking with this index in gamma_list
    gamma_list = np.logspace(np.log10(gamma_min),np.log10(gamma_max),n_gamma);
    gam_PPCC_list = np.zeros(np.size(gamma_list));
    count = 0;
    # first try decreasing gamma: