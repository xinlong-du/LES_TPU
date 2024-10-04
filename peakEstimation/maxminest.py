# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 15:00:09 2024

@author: xinlo
"""
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
from stdgaminv import stdgaminv

#%% input
recordDF=pd.read_csv('../postProcessing/LES_data/Cp/lesCpF0.csv',header=None);
record=recordDF.to_numpy().transpose()
dur_ratio=1.0;
plot_on=0;

#%% setup
n_cdf_pk = 1000;     # number of points in CDF for peaks
cdf_pk_min = 0.0005; # minimum value of CDF for peaks
cdf_pk_max = 0.9995; # maximum value of CDF for peaks
cdf_pk = np.linspace(cdf_pk_min,cdf_pk_max, n_cdf_pk); # linearly spaced CDF values for peaks

rsize = np.size(record,0);

max_est = np.zeros([rsize,1]);
min_est = np.zeros([rsize,1]);
max_std = np.zeros([rsize,1]);
min_std = np.zeros([rsize,1]);

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
    beta_coarse_list=np.zeros(n_gamma);
    mu_coarse_list=np.zeros(n_gamma);
    for j in range(n_start-1,-1,-1):
        count = count+1;
        # Obtain the Gamma Distribution Parameters for current gamma:
        s_gam_j = stdgaminv(CDF_coarse, gamma_list[j]); # standard variate
        mean_s_gam_j = np.mean(s_gam_j);
        # linear regression:
        beta_coarse_list[j] = (sum(s_gam_j*X_coarse)-n_coarse*mean_s_gam_j*mean_X_coarse)/(sum(s_gam_j**2)-n_coarse*mean_s_gam_j**2);
        mu_coarse_list[j] = mean_X_coarse - beta_coarse_list[j]*mean_s_gam_j;
        # Probability Plot Correlation Coefficient:
        gam_PPCC_list[j] = beta_coarse_list[j]*np.std(s_gam_j)/std_X_coarse;
        X_coarse_fit_j = mu_coarse_list[j] + beta_coarse_list[j]*s_gam_j;
        if plot_on:
            plt.figure();
            plt.plot(s_gam_j,X_coarse,'.',s_gam_j,X_coarse_fit_j,'-');
            plt.title('gamma: ');
            plt.show();
        if gam_PPCC_list[j]==max(gam_PPCC_list):
            gam = gamma_list[j];
            gam_PPCC_max = gam_PPCC_list[j];
        else:
            break; # stop searching once the PPCC starts to decrease
    if gam_PPCC_list[n_start-2]<gam_PPCC_list[n_start-1]:
        # if the PPCC decreased with decreasing gamma, try increasing gamma: 
        for j in range(n_start,n_gamma):
            count = count+1;
            # Obtain the Gamma Distribution Parameters for current gamma:
            s_gam_j = stdgaminv(CDF_coarse, gamma_list[j]); # standard variate
            mean_s_gam_j = np.mean(s_gam_j);
            # linear regression:
            beta_coarse_list[j] = (sum(s_gam_j*X_coarse)-n_coarse*mean_s_gam_j*mean_X_coarse)/(sum(s_gam_j**2)-n_coarse*mean_s_gam_j**2);
            mu_coarse_list[j] = mean_X_coarse - beta_coarse_list[j]*mean_s_gam_j;
            # Probability Plot Correlation Coefficient:
            gam_PPCC_list[j] = beta_coarse_list[j]*np.std(s_gam_j)/std_X_coarse;
            X_coarse_fit_j = mu_coarse_list[j] + beta_coarse_list[j]*s_gam_j;
            if plot_on:
                plt.figure()
                plt.plot(s_gam_j,X_coarse,'.',s_gam_j,X_coarse_fit_j,'-');
                plt.title('gamma: ');
                plt.show();
            if gam_PPCC_list[j]==max(gam_PPCC_list):
                gam = gamma_list[j];
                gam_PPCC_max = gam_PPCC_list[j];
            else:
                break; # stop searching once the PPCC starts to decrease
                
    # Obtain the Gamma Distribution Parameters for the best-fit gamma using all the data:
    s_gam = stdgaminv(CDF_X, gam); # standard variate
    mean_s_gam = np.mean(s_gam);
    
    beta = (sum(s_gam*sort_X)-n*mean_s_gam*mean_X)/(sum(s_gam**2)-n*mean_s_gam**2);
    mu = mean_X - beta*mean_s_gam;
    gam_PPCC = beta*np.std(s_gam)/std_X;

    X_fit = mu + beta*s_gam;
    if plot_on:
        plt.figure()
        plt.plot(gamma_list[np.where(gam_PPCC_list)],gam_PPCC_list[np.where(gam_PPCC_list)],'.',gam,gam_PPCC_max,'o')
        plt.xlabel('Shape parameter, gamma');
        plt.ylabel('PPCC');
        plt.show()
        
        plt.figure()
        plt.plot(s_gam,sort_X,'.',s_gam,X_fit,'-');
        plt.xlabel('Standard variate for gamma distribution');
        plt.ylabel('Value of time series');
        plt.title('gamma: ');
        plt.legend(['All data','Best-fit gamma distribution'])
        plt.show()
        
    # Obtain the Normal Distribution Parameters for lower portion of CDF

    CDF_split = .25; # fit normal distribution to data below this cumulative probability
    temp = sp.interpolate.interp1d(CDF_X, sort_X, );
    X_split = temp(CDF_split);
    ind_low = np.where(sort_X<X_split);
    X_low = sort_X[ind_low];
    n_low = len(X_low);
    CDF_low = CDF_X[ind_low];

    s_norm_low = -np.sqrt(2)*sp.special.erfcinv(2*CDF_low); # standard normal variate
    mean_s_norm_low = np.mean(s_norm_low);
    mean_X_low = np.mean(X_low);
    # linear regression:
    sigma_low = (sum(s_norm_low*X_low)-n_low*mean_s_norm_low*mean_X_low)/(sum(s_norm_low**2)-n_low*mean_s_norm_low**2);
    mu_low = mean_X_low - sigma_low*mean_s_norm_low;
    X_low_fit = mu_low + sigma_low*s_norm_low;
    # Probability Plot Correlation Coefficient:
    norm_PPCC = sigma_low*np.std(s_norm_low)/np.std(X_low);

    if plot_on:
        plt.figure()
        plt.plot(s_norm_low,X_low,'.',s_norm_low,X_low_fit,'-');
        plt.xlabel('Standard normal variate');
        plt.ylabel('Value of time series');
        plt.title('Normal distribution fit to lower tail');
        plt.legend(['Lower tail data','Best-fit Normal distribution'])
        plt.show()
        
        plt.figure();
        plt.plot(sort_X,CDF_X,'k.',X_fit,CDF_X,'r-',X_low_fit,CDF_low,'g-');
        plt.xlabel('Value of time series');
        plt.ylabel('Cumulative Probability');
        plt.legend(['Empirical CDF: all data','Gamma distribution fit to all data','Normal distribution fit to lower tail']);
        plt.show()
   
    # Compute the mean zero upcrossing rate of process y(t) with standardized
    # normal probability distribution.
    X_u = np.mean(sort_X[ np.where(abs(CDF_X - 0.5) == min(abs(CDF_X - 0.5))) ]); # upcrossing level
    upcross=np.where( (X[1:]>=X_u) & (X[0:-1]<X_u) );
    Nupcross = len(upcross[0]); # number of upcrossings
    if Nupcross<100:
        print('The number of median upcrossings is low. The record may be too short for accurate peak estimation.');
    y_pk = np.sqrt(2.0*np.log(-dur_ratio*Nupcross / np.log(cdf_pk))); # maximum peaks corresponding to specified cumulative probabilities
    CDF_y = 0.5*sp.special.erfc(-y_pk/np.sqrt(2));

    # Perform the mapping procedure to compute the CDF of largest peak
    # for X(t) from y(t)
    X_max = stdgaminv(CDF_y,gam)*beta + mu;
    X_min = -np.sqrt(2)*sp.special.erfcinv(2*(1-CDF_y))*sigma_low + mu_low;
    
    # Compute the Mean of the Peaks for process X(t)
    pdf_pk = -y_pk * cdf_pk * np.log(cdf_pk);
    
    if np.sign(skew_x)>0:
        max_est[i] = np.trapz(pdf_pk*X_max,y_pk);
        min_est[i] = np.trapz(pdf_pk*X_min,y_pk);
        max_std[i] = np.trapz((X_max-max_est[i])**2*pdf_pk,y_pk);
        min_std[i] = np.trapz((X_min-min_est[i])**2*pdf_pk,y_pk);
    else:
        max_est[i] = -np.trapz(pdf_pk*X_min,y_pk);
        min_est[i] = -np.trapz(pdf_pk*X_max,y_pk);
        max_std[i] = np.trapz((-X_min-max_est[i])**2*pdf_pk,y_pk);
        min_std[i] = np.trapz((-X_max-min_est[i])**2*pdf_pk,y_pk);