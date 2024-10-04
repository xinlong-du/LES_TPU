close all; clear; clc;

Cp=readmatrix('../postProcessing/LES_data/Cp/lesCpF0.csv');
Cp1=Cp';

[max_est min_est max_std min_std]=maxminest(Cp1,1);