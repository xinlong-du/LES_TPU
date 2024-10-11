# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 13:49:43 2024

@author: xinlo
"""

#wind field interpolation from WRF to CFD grid
# Step 1:
from pathlib import Path
import math
import numpy as np


# Step 2: Parameter Settings
wrf_file_path = './input/wrfout_d03_2008-08-06_020000.nc'; # WRF output file path
output_dir = './output';    # Output directory
cfd_grid_file = './input/test_CFD_grid_10m_spacing.csv'; # CFD grid coordinates file

# Ensure output directory exists
Path("./output").mkdir(parents=True, exist_ok=True)

# Step 3: Input CFD Boundary Parameters
start_lat = 22.29;   #input('Enter starting latitude of CFD boundary:');
start_lon = 114.1977; #input('Enter starting longitude of CFD boundary:');
boundary_length = 2000; #input('Enter boundary length (meters):');
grid_resolution = 444.4;   #input('Enter WRF grid resolution (in meters):');
boundary_direction = 1;#input('Enter boundary direction (1 for east-west, 2 for north-south):');

R = 6371000; # Earth's average radius (meters)

# Step 4: Generate CFD Boundary Points
if boundary_direction == 1:
    # East-west boundary - calculate latitude
    delta_lat_per_grid = grid_resolution / R * (180/math.pi);
    num_points = round(boundary_length / grid_resolution) + 1;
    cfd_lat = start_lat + np.linspace(0,(num_points-1),num_points) * delta_lat_per_grid;
    cfd_lon = np.tile(start_lon, (num_points));  # Longitude remains constant
    
elif boundary_direction == 2:
    # North-south boundary - calculate longitude
    delta_lon_per_grid = grid_resolution / (R * math.cos(math.radians(start_lat))) * (180/math.pi);
    num_points = round(boundary_length / grid_resolution) + 1;
    cfd_lon = start_lon + np.linspace(0,(num_points-1),num_points) * delta_lon_per_grid;
    cfd_lat = np.tile(start_lat, (num_points));  # Latitude remains constant
    
else:
    print('Invalid boundary direction input. Enter 1 for east-west or 2 for north-south.');

# Step 5: Output CFD Boundary Information
print('Number of generated CFD boundary points:', num_points);
print('Starting point: Latitude, Longitude', start_lat, start_lon);
print('Ending point: Latitude, Longitude', cfd_lat[-1], cfd_lon[-1]);

# Step 6: Read WRF data
from scipy.io import netcdf
file2read = netcdf.NetCDFFile(wrf_file_path,'r');
times = file2read.variables['Times'];
num_files = times.shape[0];

lonGrid = file2read.variables['XLONG'];
latGrid = file2read.variables['XLAT'];
ph_data = file2read.variables['PH'];   # Perturbation geopotential
phb_data = file2read.variables['PHB']; # Base-state geopotential
u_data = file2read.variables['U'];
v_data = file2read.variables['V'];
w_data = file2read.variables['W'];
pressure_data = file2read.variables['P'];
terrain_data = file2read.variables['HGT'];

import netCDF4
file2read2 = netCDF4.Dataset(wrf_file_path,'r')
times2 = file2read2.variables['Times']  # access a variable in the file
lonGrid2 = file2read2.variables['XLONG'];
latGrid2 = file2read2.variables['XLAT'];
