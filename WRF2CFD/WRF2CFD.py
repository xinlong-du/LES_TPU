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
import pandas as pd
from scipy.io import netcdf_file
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

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
file2read = netcdf_file(wrf_file_path,'r');
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

# Step 7: Calculate Height (in meters)
g = 9.81; # Gravitational acceleration (m/s²)
z_full = (ph_data.data + phb_data.data) / g; # Height calculation
z_centered = 0.5 * (z_full[:, 0:-1, :, :] + z_full[:, 1:, :, :]); # Interpolated to grid center

# Process U and V data to align with the grid center
u_data_centered = 0.5 * (u_data[:, :, :, 0:-1] + u_data[:, :, :, 1:]);
v_data_centered = 0.5 * (v_data[:, :, 0:-1, :] + v_data[:, :, 1:, :]);

# Step 8: Read CFD Grid Data
cfd_table = pd.read_csv(cfd_grid_file); # Read CSV file using readtable, skipping the header
cfd_coords = cfd_table.to_numpy();  # Convert table data to matrix
face_id = cfd_coords[:, 0]; # Extract Face ID
cfd_x = cfd_coords[:, 1];   # Extract X coordinate (longitude)
cfd_y = cfd_coords[:, 2];   # Extract Y coordinate (latitude)
cfd_z = cfd_coords[:, 3];   # Extract Z coordinate (height)

# Output the size of the read data to confirm correct reading
print('Successfully read CFD grid data: %d points\n');

# Sort CFD grid points if necessary
sort_idx=np.argsort(cfd_x,kind='stable');
cfd_x = cfd_x[sort_idx];
cfd_y = cfd_y[sort_idx];
cfd_z = cfd_z[sort_idx];

# Step 9: Align WRF and CFD Boundaries Using Corner Points
cfd_corners_x = [min(cfd_x), min(cfd_x), max(cfd_x), max(cfd_x)];
cfd_corners_y = [min(cfd_y), max(cfd_y), min(cfd_y), max(cfd_y)];

if boundary_direction == 1:
    wrf_x = np.tile(cfd_corners_x[0], (num_points));  # X coordinate remains constant
    wrf_y = cfd_corners_y[0] + np.linspace(0,num_points-2,num_points-1) * grid_resolution;  # Generate Y coordinates
    wrf_y = np.append(wrf_y,boundary_length);
    
elif boundary_direction == 2:
    wrf_y = np.tile(cfd_corners_y[0], (num_points));  # Y coordinate remains constant
    wrf_x = cfd_corners_x[0] + np.linspace(0,num_points-2,num_points-1) * grid_resolution;  # Generate X coordinates
    wrf_x = [wrf_x, boundary_length];
    
else:
    print('Invalid boundary direction input. Enter 1 for east-west boundary or 2 for north-south boundary.');

z_full_point = z_full[0, :, 2, 0]; 
wrf_z = z_full_point - terrain_data[0, 2, 0];

# Step 10: Interpolation and Data Output
results_matrix = np.zeros([len(face_id), 6]);  # Matrix to store interpolated results
results_matrix[:, 0] = face_id;  # First column is face_id
results_matrix[:, 1] = cfd_x;    # Second column is x coordinate
results_matrix[:, 2] = cfd_y;    # Third column is y coordinate
results_matrix[:, 3] = cfd_z;    # Fourth column is z coordinate

for j in range(0,num_files):
    time_str = str(times[j,:], 'UTF-8');  # Extract and clean the time string
    time_str_safe = time_str.replace(':', '-');  # Replace colons with hyphens
    
    u_wrf=np.zeros(shape=(np.size(u_data_centered,1), num_points));
    v_wrf=np.zeros(shape=(np.size(v_data_centered,1), num_points));
    for i in range(0,num_points):

        # Get the WRF grid index for each CFD boundary point
        locXtemp = abs(lonGrid[0,0,:] - cfd_lon[i]);
        locYtemp = abs(latGrid[0,:,0] - cfd_lat[i]);
        locX = np.argmin(locXtemp);
        locY = np.argmin(locYtemp);

        # Ensure locX and locY are within valid ranges
        locX = max(0, min(locX, np.size(u_data_centered, 3)-1));
        locY = max(0, min(locY, np.size(v_data_centered, 2)-1));

        # Extract wind profile data for this point
        hgt = z_centered[j, :, locY, locX];  # Use centered height data
        ter_point = terrain_data[0, locY, locX];

        u_point = u_data_centered[j, :, locY, locX];
        v_point = v_data_centered[j, :, locY, locX];
        w_point = w_data[j, :, locY, locX];         # Height dimension
        p_point = pressure_data[j, :, locY, locX];  # Height dimension
        spd_point = np.sqrt(u_point**2 + v_point**2);

        # Output wind profile data to a file
        output_file = '%s/D03_%s_%d.txt' % (output_dir, time_str_safe, i);
        alist_1 = np.stack((hgt, u_point, v_point, spd_point, p_point),axis=1);
        np.savetxt(output_file, alist_1, fmt='%5.4f', delimiter='\t');

        # Output terrain height data to a file
        output_file_terrain = '%s/D03_terrain_height.txt' % (output_dir);
        alist_2 = np.vstack((cfd_lon[i], cfd_lat[i], ter_point)).transpose();
        f=open(output_file_terrain,'a');
        np.savetxt(f, alist_2, fmt='%5.4f', delimiter='\t');
        f.close();
        
        u_wrf[:,i] = u_point;
        v_wrf[:,i] = v_point;

    wrf_z = hgt - ter_point; 
    
    # Add wind speed data for z = 0
    z_0 = 0;  
    u_0 = np.zeros((1,np.size(u_wrf[1, :])));  
    v_0 = np.zeros((1,np.size(v_wrf[1, :])));  

    # Extend wrf_z, u_wrf, v_wrf to include data for z = 0
    wrf_z_extended = np.append(z_0, wrf_z);  
    u_wrf_extended = np.vstack((u_0, u_wrf));  
    v_wrf_extended = np.vstack((v_0, v_wrf));  

    # Get WRF grid dimensions
    num_levels = np.size(u_wrf_extended,0);
    num_horiz = np.size(u_wrf_extended,1);

    if boundary_direction == 1:
        wrf_y_2d, wrf_z_2d = np.meshgrid(wrf_y, wrf_z_extended, indexing='ij');  

        u_wrf_2d = u_wrf_extended.transpose();  
        v_wrf_2d = v_wrf_extended.transpose();

        wrf_y_flat = wrf_y_2d.flatten(order='F');
        wrf_z_flat = wrf_z_2d.flatten(order='F');
        u_wrf_flat = u_wrf_2d.flatten(order='F');
        v_wrf_flat = v_wrf_2d.flatten(order='F');

        points=np.stack((wrf_y_flat,wrf_z_flat),axis=1);
        for i in range(0,len(cfd_y)):
            y_cfd = cfd_y[i];
            z_cfd = cfd_z[i];

            u_cfd_interp = griddata(points,u_wrf_flat,(y_cfd, z_cfd),method='linear');
            v_cfd_interp = griddata(points,v_wrf_flat,(y_cfd, z_cfd),method='linear');

            results_matrix[i, 4] = u_cfd_interp;  
            results_matrix[i, 5] = v_cfd_interp;
    
    elif boundary_direction == 2:
        wrf_x_2d, wrf_z_2d = np.meshgrid(wrf_x, wrf_z_extended, indexing='ij');  

        u_wrf_2d = u_wrf_extended.transpose();  
        v_wrf_2d = v_wrf_extended.transpose();

        wrf_x_flat = wrf_x_2d.flatten(order='F');
        wrf_z_flat = wrf_z_2d.flatten(order='F');
        u_wrf_flat = u_wrf_2d.flatten(order='F');
        v_wrf_flat = v_wrf_2d.flatten(order='F');

        points=np.stack((wrf_x_flat,wrf_z_flat),axis=1);
        for i in range(0,len(cfd_x)):
            x_cfd = cfd_x(i);
            z_cfd = cfd_z(i);

            u_cfd_interp = griddata(points,u_wrf_flat,(x_cfd, z_cfd),method='cubic');
            v_cfd_interp = griddata(points,v_wrf_flat,(x_cfd, z_cfd),method='cubic');

            results_matrix[i, 4] = u_cfd_interp;  
            results_matrix[i, 5] = v_cfd_interp;  

    # Save Results to a File
    column_names = np.array(['faceID', 'x', 'y', 'z', 'xvelocity', 'yvelocity']);
    column_names = column_names.reshape(1,6);
    output_file = '%s/D03_velocity_%s.txt'% (output_dir, time_str_safe);
    f=open(output_file,'w');
    np.savetxt(f, column_names, fmt='%s', delimiter='\t');
    f.close();
    f=open(output_file,'a');
    np.savetxt(f, results_matrix, fmt='%5.4f', delimiter='\t');
    f.close();

    # Plot - Compare Interpolated Results and Original Results
    fig=plt.figure(figsize=(6,6));
    ax=fig.add_subplot(projection='3d');
    fig_font_size = 8;
    
    if boundary_direction == 1:
        ax.scatter(wrf_y_2d, u_wrf_2d, wrf_z_2d, s=10, label='Original WRF Data');
        ax.scatter(results_matrix[:, 2], results_matrix[:, 4], results_matrix[:, 3], s=5, label='Interpolated CFD Data');
        ax.set_xlabel('Y (m)');
        ax.set_ylabel('U Velocity (m/s)');
        ax.set_zlabel('Z (m)');
    elif boundary_direction == 2:
        ax.scatter(wrf_x_2d, u_wrf_2d, wrf_z_2d, s=10, label='Original WRF Data');
        ax.scatter(results_matrix[:, 1], results_matrix[:, 4], results_matrix[:, 3], s=5, label='Interpolated CFD Data');
        ax.set_xlabel('X (m)');
        ax.set_ylabel('U Velocity (m/s)');
        ax.set_zlabel('Z (m)');
    
    plt.title('Comparison of Original WRF and Interpolated CFD Results at Time %s' % time_str_safe);
    plt.legend(loc='upper right')

    # Set z-axis range
    z_max = max(cfd_z);
    ax.set_zlim([0, z_max]);    

    # Save the plot
    figName='%s/Comparison_%s.png' % (output_dir, time_str_safe);
    plt.savefig(figName, transparent=False, bbox_inches='tight', dpi=100)
    plt.show()