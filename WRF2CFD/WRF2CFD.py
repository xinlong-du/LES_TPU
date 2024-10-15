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
from scipy.io import netcdf_file
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

#file2read.close();

import netCDF4
file2read2 = netCDF4.Dataset(wrf_file_path,'r')
times2 = file2read2.variables['Times']  # access a variable in the file
lonGrid2 = file2read2.variables['XLONG'];
latGrid2 = file2read2.variables['XLAT'];

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
results_matrix[:, 1] = face_id;  # First column is face_id
results_matrix[:, 2] = cfd_x;    # Second column is x coordinate
results_matrix[:, 3] = cfd_y;    # Third column is y coordinate
results_matrix[:, 4] = cfd_z;    # Fourth column is z coordinate

for j in range(1,num_files):
    time_str = strtrim(times[:, j]);  # Extract and clean the time string
    time_str_safe = strrep(time_str, ':', '-');  # Replace colons with hyphens
    
    for i in range(1,num_points):

        # Get the WRF grid index for each CFD boundary point
        locX = min(abs(lonGrid[:,1,1] - cfd_lon[i]));
        locY = min(abs(latGrid[1,:,1] - cfd_lat[i]));

        # Ensure locX and locY are within valid ranges
        locX = max(1, min(locX, size(u_data_centered, 1)));
        locY = max(1, min(locY, size(v_data_centered, 2)));

        # Extract wind profile data for this point
        hgt = squeeze(z_centered[locX, locY, :, j]);  # Use centered height data
        ter_point = terrain_data(locX, locY);

        u_point = squeeze(u_data_centered[locX, locY, :, j]);
        v_point = squeeze(v_data_centered[locX, locY, :, j]);
        w_point = squeeze(w_data[locX, locY, :, j]);         # Height dimension
        p_point = squeeze(pressure_data[locX, locY, :, j]);  # Height dimension
        spd_point = sqrt(u_point**2 + v_point**2);

        # Output wind profile data to a file
        output_file = sprintf('%sD03_%s_%d_%d.txt', output_dir, time_str_safe, i);
        alist_1 = [hgt, u_point, v_point, spd_point, p_point];
        writematrix(alist_1, output_file, 'Delimiter', 'tab');

        # Output terrain height data to a file
        output_file_terrain = sprintf('%sD03_terrain_height.txt', output_dir);
        alist_2 = [cfd_lon(i), cfd_lat(i), ter_point];
        writematrix(alist_2, output_file_terrain, 'Delimiter', 'tab', 'WriteMode', 'append');
        
        u_wrf[:,i] = u_point;
        v_wrf[:,i] = v_point;

    end
    
    wrf_z = hgt - ter_point; 
    
    # Add wind speed data for z = 0
    z_0 = 0;  
    u_0 = zeros(size(u_wrf[1, :]));  
    v_0 = zeros(size(v_wrf[1, :]));  

    # Extend wrf_z, u_wrf, v_wrf to include data for z = 0
    wrf_z_extended = np.append(z_0, wrf_z);  
    u_wrf_extended = np.append(u_0, u_wrf);  
    v_wrf_extended = np.append(v_0, v_wrf);  

    # Get WRF grid dimensions
    [num_levels, num_horiz] = size(u_wrf_extended); 

    if boundary_direction == 1:
        [wrf_y_2d, wrf_z_2d] = ndgrid(wrf_y, wrf_z_extended);  

        u_wrf_2d = permute(u_wrf_extended, [2, 1]);  
        v_wrf_2d = permute(v_wrf_extended, [2, 1]);

        wrf_y_flat = double(wrf_y_2d);
        wrf_z_flat = double(wrf_z_2d);
        u_wrf_flat = double(u_wrf_2d);
        v_wrf_flat = double(v_wrf_2d);

        F_u = scatteredInterpolant(wrf_y_flat, wrf_z_flat, u_wrf_flat, 'linear', 'none');
        F_v = scatteredInterpolant(wrf_y_flat, wrf_z_flat, v_wrf_flat, 'linear', 'none');    
        
        for i in range(1,length(cfd_y)):
            y_cfd = cfd_y(i);
            z_cfd = cfd_z(i);

            u_cfd_interp = F_u(y_cfd, z_cfd);
            v_cfd_interp = F_v(y_cfd, z_cfd);

            results_matrix[i, 5] = u_cfd_interp;  
            results_matrix[i, 6] = v_cfd_interp;  
        end
    
    elif boundary_direction == 2:
        [wrf_x_2d, wrf_z_2d] = ndgrid(wrf_x, wrf_z_extended);  

        u_wrf_2d = permute(u_wrf_extended, [2, 1]);  
        v_wrf_2d = permute(v_wrf_extended, [2, 1]);

        wrf_x_flat = double(wrf_x_2d);
        wrf_z_flat = double(wrf_z_2d);
        u_wrf_flat = double(u_wrf_2d);
        v_wrf_flat = double(v_wrf_2d);

        F_u = scatteredInterpolant(wrf_x_flat, wrf_z_flat, u_wrf_flat, 'natural', 'none');
        F_v = scatteredInterpolant(wrf_x_flat, wrf_z_flat, v_wrf_flat, 'natural', 'none');    
        
        for i in range(1,length(cfd_x)):
            x_cfd = cfd_x(i);
            z_cfd = cfd_z(i);

            u_cfd_interp = F_u(x_cfd, z_cfd);
            v_cfd_interp = F_v(x_cfd, z_cfd);

            results_matrix[i, 5] = u_cfd_interp;  
            results_matrix[i, 6] = v_cfd_interp;  
        end
        
    end

    # Save Results to a File
    column_names = {'faceID', 'x', 'y', 'z', 'xvelocity', 'yvelocity'};
    output_file = sprintf('%sD03_velocity_%s.txt', output_dir, time_str_safe);
    results_table = array2table(results_matrix, 'VariableNames', column_names);
    writetable(results_table, output_file, 'Delimiter', '\t');

    # Plot - Compare Interpolated Results and Original Results
    figure;
    
    if boundary_direction == 1:
        scatter3(wrf_y_2d, u_wrf_2d, wrf_z_2d, 50, 'k', 'filled');
        #hold on;
        scatter3(results_matrix[:, 3], results_matrix[:, 5], results_matrix[:, 4], 'r');  
        colorbar;
        xlabel('Y (m)');
        ylabel('U Velocity (m/s)');
        zlabel('Z (m)');
    elif boundary_direction == 2:
        scatter3(wrf_x_2d, u_wrf_2d, wrf_z_2d, 50, 'k', 'filled');
        #hold on;
        scatter3(results_matrix[:, 2], results_matrix[:, 5], results_matrix[:, 4], 'r');  
        colorbar;
        xlabel('X (m)');
        ylabel('U Velocity (m/s)');
        zlabel('Z (m)');
    end
    

    title(sprintf('Comparison of Original WRF and Interpolated CFD Results at Time %s', time_str_safe));
    legend('Original WRF Data', 'Interpolated CFD Data');

    # Set z-axis range
    z_max = max(cfd_z);
    zlim([0, z_max]);    

    # Save the plot
    saveas(gcf, sprintf('%sComparison_%s.png', output_dir, time_str_safe));
    
end