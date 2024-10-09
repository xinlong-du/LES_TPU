% -------------------------------------------------------------------------
% MATLAB Script for Wind Field Interpolation from WRF to CFD Grid
% -------------------------------------------------------------------------
% Version: 1.0
% Developed by: Zhejiang University in 2024
% -------------------------------------------------------------------------
% Description:
% This script performs interpolation of wind field data from WRF model 
% output to a specified CFD grid. It supports handling east-west and 
% north-south boundary setups, ensuring proper alignment of the coordinate 
% system. The results include wind profile data, terrain height, and 
% comparison plots between the original WRF data and the interpolated CFD 
% data.
%
% Note: This script is for educational and research purposes only and is 
% not intended for commercial use.
%
% Citation: Huang MF, Liao Sunce, Lou Wenjuan, Lin Wei, and Kareem Ahsan. 2024. 
% Multi-scale simulation of typhoon wind field at building scale utilizing mesoscale model with nested large eddy simulation. 
% Journal of Wind Engineering and Industrial Aerodynamics, 249, 105733. 
% -------------------------------------------------------------------------

% Step 1: Initialization
clc; clear; close all;

% Step 2: Parameter Settings
wrf_file_path = ''; % WRF output file path
output_dir = ''; % Output directory
cfd_grid_file = ''; % CFD grid coordinates file

% Ensure output directory exists
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Step 3: Input CFD Boundary Parameters
start_lat = input('Enter starting latitude of CFD boundary:');
start_lon = input('Enter starting longitude of CFD boundary:');
boundary_length = input('Enter boundary length (meters):');
grid_resolution = input('Enter WRF grid resolution (in meters):');
boundary_direction = input('Enter boundary direction (1 for east-west, 2 for north-south):');

R = 6371000; % Earth's average radius (meters)

% Step 4: Generate CFD Boundary Points
if boundary_direction == 1
    % East-west boundary - calculate latitude
    delta_lat_per_grid = grid_resolution / R * (180/pi);
    num_points = round(boundary_length / grid_resolution) + 1;
    cfd_lat = start_lat + (0:(num_points-1)) * delta_lat_per_grid;
    cfd_lon = repmat(start_lon, 1, num_points);  % Longitude remains constant
    
elseif boundary_direction == 2
    % North-south boundary - calculate longitude
    delta_lon_per_grid = grid_resolution / (R * cos(deg2rad(start_lat))) * (180/pi);
    num_points = round(boundary_length / grid_resolution) + 1;
    cfd_lon = start_lon + (0:(num_points-1)) * delta_lon_per_grid;
    cfd_lat = repmat(start_lat, 1, num_points);  % Latitude remains constant
    
else
    error('Invalid boundary direction input. Enter 1 for east-west or 2 for north-south.');
end

% Step 5: Output CFD Boundary Information
fprintf('Number of generated CFD boundary points: %d\n', num_points);
fprintf('Starting point: Latitude %.6f, Longitude %.6f\n', start_lat, start_lon);
fprintf('Ending point: Latitude %.6f, Longitude %.6f\n', cfd_lat(end), cfd_lon(end));

% Step 6: Read WRF Data
times = ncread(wrf_file_path, 'Times');
num_files = size(times, 2);

lonGrid = ncread(wrf_file_path, 'XLONG');
latGrid = ncread(wrf_file_path, 'XLAT');
ph_data = ncread(wrf_file_path, 'PH');   % Perturbation geopotential
phb_data = ncread(wrf_file_path, 'PHB'); % Base-state geopotential
u_data = ncread(wrf_file_path, 'U');
v_data = ncread(wrf_file_path, 'V');
w_data = ncread(wrf_file_path, 'W');
pressure_data = ncread(wrf_file_path, 'P');
terrain_data = ncread(wrf_file_path, 'HGT');

% Step 7: Calculate Height (in meters)
g = 9.81; % Gravitational acceleration (m/s²)
z_full = (ph_data + phb_data) / g; % Height calculation
z_centered = 0.5 * (z_full(:, :, 1:end-1, :) + z_full(:, :, 2:end, :)); % Interpolated to grid center

% Process U and V data to align with the grid center
u_data_centered = 0.5 * (u_data(1:end-1, :, :, :) + u_data(2:end, :, :, :));
v_data_centered = 0.5 * (v_data(:, 1:end-1, :, :) + v_data(:, 2:end, :, :));

% Step 8: Read CFD Grid Data
cfd_table = readtable(cfd_grid_file); % Read CSV file using readtable, skipping the header
cfd_coords = table2array(cfd_table); % Convert table data to matrix
face_id = cfd_coords(:, 1); % Extract Face ID
cfd_x = cfd_coords(:, 2);   % Extract X coordinate (longitude)
cfd_y = cfd_coords(:, 3);   % Extract Y coordinate (latitude)
cfd_z = cfd_coords(:, 4);   % Extract Z coordinate (height)

% Output the size of the read data to confirm correct reading
fprintf('Successfully read CFD grid data: %d points\n', size(cfd_coords, 1));

% Sort CFD grid points if necessary
[cfd_x, sort_idx] = sort(cfd_x);
cfd_y = cfd_y(sort_idx);
cfd_z = cfd_z(sort_idx);

% Step 9: Align WRF and CFD Boundaries Using Corner Points
cfd_corners_x = [min(cfd_x), min(cfd_x), max(cfd_x), max(cfd_x)];
cfd_corners_y = [min(cfd_y), max(cfd_y), min(cfd_y), max(cfd_y)];

if boundary_direction == 1
    wrf_x = repmat(cfd_corners_x(1), 1, num_points);  % X coordinate remains constant
    wrf_y = cfd_corners_y(1) + (0:(num_points-2)) * grid_resolution;  % Generate Y coordinates
    wrf_y = [wrf_y, boundary_length];
    
elseif boundary_direction == 2
    wrf_y = repmat(cfd_corners_y(1), 1, num_points);  % Y coordinate remains constant
    wrf_x = cfd_corners_x(1) + (0:(num_points-2)) * grid_resolution;  % Generate X coordinates
    wrf_x = [wrf_x, boundary_length];
    
else
    error('Invalid boundary direction input. Enter 1 for east-west boundary or 2 for north-south boundary.');
end

z_full_point = z_full(1, 3, :, 1); 
wrf_z = z_full_point(:) - terrain_data(1, 3, 1);

% Step 10: Interpolation and Data Output
results_matrix = zeros(length(face_id), 6);  % Matrix to store interpolated results
results_matrix(:, 1) = face_id;  % First column is face_id
results_matrix(:, 2) = cfd_x;    % Second column is x coordinate
results_matrix(:, 3) = cfd_y;    % Third column is y coordinate
results_matrix(:, 4) = cfd_z;    % Fourth column is z coordinate

for j = 1:num_files
    time_str = strtrim(times(:, j)');  % Extract and clean the time string
    time_str_safe = strrep(time_str, ':', '-');  % Replace colons with hyphens
    
    for i = 1:num_points

        % Get the WRF grid index for each CFD boundary point
        [~, locX] = min(abs(lonGrid(:,1,1) - cfd_lon(i)));
        [~, locY] = min(abs(latGrid(1,:,1) - cfd_lat(i)));

        % Ensure locX and locY are within valid ranges
        locX = max(1, min(locX, size(u_data_centered, 1)));
        locY = max(1, min(locY, size(v_data_centered, 2)));

        % Extract wind profile data for this point
        hgt = squeeze(z_centered(locX, locY, :, j));  % Use centered height data
        ter_point = terrain_data(locX, locY);

        u_point = squeeze(u_data_centered(locX, locY, :, j));
        v_point = squeeze(v_data_centered(locX, locY, :, j));
        w_point = squeeze(w_data(locX, locY, :, j));  % Height dimension
        p_point = squeeze(pressure_data(locX, locY, :, j));  % Height dimension
        spd_point = sqrt(u_point.^2 + v_point.^2);

        % Output wind profile data to a file
        output_file = sprintf('%sD03_%s_%d_%d.txt', output_dir, time_str_safe, i);
        alist_1 = [hgt, u_point, v_point, spd_point, p_point];
        writematrix(alist_1, output_file, 'Delimiter', 'tab');

        % Output terrain height data to a file
        output_file_terrain = sprintf('%sD03_terrain_height.txt', output_dir);
        alist_2 = [cfd_lon(i), cfd_lat(i), ter_point];
        writematrix(alist_2, output_file_terrain, 'Delimiter', 'tab', 'WriteMode', 'append');
        
        u_wrf(:,i) = u_point;
        v_wrf(:,i) = v_point;

    end
    
    wrf_z = hgt - ter_point; 
    
    % Add wind speed data for z = 0
    z_0 = 0;  
    u_0 = zeros(size(u_wrf(1, :)));  
    v_0 = zeros(size(v_wrf(1, :)));  

    % Extend wrf_z, u_wrf, v_wrf to include data for z = 0
    wrf_z_extended = [z_0; wrf_z];  
    u_wrf_extended = [u_0; u_wrf];  
    v_wrf_extended = [v_0; v_wrf];  

    % Get WRF grid dimensions
    [num_levels, num_horiz] = size(u_wrf_extended); 

    if boundary_direction == 1
        [wrf_y_2d, wrf_z_2d] = ndgrid(wrf_y, wrf_z_extended);  

        u_wrf_2d = permute(u_wrf_extended, [2, 1]);  
        v_wrf_2d = permute(v_wrf_extended, [2, 1]);

        wrf_y_flat = double(wrf_y_2d(:));
        wrf_z_flat = double(wrf_z_2d(:));
        u_wrf_flat = double(u_wrf_2d(:));
        v_wrf_flat = double(v_wrf_2d(:));

        F_u = scatteredInterpolant(wrf_y_flat, wrf_z_flat, u_wrf_flat, 'linear', 'none');
        F_v = scatteredInterpolant(wrf_y_flat, wrf_z_flat, v_wrf_flat, 'linear', 'none');    
        
        for i = 1:length(cfd_y)
            y_cfd = cfd_y(i);
            z_cfd = cfd_z(i);

            u_cfd_interp = F_u(y_cfd, z_cfd);
            v_cfd_interp = F_v(y_cfd, z_cfd);

            results_matrix(i, 5) = u_cfd_interp;  
            results_matrix(i, 6) = v_cfd_interp;  
        end
    
    elseif boundary_direction == 2
        [wrf_x_2d, wrf_z_2d] = ndgrid(wrf_x, wrf_z_extended);  

        u_wrf_2d = permute(u_wrf_extended, [2, 1]);  
        v_wrf_2d = permute(v_wrf_extended, [2, 1]);

        wrf_x_flat = double(wrf_x_2d(:));
        wrf_z_flat = double(wrf_z_2d(:));
        u_wrf_flat = double(u_wrf_2d(:));
        v_wrf_flat = double(v_wrf_2d(:));

        F_u = scatteredInterpolant(wrf_x_flat, wrf_z_flat, u_wrf_flat, 'natural', 'none');
        F_v = scatteredInterpolant(wrf_x_flat, wrf_z_flat, v_wrf_flat, 'natural', 'none');    
        
        for i = 1:length(cfd_x)
            x_cfd = cfd_x(i);
            z_cfd = cfd_z(i);

            u_cfd_interp = F_u(x_cfd, z_cfd);
            v_cfd_interp = F_v(x_cfd, z_cfd);

            results_matrix(i, 5) = u_cfd_interp;  
            results_matrix(i, 6) = v_cfd_interp;  
        end
        
    end

    % Save Results to a File
    column_names = {'faceID', 'x', 'y', 'z', 'xvelocity', 'yvelocity'};
    output_file = sprintf('%sD03_velocity_%s.txt', output_dir, time_str_safe);
    results_table = array2table(results_matrix, 'VariableNames', column_names);
    writetable(results_table, output_file, 'Delimiter', '\t');

    % Plot - Compare Interpolated Results and Original Results
    figure;
    
    if boundary_direction == 1
        scatter3(wrf_y_2d(:), u_wrf_2d(:), wrf_z_2d(:), 50, 'k', 'filled');
        hold on;
        scatter3(results_matrix(:, 3), results_matrix(:, 5), results_matrix(:, 4), 'r');  
        colorbar;
        xlabel('Y (m)');
        ylabel('U Velocity (m/s)');
        zlabel('Z (m)');
    elseif boundary_direction == 2
        scatter3(wrf_x_2d(:), u_wrf_2d(:), wrf_z_2d(:), 50, 'k', 'filled');
        hold on;
        scatter3(results_matrix(:, 2), results_matrix(:, 5), results_matrix(:, 4), 'r');  
        colorbar;
        xlabel('X (m)');
        ylabel('U Velocity (m/s)');
        zlabel('Z (m)');
    end
    

    title(sprintf('Comparison of Original WRF and Interpolated CFD Results at Time %s', time_str_safe));
    legend('Original WRF Data', 'Interpolated CFD Data');

    % Set z-axis range
    z_max = max(cfd_z);
    zlim([0 z_max]);    

    % Save the plot
    saveas(gcf, sprintf('%sComparison_%s.png', output_dir, time_str_safe));
    
end
