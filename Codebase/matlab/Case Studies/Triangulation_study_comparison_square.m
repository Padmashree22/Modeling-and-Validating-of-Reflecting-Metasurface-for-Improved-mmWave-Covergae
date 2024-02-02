% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, workspace, and close all figures
clc; 
clear all; 
close all;

% Define file paths for simulation data
file_name = 'path_to_file_with_square_plate.csv'; 
file_name1 ='path_to_file_with_square_plate.csv';
file_name2 ='path_to_file_with_square_plate.csv';

% Check if files exist, throw an error if not
assert(exist(file_name, 'file') == 2, 'File does not exist')
assert(exist(file_name1, 'file') == 2, 'File does not exist')
assert(exist(file_name2, 'file') == 2, 'File does not exist')

% Load simulation data for each model
[RCS_mat, phi_vec, theta_vec] = load_simulation_data(file_name);
[RCS_mat1, phi_vec1, theta_vec1] = load_simulation_data(file_name1);
[RCS_mat2, phi_vec2, theta_vec2] = load_simulation_data(file_name2);

%% Upper limit of RCS
% Compute the maximum RCS value for scaling the plots
RCS_max = max(RCS_mat, [], 'all'); % dB
RCS_max = ceil(RCS_max./10).*10;
RCS_min = RCS_max - 50; % dB

% Combine RCS matrices for triangles and calculate the total RCS
RCS = (sum(sqrt(RCS_mat1 + RCS_mat2), [3, 4])).^2;
Max_RCS = max(max(RCS));

%% Plotting
% Plot pattern slices for square and combined triangle models
figure;
theta_slice = 0;
k = find(theta_vec == theta_slice);
plot_pattern_slice(RCS_mat(:, k), theta_vec, ...
    'Gain Square[dB]', 'Elevation angle at receiver [deg]', ...
    RCS_min, RCS_max, ...
    -90, 90, 'r');
hold on;
plot_pattern_slice(RCS(:, k), theta_vec, ...
    'Gain Triangle Together[dB]', 'Elevation angle at receiver [deg]', ...
    RCS_min, RCS_max, ...
    -90, 90, 'b');
legend('Square', 'Triangle Together');

% Create a 2x1 subplot for heatmaps of square and combined triangle models
figure;
subplot(1, 2, 1);
plot_heatmap_simulation(RCS_mat, phi_vec, theta_vec, ...
    'Gain Square [dB]', 'Phi [deg]', 'Theta [deg]', ...
    RCS_min, RCS_max);
ylim([-20 20]);
xlim([-20 20]);

subplot(1, 2, 2);
plot_heatmap_simulation(RCS, phi_vec, theta_vec, ...
    'Gain Triangle Together [dB]', 'Phi [deg]', 'Theta [deg]', ...
    RCS_min, RCS_max);
ylim([-20 20]);
xlim([-20 20]);
