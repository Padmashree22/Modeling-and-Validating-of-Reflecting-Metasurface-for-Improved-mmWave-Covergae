% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, close all figures, and clear workspace
clc;
close all;
clear all;

%% Add necessary paths and load data
addpath('Analytical_Model_functions')  % Load helper functions
addpath('Plotting functions\')
load('LOS_scenario1_grid','LOS_grid')

%% Specify folder and file information
folder = 'folder_path'; % Replace with the folder path

file1 = 'Urban_narrow_K010.csv'; % Replace with the actual file name
file2 = 'Urban_narrow_K030.csv'; % Replace with the actual file name
file3 = 'Urban_narrow_K050.csv'; % Replace with the actual file name
file4 = 'Urban_narrow_K070.csv'; % Replace with the actual file name

%% Import data from CSV files
data1 = readtable(fullfile(folder, file1)); % Import data from the first CSV file
data2 = readtable(fullfile(folder, file2)); % Import data from the second CSV file
data3 = readtable(fullfile(folder, file3)); % Import data from the first CSV file
data4 = readtable(fullfile(folder, file4)); % Import data from the second CSV file

% Extract specific columns from the data
column1_data = table2array(data1); % Extract the 5th column 
column2_data = table2array(data2); % Extract the 2nd column
column3_data = table2array(data3); % Extract the 5th column 
column4_data = table2array(data4); % Extract the 2nd column

%% Define simulation parameters
freq = 28e9;  % Frequency in Hz (28 GHz)
c = physconst('Light');  % Speed of light
phi_in = 0; % Incident azimuth angle [deg]
theta_in = 0; % Incident elevation angle [deg]

%% Define pattern resolution and output angles
pattern_resolution = 0.5; % Resolution of RCS [deg]
phi_out_min = -90; % Minimum reflected azimuth wave [deg] 
phi_out_max = 90; % Maximum reflected azimuth wave [deg]
theta_out_min = -90; % Minimum reflected elevation wave [deg]
theta_out_max = 90; % Maximum reflected elevation wave [deg]

% Create vectors for reflected angles
theta_out = theta_out_min:pattern_resolution:theta_out_max; % Reflected elevation wave [deg]
phi_out = phi_out_min:pattern_resolution:phi_out_max; % Reflected azimuth wave [deg]
lambda = c / freq; % Wavelength [m]

%% Loop through different configurations and calculate RCS
k = 0.3;
N = 1:32; % Modules along Y axis 
M = 1; % Modules along Z axis
c1 = [1/2, length(M)/2, length(M)*length(N)/2, (length(M)*sqrt(length(N)))/2]; %beamwidth correction factors
d1 = [1/2, length(N)/2, length(M)*length(N)/2, (length(N)*sqrt(length(M)))/2];
i = 1;
for i = 1:length(c1)
    c = c1(i);
    d = d1(i);
    a_dim = k * (3/32); % Dimensions along Y axis [m]
    b_dim = k * (3); % Dimensions along Z axis [m]
    alpha = (27.6) * ones(length(M), length(N)); % Horizontal slope angle [deg]
    beta = (0) * ones(length(M), length(N)); % Vertical slope angle [deg]
    
    %% Calculate Gain_RCS
    [sigma_sum_max(i), sigma_sum_mat, eta_a, eta_b] = fn_SF_HELIOS_array_eff(a_dim, b_dim, N, M, alpha, beta, phi_in, theta_in, phi_out, theta_out, lambda, c, d);
    sigma_sum_mat = (10 * log10(sigma_sum_mat(:,:,end).^2));
    sigma_sum_mat1(:,i) = (sigma_sum_mat(:));
    Mean_RCS = mean(sigma_sum_mat(:));
    i = i + 1;
end

%% Violin plot to compare simulation data with calculated RCS
sim_data = column4_data(:,end);
Data_all = [(sim_data - sigma_sum_mat1(:,1)), (sim_data - sigma_sum_mat1(:,2)), (sim_data - sigma_sum_mat1(:,3)), (sim_data - sigma_sum_mat1(:,4))];
violin(Data_all, 'xlabel', {'(1/2,1/2)', '(M/2,N/2)', '(MN/2,MN/2)', '(M*sqrt(N))/2,(N*sqrt(M)/2'}, 'edgecolor', 'k');
title('K=0.7');
ylabel('Delta Gain [dB]');
ylim([-90 130]);
