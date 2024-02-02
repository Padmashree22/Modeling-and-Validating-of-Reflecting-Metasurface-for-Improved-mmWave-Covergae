% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 


% Clear command window, close all figures, and clear workspace
clc; 
close all; 
clear all;

% Add path to helper function
addpath('Helper_function') % load helper_function

% Load LOS (Line-of-Sight) scenario grid
load('LOS_scenario1_grid','LOS_grid')

%% Define simulation parameters

% Frequency in Hz (28 GHz)
freq = 28e9;  

% Speed of light
c = physconst('Light');  

% Incident azimuth angle [deg]
phi_in = 0; 

% Incident elevation angle [deg]
theta_in = 0;

%% Define pattern resolution and output angles

% Resolution of RCS [deg]
pattern_resolution = 1; 

% Minimum and maximum reflected azimuth wave [deg]
phi_out_min = -90; 
phi_out_max = 90; 

% Minimum and maximum reflected elevation wave [deg]
theta_out_min = -90; 
theta_out_max = 90; 

%% Define array parameters

% Modules along Y axis
N = 1:32; 

% Modules along Z axis
M = 1; 

% Horizontal slope angle [deg]
alpha = (0.2) * ones(length(M), length(N)); 
alpha = flip(alpha);

% Vertical slope angle [deg]
beta = (0.1) * ones(length(M), length(N)); 

%% Define reflected angles

% Reflected elevation wave [deg]
theta_out = theta_out_min:pattern_resolution:theta_out_max; 

% Reflected azimuth wave [deg]
phi_out = phi_out_min:pattern_resolution:phi_out_max; 

% Wavelength [m]
lambda = c / freq; 

% Dimensions along Y axis [m]
a_dim = 3 / 32; 

% Dimensions along Z axis [m]
b_dim = 3; 

% Coefficients for varying reflector dimensions
k = [0.1:0.1:2.2, 2.5]; 

%% Gain_RCS

% Loop over different coefficients to vary reflector dimensions
for i = 1:length(k)
    a_dim(i) = k(i) / 32; 
    b_dim(i) = k(i); 
    
    % Call the function to calculate RCS and array efficiency
    [sigma_sum_max(i), sigma_sum_mat, eta_a, eta_b] = fn_SF_HELIOS_array_eff(a_dim(i), b_dim(i), N, M, alpha, beta, phi_in, theta_in, phi_out, theta_out, lambda);
    
    % Square of the last slice of sigma_sum_mat
    sigma_sum_mat = (sigma_sum_mat(:,:,end).^2);
    
    % Increment counter
    i = i + 1;
end

% Plot analytical RCS calculation results
plot(k, sigma_sum_max, 'r', 'LineWidth', 2);
hold on;

%% Max gain for A=B in urban scenario

% Folder containing simulation data files
folderPaths = 'folder_path';

% Initialize loop counter
n = 1;

% Loop over simulations
while n <= 1
    % Read CSV files and extract RCS data
    fileList = dir(fullfile(folderPaths, '*.csv'));
    frames = cell(1, numel(fileList));
    
    for i = 1:numel(fileList)
        fileName = fileList(i).name;
        filePath = fullfile(folderPaths, fileName);
        assert(exist(filePath, 'file') == 2, 'File does not exist')
        [RCS_mat, phi_vec, theta_vec] = load_simulation_data(filePath); % -> code appended
        RCS_Max(i) = max(max(RCS_mat(:)));
        Mean_Sim(i) = mean(RCS_mat(:));
    end

    % Plot simulation results
    plot(k, RCS_Max, 'b--o', 'LineWidth', 2);
    
    % Increment loop counter
    n = n + 1;
end

% Label the axes and add a legend
xlabel('Reflector Dimensions A=B [cm]');
ylabel('Gain [dB]');
legend('RCS Analytical Calculation', 'Simulation Results');
