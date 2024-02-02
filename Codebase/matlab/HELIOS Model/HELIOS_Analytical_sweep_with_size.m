% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, workspace, and close all figures
clc; 
clear all; 
close all;

% Add paths to helper functions
addpath('Analytical_Model_functions\') 
addpath('Plotting_functions\')

%% Parameters

% Speed of light [m/s]
c = physconst('Light');

% Frequency [Hz]
f = 28e9; 

% Incident azimuth and elevation angles [deg]
phi_in = 0; 
theta_in = 0;

%% Pattern Resolution Parameters

% Resolution of RCS [deg]
pattern_resolution = 0.1; 

% Minimum and maximum reflected azimuth waves [deg]
phi_out_min = -90; 
phi_out_max = 90; 

% Minimum and maximum reflected elevation waves [deg]
theta_out_min = -90; 
theta_out_max = 90; 

%% Antenna Module Configuration

% Number of modules along Y and Z axes
N = 1; 
M = 1; 

% Dimensions of the antenna module [m]
a = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]; 

% Horizontal and vertical slope angles [deg]
alpha = 10 .* ones(length(M), length(N)); 
beta = 10 .* ones(length(M), length(N)); 

%% Reflected Wave Parameters

% Reflected elevation and azimuth waves [deg]
theta_out = theta_out_min:pattern_resolution:theta_out_max; 
phi_out = phi_out_min:pattern_resolution:phi_out_max; 

% Wavelength [m]
lambda = c/f; 

% Initialize sigma_sum_max vector
sigma_sum_max = zeros(size(a));

% Loop over antenna module dimensions
for i = 1:length(a)
    b = a(i);
    % Call the analytical model function to calculate RCS
    [sigma_sum_max(i), ~, ~, ~] = fn_SF_HELIOS_array_eff_modl(a(i), b, N, M, alpha, beta, phi_in, theta_in, phi_out, theta_out, lambda);
end

% Plot analytical RCS calculation results
p_RCS = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
plot(p_RCS, sigma_sum_max, 'r', 'LineWidth', 2);
hold on;

%% Simulation Results

% Specify folder path for simulation results
folderPaths = ["folder_path"];

% Loop over folder paths
n = 1;
while n <= 1
    folderPath = folderPaths(n);
    fileList = dir(fullfile(folderPath, '*.csv'));
    RCS_max = zeros(1, numel(fileList));

    % Iterate through each file in the list
    for i = 1:numel(fileList)
        fileName = fileList(i).name;     
        filePath = fullfile(folderPath, fileName);  
        
        % Load simulation data using a helper function (not provided here)
        [RCS_mat, ~, ~] = load_simulation_data(filePath); 
        RCS_max(i) = max(max(RCS_mat));
    end
    
    % Plot simulation results
    b = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
    plot(b, RCS_max, 'b--o', 'LineWidth', 2);
    n = n + 1;
end

% Labeling and legend
xlabel('Width b (cm) of reflector module');
ylabel('Gain (dB)')
legend('RCS Analytical Calculation', 'Simulation Results');
