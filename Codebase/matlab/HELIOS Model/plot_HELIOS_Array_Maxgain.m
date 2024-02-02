% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, close all figures, and clear all workspace variables
clc; 
close all; 
clear all;

% Add paths to helper functions
addpath('Analytical_Model_functions\')
addpath('Plotting_functions\')

%% Parameters

% Constants
c = physconst('Light'); % Speed of light [m/s]
f = 28e9; % Frequency [Hz]
phi_in = 0; % Incident azimuth angle [deg]
theta_in = 0; % Incident elevation angle [deg]

%% Resolution and Range

% Resolution of RCS [deg]
pattern_resolution = 0.1;

% Reflected azimuth and elevation wave ranges [deg]
phi_out_min = -90; 
phi_out_max = 90; 
theta_out_min = -90; 
theta_out_max = 90;

%% Antenna Array Configuration

% Number of modules along Y axis and Z axis
N = 1:8; 
M = 1; 

% Dimensions along Y and Z axes [m]
a = 0.1; 
b = 0.1; 

% Horizontal and vertical slope angles [deg]
alpha = 10 .* ones(length(M),length(N)); 
beta = 10 .* ones(length(M),length(N));

%% Reflected Wave Range and Wavelength

% Reflected elevation and azimuth wave ranges [deg]
theta_out = theta_out_min:pattern_resolution:theta_out_max; 
phi_out = phi_out_min:pattern_resolution:phi_out_max; 

% Wavelength [m]
lambda = c / f; 

%% Upper limit of RCS

% Maximum RCS [dB]
RCS_max = 10 .* log10(4 .* pi .* (a .* b ./ lambda).^2); 
RCS_max = ceil(RCS_max ./ 10) .* 10;
RCS_min = RCS_max - 50; % Minimum RCS [dB]

% Calculate RCS using analytical model function
[sigma_sum_max, sigma_sum_mat, eta_a, eta_b] = fn_SF_HELIOS_array_eff_modl(a, b, N, M, alpha, beta, phi_in, theta_in, phi_out, theta_out, lambda);

%% HFSS Simulations for max values

% Folder paths for HFSS simulation results
folderPaths = ["folderPaths"];

% Loop through folders and extract RCS values
n = 1;
while n <= 1
    folderPath = folderPaths(n);
    fileList = dir(fullfile(folderPath, '*.csv'));
    fileName = ["file_name.csv"];
    i = 1;
    while i <= length(fileName)
        % Iterate through each file in the list
        filePath = fullfile(folderPath, fileName(i));
        [RCS_mat, phi_vec, theta_vec] = load_simulation_data(filePath); % Function to load simulation data (not provided)
        RCS_max(n) = max(max(RCS_mat));
        i = i + 1;
    end
    n = n + 1;
end

%% Plotting

% Diagonal of the sigma_sum_max matrix
sigma_diag = diag(sigma_sum_max);

% X-axis values for plotting
p_HFSS = [1,2,3,4,5,6,7,8];

% Plot RCS and Simulation results
plot(p_HFSS, sigma_sum_max, 'r', p_HFSS, RCS_max, 'b--o', 'LineWidth', 2);
hold on;

% X-axis labels
N1 = {'1x1','2x2', '3x3', '4x4', '5x5', '6x6','7x7','8x8'};
xticks(p_HFSS);
xticklabels(N1);

% Y-axis limits and labels
ylim([5,50]);
xlabel('Number of Reflector Array');
ylabel('Gain (dB)');

% Legend
legend('RCS with shadowing factor', 'Simulation results');
