% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, close all figures, and clear workspace
clc;
close all;
clear all;

% Add paths to helper functions
addpath('Analytical_Model_functions\')
addpath('Plotting_functions\')

%% Parameters

% Physical constants
c = physconst('Light'); % Speed of light [m/s]
f = 28e9; % Frequency [Hz]
phi_in = 0; % Incident azimuth angle [deg]
theta_in = 0; % Incident elevation angle [deg]

%% Resolution and Range for RCS

% Resolution for RCS calculation
pattern_resolution = 0.2; % Resolution of RCS [deg]

% Ranges for reflected azimuth and elevation waves
phi_out_min = -90; % Minimum reflected azimuth wave [deg]
phi_out_max = 90; % Maximum reflected azimuth wave [deg]
theta_out_min = -90; % Minimum reflected elevation wave [deg]
theta_out_max = 90; % Maximum reflected elevation wave [deg]

%% Module Configuration and Dimensions

% Number of modules along Y and Z axes
N = 1:2; % Modules along Y axis 
M = 1:2; % Modules along Z axis

% Dimensions of modules
a = 0.1; % Dimensions along Y axis [m]
b = 0.1; % Dimensions along Z axis [m]

% Slope angles
alpha = (10) .* ones(length(M), length(N)); % Horizontal Slope angle [deg]
beta = (10) .* ones(length(M), length(N)); % Vertical Slope angle [deg]

%% Wavelength Calculation

lambda = c / f; % Wavelength [m]

%% Analytical RCS Calculation

% Calculate effective RCS using analytical model
[sigma_sum_max, sigma_sum_mat, eta_a, eta_b] = fn_SF_HELIOS_array_eff(a, b, N, M, alpha, beta, phi_in, theta_in, phi_out, theta_out, lambda);
sigma_sum_mat = (10 * log10(sigma_sum_mat(:,:,end).^2));

%% Simulations

% Load simulation data from a CSV file
file_name = ['C:\Users\padma\OneDrive\Desktop\Matlab_codes\Results\Results_2x2\Alphabeta_10.csv'];
assert(exist(file_name, 'file') == 2, 'File does not exist')
[RCS_mat, phi_vec, theta_vec] = load_simulation_data(file_name);

% Find
