% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, workspace, and close all figures
clc; 
clear all; 
close all;

% Add path to the folder containing IRS model functions
addpath('IRS_model_functions\') % Load helper_function

%% Parameters

% Constants
c = 3e8;          % Speed of light in m/s
f = 28e9;         % Frequency in Hz
lambda = c / f;   % Wavelength in meters
gamma = 0.9;      % Reflection coefficient

% Angles for receiver (rx) and transmitter (tx)
rx_phi = -90:0.1:90;
rx_theta = 0.000001;
tx_phi = 0.00001;
tx_theta = 0.0001;
tx_gamma = 0.000001;

% Compute the reflection coefficient for the receiver
[rx_gamma] = gamma_rx(rx_theta, rx_phi);

%% Size and dimensions

% Number of elements in the IRS array
N = 1;
M = 1;

% Dimensions of each IRS element
b = 0.366;
a = b / N;

%% IRS Models

% Compute gains for different IRS models
[IRS_gain1] = IRS_model1(a, b, N, M, tx_phi, tx_theta, rx_phi, rx_theta, lambda); % IRS MODEL 1 (𝛾_Rx=𝛾_dev)
[IRS_gain2] = IRS_model2(a, b, N, M, tx_phi, tx_theta, rx_phi, rx_theta, lambda); % IRS MODEL 2
[IRS_gain3] = IRS_model3(a, b, N, M, tx_phi, tx_theta, rx_phi, rx_theta, lambda); % IRS MODEL 3

% Plot the results
plot(rx_phi, IRS_gain1, rx_phi, IRS_gain2, rx_phi, IRS_gain3, 'LineWidth', 2)
xlabel('Azimuth angle of receiver [deg]');
ylabel('Gain at IRS [dB]');
legend('Model 1-Özdogan', 'Model 2-Ntontin', 'Model 3-Tang');
