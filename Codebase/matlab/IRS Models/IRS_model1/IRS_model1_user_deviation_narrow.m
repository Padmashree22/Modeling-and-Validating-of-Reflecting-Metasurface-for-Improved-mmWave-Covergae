% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear the command window, workspace, and close all figures
clc;
clear all;
close all;

% Add the path to the helper functions folder
addpath('IRS_model_functions\') % Load helper_function

%% Parameters
c = 3e8;      % Speed of light in m/s
f = 28e9;     % Frequency in Hz
lambda = c/f; % Wavelength in meters
gamma = 0.9;  % Reflection coefficient

phi_dev = 0:0.1:90;          % Azimuth deviation in degrees
theta_dev = [0.00000001,30,45,60,90]; % Elevation deviation in degrees

tx_phi = 0.00001;   % Transmitter azimuth angle in degrees
tx_theta = 0.0001;  % Transmitter elevation angle in degrees

rx_theta = 45;      % Receiver elevation angle in degrees
rx_phi = 0.0001;    % Receiver azimuth angle in degrees

a = 10 * lambda; % IRS element width in meters
b = 10 * lambda; % IRS element height in meters
N = 1;           % Number of IRS elements in the horizontal direction
M = 1;           % Number of IRS elements in the vertical direction

% Calculate IRS gain using the specified parameters
[IRS_gain] = IRS_model1_dev(a, b, N, M, tx_phi, tx_theta, phi_dev, theta_dev, rx_phi, rx_theta, lambda);

% Plot the results
plot(phi_dev, IRS_gain, 'LineWidth', 2)
xlabel('Azimuth Deviation [deg]');
ylabel('Gain at IRS model 1 [dB]')
legend('𝜽_{dev}=0','𝜽_{dev}=30','𝜽_{dev}=45','𝜽_{dev}=60','𝜽_{dev}=90')
