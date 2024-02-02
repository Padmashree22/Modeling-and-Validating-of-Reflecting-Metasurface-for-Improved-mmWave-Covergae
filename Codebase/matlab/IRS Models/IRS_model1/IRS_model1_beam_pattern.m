% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, workspace, and close all figures
clc;
clear all;
close all;

% Add path to the folder containing IRS model functions
addpath('IRS_model_functions\') % load helper_function

%% Parameters

% Constants
c = 3e8;         % Speed of light in meters/second
f = 28e9;        % Frequency in Hertz
lambda = c/f;    % Wavelength in meters
gamma = 0.9;     % Reflection coefficient

% Deviations for theta and phi
theta_dev = 0.0000000001; 
phi_dev = -90:0.1:90;

% Transmit (TX) angles
tx_phi = 0.00001;
tx_theta = 0.0001;

% Receive (RX) angles
rx_theta = 45; 
rx_phi = 0.0001;

% Dimensions of the IRS element (assumed square for simplicity)
a = 0.366; % in centimeters
b = 0.366;

% Number of elements in the horizontal and vertical directions
N = [50];
M = N;

% Calculate IRS gain using the specified parameters
[IRS_gain] = IRS_model1_dev(a, b, N, M, tx_phi, tx_theta, phi_dev, theta_dev, rx_phi, rx_theta, lambda);

% Plot the results
plot(phi_dev, IRS_gain, 'LineWidth', 2)
hold on;
xlabel('α_{dev} variation [deg]');
ylabel('Gain at IRS model 1 [dB]')
legend(sprintf('𝛾_{Rx} %d', rx_theta));
