% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear the command window, workspace, and close all figures
clc; 
clear all; 
close all;

% Add the path to the folder containing helper functions
addpath('IRS_model_functions\') % Load helper_function

%% Parameters
c = 3e8;        % Speed of light in m/s
f = 28e9;       % Frequency in Hz
lambda = c / f; % Wavelength in meters
gamma = 0.9;    % Reflection coefficient
r_tx_RIS = 100; % Distance from transmitter to RIS in meters
r_rx_RIS = 100; % Distance from receiver to RIS in meters

tx_gamma = 60; % Transmitter reflection angle in degrees
rx_theta = 0.00000001; % Receiver elevation angle in degrees
rx_phi = -90:0.001:90;  % Receiver azimuth angles in degrees

% Compute the reflection coefficient at the receiver for varying azimuth angles
[rx_gamma] = gamma_rx(rx_theta, rx_phi);

% Compute transmitter and receiver antenna patterns
for x = 1:length(tx_gamma)
    F_t(x) = (cosd(tx_gamma(x) / sqrt(2))).^3;
end

for x = 1:length(rx_gamma)
    F_r(x) = (cosd(rx_gamma(x) / sqrt(2))).^3;
end

% Vary the number of RIS elements (N) and compute RIS gain
N_var = 10:10:50;
figure;
y = 1;

for N = N_var
    % RIS dimensions
    a = N * lambda; % in meters
    b = N * lambda; % in meters
    
    % Calculate RIS gain
    RIS_gain1 = ((4 * pi) / (lambda^2))^2 * gamma^2 * (a * b)^2 ...
        * abs(sqrt(F_t) .* sqrt(F_r) .* (exp((-1j * 2 * pi * (r_tx_RIS + r_rx_RIS)) / lambda))).^2;

    % Convert RIS gain to dB
    RIS_gain1 = 10 * log10(RIS_gain1);
    RIS_gain1_min = min(min(RIS_gain1));
    RIS_gain1_max(y) = max(max(RIS_gain1));

    % Plot RIS gain for each N
    plot(rx_phi, RIS_gain1, 'LineWidth', 2)
    hold on;
    
    xlabel('Azimuth at Receiver [degrees]');
    ylabel('Gain at RIS [dB]')
    legends{y} = sprintf('a=b= %d lambda', N);
    y = y + 1;
end

% Display legend
legend(legends)

