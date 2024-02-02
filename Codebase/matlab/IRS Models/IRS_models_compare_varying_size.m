% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, workspace, and close all figures
clc;
clear all;
close all;

% Add path to helper functions folder
addpath('IRS_model_functions\') % load helper_function

%% Parameters

% Constants
c = 3e8;        % Speed of light in m/s
f = 28e9;       % Frequency in Hz
lambda = c/f;   % Wavelength in meters
gamma = 0.9;    % Reflection coefficient

% Angles
rx_phi = -90:0.1:90;
rx_theta = 0.000001;
tx_phi = 0.00001;
tx_theta = 0.0001;
tx_alpha = 0.0001;

% IRS Dimensions
N = 1;
M = 1;
b1 = [1, (5:5:55)]; % Reflectors lengths in cm

%% IRS models

% Loop through different reflector lengths
for i = 1:length(b1)
    b(i) = b1(i);
    a(i) = b1(i) / N;
    
    % Compute gains for different IRS models
    [IRS_gain1] = IRS_model1(a(i)/100, b(i)/100, N, M, tx_phi, tx_theta, rx_phi, rx_theta, lambda);
    [IRS_gain2] = IRS_model2(a(i)/100, b(i)/100, N, M, tx_phi, tx_theta, rx_phi, rx_theta, lambda);
    [IRS_gain3] = IRS_model3(a(i)/100, b(i)/100, N, M, tx_phi, tx_theta, rx_phi, rx_theta, lambda);
    
    % Find maximum gain for each model
    IRS_gain_max1(i) = max(IRS_gain1);
    IRS_gain_max2(i) = max(IRS_gain2);
    IRS_gain_max3(i) = max(IRS_gain3);
end

% Plot the results
plot(b1, IRS_gain_max1, '-.', b1, IRS_gain_max2, '--', b1, IRS_gain_max2, ':', 'LineWidth', 2)
xlabel('Length of reflector (m)');
ylabel('Gain at IRS [dB]');
legend('Model 1 - Özdogan', 'Model 2 - Ntontin', 'Model 3 - Tang');
