% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, workspace, and close all figures
clc;
clear all;
close all;

% Add path to the folder containing helper functions
addpath('IRS_model_functions\') % Load helper_function

%% Parameters
c = 3e8;            % Speed of light in m/s
f = 28e9;           % Frequency in Hz
lambda = c / f;     % Wavelength in meters
gamma = 0.9;        % Reflection coefficient
r_tx_RIS = 10;      % Distance from transmitter to RIS in meters
r_rx_RIS = 10;      % Distance from receiver to RIS in meters

tx_gamma = 0;               % Transmitter reflection coefficient
tx_theta = 0.00000001;      % Transmitter elevation angle in radians
tx_phi = 0.00000001;        % Transmitter azimuth angle in radians

rx_theta1 = [0.00000001, 30, 45, 60, 90]; % Receiver elevation angles in degrees
rx_phi = -90:0.001:90;      % Receiver azimuth angles in degrees
N_var = 50;                 % Variable for loop
figure;

z = 1; % Counter for subplot
for rx_theta = rx_theta1
    % Calculate receiver reflection coefficient based on elevation angle
    [rx_gamma] = gamma_rx(rx_theta, rx_phi);
    
    y = 1; % Counter for legend
    for N = N_var
        a = N * lambda; % Length in meters
        b = N * lambda; % Width in meters
        M = N;
        
        % Calculate IRS gain using the specified parameters
        [IRS_gain] = IRS_model3(a, b, N, M, tx_phi, tx_theta, rx_phi, rx_theta, lambda);
        
        % Plot the IRS gain
        plot(rx_phi, IRS_gain, 'LineWidth', 2)
        hold on;
        xlabel('Azimuth at Receiver [degrees]');
        ylabel('Gain at IRS Model 3 [dB]');
        legend_str{y} = ['\theta_{rx}=' num2str(rx_theta)]; % Create legend string
        y = y + 1;
    end
    
    z = z + 1;
end

% Display legend with appropriate labels
legend(legend_str);

% Set the title for the plot
title('IRS Model 3 Gain Variation');

% Reset the path to the original path
rmpath('IRS_model_functions\');
