% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, workspace, and close all figures
clc;
clear all;
close all;

% Add path to helper functions
addpath('IRS_model_functions\') % load helper_function

%% Parameters

% Constants
c = 3e8;         % Speed of light in m/s
f = 28e9;        % Frequency in Hz
lambda = c/f;    % Wavelength in meters
gamma = 0.9;     % Reflection coefficient

% Receiver angles
rx_theta = 10;                % Elevation angle at the receiver
rx_phi = -90.1:0.5:90;        % Azimuth angles at the receiver
[rx_gamma] = gamma_rx(rx_theta, rx_phi);

% Transmitter angles
tx_phi = 0.00001;             % Azimuth angle at the transmitter
tx_theta = 0.0001;            % Elevation angle at the transmitter

% Dimensions of the IRS
a = lambda;      % in meters
b = lambda;      % in meters

% Number of reflecting elements
N_var = [1, (10:10:100)];   % Varying the number of elements
M = N_var;

% Plotting the results
figure;

% Loop over different values of N (number of elements)
y = 1;
for N = N_var
    x = 1;

    % Calculate IRS gain using the specified parameters
    [IRS_gain] = IRS_model1(a, b, N, M, tx_phi, tx_theta, rx_phi, rx_theta, lambda);

    % Plotting the results
    plot(rx_phi, IRS_gain, 'LineWidth', 2)
    hold on
    xlabel('Azimuth at Receiver [deg]');
    ylabel('Gain at IRS model 1 [dB]')
    legends{y} = sprintf('a=b= %d lambda', N);
    y = y + 1;
end

legend(legends)
