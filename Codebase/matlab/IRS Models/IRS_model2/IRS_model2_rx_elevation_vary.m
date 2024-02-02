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
c = 3e8;         % Speed of light in m/s
f = 28e9;        % Frequency in Hz
lambda = c/f;    % Wavelength in meters
gamma = 0.9;     % Reflection coefficient

rx_theta = [0.0000001, 30, 45, 60, 90]; 
rx_phi = -90:0.01:90;
rx_gamma = gamma_rx(rx_theta, rx_phi);  % Function call to get gamma_rx values

tx_phi = 0.00001;
tx_theta = 0.0001;
tx_gamma = 0;
N = 50;

%% Plotting
figure;

y = 1;
for rx_theta_val = rx_theta
    x = 1;
    for gamma_rx_val = rx_gamma(y,:)
        % Set the dimensions of the RIS
        a = N * lambda; % in meters
        b = N * lambda; % in meters
        
        % Calculate RIS gain based on given formula
        RIS_gain1(x) = ((4 * pi) / (lambda^2))^2 * gamma^2 * (a * b)^2 * cosd(tx_gamma / sqrt(2)) * cosd(gamma_rx_val / sqrt(2));
        x = x + 1;
    end
    
    % Convert RIS gain to dB
    RIS_gain1 = 10*log10(abs(RIS_gain1));
    RIS_gain1_min = min(min(RIS_gain1));
    RIS_gain1_max(y) = max(max(RIS_gain1));
    
    % Find peaks in the gain profile
    [pks1, locs1] = findpeaks(RIS_gain1, rx_phi);
    [g, h] = size(pks1);
    
    % Calculate gain loss
    if h == 0 || h == 1
        gain_loss1(y) = 0;
    else
        [max_RIS1(y), I1(y)] = max(pks1);
        gain_loss1(y) = max_RIS1(y) - pks1(I1(y) + 1);
    end
    
    % Plot the gain profile
    plot(rx_phi, RIS_gain1, 'LineWidth', 2)
    hold on;
    
    % Increment for the legend
    y = y + 1;
end

% Labeling and legend
xlabel('Azimuth at Receiver [deg]');
ylabel('Gain at IRS model 2 [dB]')
legend('\theta_{rx}=0', '\theta_{rx}=30', '\theta_{rx}=45', '\theta_{rx}=60', '\theta_{rx}=90')
