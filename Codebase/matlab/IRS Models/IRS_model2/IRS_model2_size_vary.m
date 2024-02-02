% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, workspace, and close all figures
clc;
clear all;
close all;

% Add path to the folder containing helper functions
addpath('IRS_model_functions\')

%% Parameters
c = 3e8;            % Speed of light in m/s
f = 28e9;           % Frequency in Hz
lambda = c / f;     % Wavelength in meters
gamma = 0.9;        % Reflection coefficient
rx_theta = 0.00000001;
rx_phi = -90:0.001:90;

% Compute reflection coefficient for various receiver angles
[rx_gamma] = gamma_rx(rx_theta, rx_phi);

tx_gamma = 90;

%% Variable Parameters
N_var = [1, (10:10:100)];

%% Plotting
figure;

y = 1; % Counter for legends
for N = N_var
    x = 1;
    for gamma_rx = rx_gamma
        % Dimensions of the RIS (Reflective Intelligent Surface)
        a = N * lambda; % in meters
        b = N * lambda; % in meters
        
        % Compute RIS gain
        RIS_gain1(x) = ((4 * pi) / (lambda^2))^2 * gamma^2 * (a * b)^2 * cosd(tx_gamma / sqrt(2)) * cosd(gamma_rx / sqrt(2));
        
        x = x + 1;
    end
    
    % Convert gain to dB
    RIS_gain1 = 10 * log10(abs(RIS_gain1));
    
    % Find minimum and maximum gain values
    RIS_gain1_min = min(min(RIS_gain1));
    RIS_gain1_max(y) = max(max(RIS_gain1));
    
    % Find peaks and locations of peaks in the gain
    [pks1, locs1] = findpeaks(RIS_gain1, rx_phi);
    [g, h] = size(pks1);
    
    % Calculate gain loss
    if h == 0 || h == 1
        gain_loss1(y) = 0;
    else
        [max_RIS1(y), I1(y)] = max(pks1);
        gain_loss1(y) = max_RIS1(y) - pks1(I1(y) + 1);
    end
    
    % Plot the gain
    plot(rx_phi, RIS_gain1, 'LineWidth', 2)
    hold on;
    
    xlabel('Azimuth at Receiver [deg]');
    ylabel('Gain at IRS model 2 [dB]');
    legends{y} = sprintf('a=b= %d lambda', N);
    
    y = y + 1;
end

legend(legends);
