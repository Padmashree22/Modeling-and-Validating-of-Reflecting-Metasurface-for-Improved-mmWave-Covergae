% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, close all figures, and clear workspace
clc;
close all;
clear all;

% Add paths to required functions
addpath('Analytical_Model_functions'); % Load analytical model functions
addpath('Plotting functions\'); % Load plotting functions

%% Parameters

% Physical constants and frequency
c = physconst('Light'); % speed of light [m/s]
f = 28e9; % frequency [Hz]

% Incident angles
phi_in = 0; % Incident azimuth angle [deg]
theta_in = 0; % Incident elevation angle [deg]

%% Pattern resolution and output angles

% Resolution of RCS [deg]
pattern_resolution = 1;

% Range of reflected azimuth and elevation waves
phi_out_min = -90; % Minimum reflected azimuth wave [deg]
phi_out_max = 90; % Maximum reflected azimuth wave [deg]
theta_out_min = -90; % Minimum reflected elevation wave [deg]
theta_out_max = 90; % Maximum reflected elevation wave [deg]

% Generate vectors for reflected azimuth and elevation waves
theta_out = theta_out_min:pattern_resolution:theta_out_max;
phi_out = phi_out_min:pattern_resolution:phi_out_max;

% Calculate wavelength
lambda = c/f; % wavelength [m]

%% Array and slope parameters

% Array parameters
N = 1; % Modules along Y axis
M = 1:4; % Modules along Z axis
a = 0.4; % Dimensions along Y axis [m]
b = 0.1; % Dimensions along Z axis [m]

% Slope angles
alpha = (0.2) .* ones(length(M), length(N)); % Hor. Slope angle [deg]
beta = (0.1) .* ones(length(M), length(N)); % Vert. Slope angle [deg]

% Calculate effective RCS using the analytical model
[sigma_sum_max, sigma_sum_mat, eta_a, eta_b] = fn_SF_HELIOS_array_eff(a, b, N, M, alpha, beta, phi_in, theta_in, phi_out, theta_out, lambda);
sigma_sum_mat = (10 .* log10(sigma_sum_mat.^2));
Mean_RCS = mean(sigma_sum_mat(:));

%% Simulations

% Load simulation data from a CSV file
file_name = 'path_to_the_file.csv';
assert(exist(file_name, 'file') == 2, 'File does not exist')
data = readtable(fullfile(file_name)); % Import data from the CSV file
column1_data = table2array(data); % Extract the column

% Load simulation data into variables
[RCS_mat, phi_vec, theta_vec] = load_simulation_data(file_name);
RCS_Max = max(max(RCS_mat));
Mean_sim = mean(RCS_mat(:));

%% Upper limit of RCS

% Calculate upper limit of RCS for plotting
RCS_max = max(RCS_mat, [], 'all'); % dB
RCS_max = ceil(RCS_max./10) .* 10;
RCS_min = RCS_max - 50; % dB

%% Plots: Heatmap and 2D plots

figure;

if length(phi_out) > 1 && length(theta_out) > 1
    % Plot analytical heatmap
    subplot(1, 2, 1);
    plot_heatmap_analytical(sigma_sum_mat, phi_out, theta_out, ...
        'Analytical gain [dB]', 'Phi [deg]', 'Theta [deg]', ...
        RCS_min, RCS_max);
    ylim([-15 15]);
    xlim([-15 15]);

    % Plot simulation heatmap
    subplot(1, 2, 2);
    plot_heatmap_simulation(RCS_mat, phi_vec, theta_vec, ...
        'Simulation gain [dB]', 'Phi [deg]', 'Theta [deg]', ...
        RCS_min, RCS_max);
    ylim([-15 15]);
    xlim([-15 15]);

else
    if length(phi_out) > 1
        % Plot azimuth slices
        figure;
        plot_pattern_slice((10 .* log10(sigma_sum_mat(:,:,end).^2)), phi_out, ...
            'Gain [dB]', 'Azimuth angle at receiver [deg]', ...
            RCS_min, RCS_max, -90, 90, 'r');
        hold on;
        theta_slice = -theta_in; % deg
        k = find(theta_vec == theta_slice);
        plot_pattern_slice(RCS_mat(:, k), phi_vec, ...
            'Gain [dB]', 'Azimuth angle at receiver [deg]', ...
            RCS_min, RCS_max, ...
            -90, 90, 'b-.');
        clearvars k theta_slice
        legend('RCS with shadowing factor', 'Simulation results');

    else
        % Plot elevation slices
        plot_pattern_slice((10 .* log10(sigma_sum_mat(:,:,end).^2)), theta_out, ...
            'Gain [dB]', 'Elevation angle at receiver [deg]', ...
            RCS_min, RCS_max, -90, 90, 'r');
        hold on;
        phi_slice = phi_in; % deg
        k = find(phi_vec == phi_slice);
        plot_pattern_slice(RCS_mat(k, :), theta_vec, ...
            'Gain [dB]', 'Elevation angle at receiver [deg]', ...
            RCS_min, RCS_max, ...
            -90, 90, 'b-.');
        legend(['Phi = ', num2str(phi_slice), ' deg'])
        clearvars k phi_slice
        legend('RCS with shadowing factor', 'Simulation results');
    end
end

%% Violin plots

figure;
Data_all = column1_data(:, end) - sigma_sum_mat(:);
violin(Data_all, 'xlabel', {'4x1 Flatplate'}, 'edgecolor', 'k');
ylabel('Delta Gain [dB]');
ylim([-40 100]);
