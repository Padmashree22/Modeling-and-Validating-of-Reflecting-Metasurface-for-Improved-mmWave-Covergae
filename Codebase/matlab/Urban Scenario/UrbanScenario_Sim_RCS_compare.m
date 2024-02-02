% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, close all figures, and clear workspace
clc; close all; clear all;

% Add path to helper function
addpath('Helper_function')

% Load LOS grid data
load('LOS_scenario1_grid', 'LOS_grid')

%% Define simulation parameters
freq = 28e9;  % Frequency in Hz (28 GHz)
c = physconst('Light');  % Speed of light
tx_position = [309.9, 170.1, 10.0];  % Transmitter position [x, y, z] in meters 
reflector_center = [195, 170, 5];  % Reflector footprint [x, y, z] in meters
reflector_size = [3, 3];  % [length, width] in meters
h_E = 1;  % Effective environment height

% Transmitter and receiver parameters
P_TX = 20.3;  % Transmitter power in dBm
G_TX = 19.7;  % Transmitter gain in dBi
G_RX = 13.5;  % Receiver gain in dBi
R_Sensitivity = 83.5;  % Receiver sensitivity in dBm

receiver_x = 120:1:240;
receiver_y = 120:1:240;
receiver_z = 1.5;
[X, Y] = meshgrid(receiver_x, receiver_y);
Y = rot90(Y, 2);

phi_in = 80.17;  % Incident azimuth angle [deg]
theta_in = -2.46;  % Incident elevation angle [deg]

%% Define pattern resolution and limits
pattern_resolution = 0.25;  % Resolution of RCS [deg]
phi_out_min = -90;  % Minimum reflected azimuth wave [deg]
phi_out_max = 90;   % Maximum reflected azimuth wave [deg]
theta_out_min = -90;  % Minimum reflected elevation wave [deg]
theta_out_max = 90;   % Maximum reflected elevation wave [deg]

%% Define array parameters
N = 1:32;  % Modules along Y axis 
M = 1;  % Modules along Z axis

% Slope angle definitions
alpha = (27.6) * ones(length(M), length(N));  % Horizontal slope angle [deg]
alpha = flip(alpha);
beta = (0) * ones(length(M), length(N));  % Vertical slope angle [deg]

%% Calculate effective environment parameters
lambda = c / freq;  % Wavelength [m]
a_dim = 3 / 32;  % Dimensions along Y axis [m]
b_dim = 3;  % Dimensions along Z axis [m]
[sigma_sum_max, sigma_sum_mat, eta_a, eta_b] = fn_SF_HELIOS_array_eff(a_dim, b_dim, N, M, alpha, beta, phi_in, theta_in, phi_out, theta_out, lambda);
sigma_sum_mat = (10 * log10(sigma_sum_mat.^2));
mean_anal = mean(sigma_sum_mat(:));

%% Simulations
file_name = 'file_name.csv';
assert(exist(file_name, 'file') == 2, 'File does not exist')
[RCS_mat, phi_vec, theta_vec] = load_simulation_data(file_name);
RCS_Max = max(max(RCS_mat));
Mean_Sim = mean(RCS_mat(:));
RCS_max = max(RCS_mat, [], 'all');
RCS_max = ceil(RCS_max./10) * 10;
RCS_min = RCS_max - 50;

%% Plots
if length(phi_out) > 1 && length(theta_out) > 1
    subplot(1, 2, 1);
    plot_heatmap_analytical(sigma_sum_mat, phi_out, theta_out, 'RCS [dB]', 'Phi [deg]', 'Theta [deg]', -30, 55);
    ylim([-20 20]);
    xlim([-40 20]);
    subplot(1, 2, 2);
    plot_heatmap_simulation(RCS_mat, phi_vec, theta_vec, 'Gain [dB]', 'Phi [deg]', 'Theta [deg]', -30, 55);
    ylim([-20 20]);
    xlim([-40 20]);
else
    if length(phi_out) > 1
        plot_pattern_slice((10 * log10(sigma_sum_mat(:,:,end).^2)), phi_out, 'Gain [dB]', 'Azimuth angle at receiver [deg]', RCS_min, RCS_max, -90, 90, 'r');
        hold on;
        theta_slice = theta_in;
        k = find(theta_vec == theta_slice);
        plot_pattern_slice(RCS_mat(:,k), phi_vec, 'Gain [dB]', 'Azimuth angle at receiver [deg]', RCS_min, RCS_max, -90, 90, 'b-.');
        legend('RCS with shadowing factor', 'Simulation results');
    else
        plot_pattern_slice((10 * log10(sigma_sum_mat(:,:,end).^2)), theta_out, 'Gain [dB]', 'Elevation angle at receiver [deg]', RCS_min, RCS_max, -90, 90, 'r');
        hold on;
        phi_slice = phi_in;
        k = find(phi_vec == phi_slice);
        plot_pattern_slice(RCS_mat(k,:), theta_vec, 'Gain [dB]', 'Elevation angle at receiver [deg]', RCS_min, RCS_max, -90, 90, 'b-.');
        legend('RCS with shadowing factor', 'Simulation results');
    end
end

%% CDF plots
figure;
vector = sigma_sum_mat;
valid_values = vector(~isnan(vector) & vector ~= -inf);
sorted_values = sort(valid_values);
[f, x] = ksdensity(sorted_values);
cdf_values = cumsum(f) / sum(f);

vector_Sim = RCS_mat;
valid_values_Sim = vector_Sim(~isnan(vector_Sim) & vector_Sim ~= -inf);
sorted_values_Sim = sort(valid_values_Sim);
[fsim, xsim] = ksdensity(sorted_values_Sim);
cdf_values_Sim = cumsum(fsim) / sum(fsim);

plot(x, cdf_values, 'r', xsim, cdf_values_Sim, 'b', 'LineWidth', 2);
xlabel('Gain [dB]');
ylabel('ECDF');
hold on;
grid on;
legend('Analytical Narrow', 'Simulation Narrow');
xlim([-100 55]);
