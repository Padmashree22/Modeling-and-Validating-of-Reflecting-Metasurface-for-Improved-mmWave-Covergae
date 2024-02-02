% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, close all figures, and clear workspace
clc; 
close all; 
clear all;

% Add paths for helper functions and plotting functions
addpath('Analytical_Model_functions\') % Load helper function
addpath('Plotting_functions\')

% Define parameters
N = 1;
M = 1:3;
a = 0.1;
b = 0.1;

% Set angles for incident and outgoing waves
phi_in = 0; % Incident angle in degrees
theta_in = -90:1:90; % Incident angle range in degrees
phi_out = 0; % Outgoing angle in degrees
pattern_resolution = 0.1; % Resolution for outgoing angle pattern in degrees
theta_out = -90:pattern_resolution:90; % Outgoing angle range in degrees

% Parameters for electromagnetic waves
c = physconst('Light'); % Speed of light in m/s
f = 28e9; % Frequency in Hz
lambda = c / f; % Wavelength in meters

% Calculate upper limit of Radar Cross Section (RCS)
RCS_max = 10 * log10(4 * pi * (a * b / lambda)^2); % dB
RCS_max = ceil(RCS_max / 10) * 10; % Round up to the nearest 10
RCS_min = RCS_max - 50; % dB

% Initialize slope angle variables for efficiency calculation
alpha = 0 .* ones(length(M), length(N));
beta = 10 .* ones(length(M), length(N));

% Loop over incident angles and calculate shadow factors
y = 1;
for x = theta_in
    phi_in = 0;
    [a_eff(:,:,y), b_eff(:,:,y), eta_a(:,:,y),  eta_b(:,:,y)] = fn_ShadowFactor_eff(a, b, N, M, alpha, beta, phi_in, x);
    y = y + 1;
end

% Plot efficiency as a function of incident angle
figure;
plot(theta_in, reshape(eta_b(1,:,:) * 100, 1, length(theta_in)), 'r', 'LineWidth', 2);
hold on
plot(theta_in, reshape(eta_b(2,:,:) * 100, 1, length(theta_in)), 'g-.', 'LineWidth', 2);
hold on
plot(theta_in, reshape(eta_b(3,:,:) * 100, 1, length(theta_in)), 'b--', 'LineWidth', 2);
xlabel('theta_{in} [deg]');
ylabel('\eta % ');
legend('Reflector 1', 'Reflector 2', 'Reflector 3');
