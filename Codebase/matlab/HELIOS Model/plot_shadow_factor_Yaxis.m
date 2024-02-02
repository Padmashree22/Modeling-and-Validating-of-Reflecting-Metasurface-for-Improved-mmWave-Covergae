% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, close all figures, and clear workspace
clc;
close all;
clear all;

% Add paths to helper functions
addpath('Analytical_Model_functions\') % load helper_function
addpath('Plotting_functions\')

% Define parameters
N = 1:3;     % Number of reflectors
M = 1;       
a = 0.1;     
b = 0.1;     

% Set alpha and beta slope angle constants
alpha = 10 .* ones(length(M), length(N));
beta = 0 .* ones(length(M), length(N));

% Set pattern resolution in degrees
pattern_resolution = 0.1; 

% Set input and output angles in degrees
phi_in = -90:1:90;
theta_in = 0;
phi_out = 0;
theta_out = -90:pattern_resolution:90;

% Define electromagnetic parameters
c = physconst('Light');   % Speed of light (m/s)
f = 28e9;                 % Frequency (Hz)
lambda = c / f;           % Wavelength (m)

% Calculate upper limit of RCS
RCS_max = 10 .* log10(4 .* pi .* (a .* b ./ lambda).^2); % dB
RCS_max = ceil(RCS_max ./ 10) .* 10;
RCS_min = RCS_max - 50; % dB

% Loop over input angles
y = 1;
for x = phi_in 
    theta_in = 0;
    % Call helper function to calculate effective parameters
    [a_eff(:,:,y), b_eff(:,:,y), eta_a(:,:,y),  eta_b(:,:,y)] = fn_ShadowFactor_eff(a, b, N, M, alpha, beta, x, theta_in);
    y = y + 1;
end

% Plot the results
figure;
plot(phi_in, (reshape(eta_a(:,1,:).*100,1,length(phi_in))), 'r', 'LineWidth', 2);
hold on;
plot(phi_in, (reshape(eta_a(:,2,:).*100,1,length(phi_in))), 'g-.', 'LineWidth', 2);
hold on;
plot(phi_in, (reshape(eta_a(:,3,:).*100,1,length(phi_in))), 'b--', 'LineWidth', 2);
xlabel('phi_{in} [deg]');
ylabel('𝜂 % ');
legend('Reflector 1','Reflector 2','Reflector 3');
