% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear Command Window, Close All Figures, and Clear Workspace
clc;
close all;
clear all;

% Add paths to helper functions and plotting functions
addpath('Analytical_Model_functions\') % load helper_function
addpath('Plotting_functions\')

% Load LOS grid data
load('LOS_scenario1_grid','LOS_grid')

%% Define simulation parameters

freq = 28e9;  % Frequency in Hz (28 GHz)
c = physconst('Light');  % Speed of light

% Transmitter position [x, y, z] in meters 
tx_position = [309.9, 170.1, 10.0]; 

% Reflector footprint center [x, y, z] in meters
reflector_center = [195, 170, 5]; 

% Reflector size [length, width] in meters
reflector_size = [3, 3];  

h_E=1;  % Effective environment height

% Parameters related to power and gain
P_TX= 20.3;         % Transmit power (dBm)
G_TX= 19.7;         % Transmit antenna gain (dBi)
G_RX= 13.5;         % Receive antenna gain (dBi)
R_Sensitivity=83.5; % Receiver sensitivity (dBm)

% Grid receiver coordinates
receiver_x= 120:1:240;
receiver_y= 120:1:240;
receiver_z=1.5;

% Create meshgrid for receiver coordinates
[X, Y] = meshgrid(receiver_x, receiver_y);
Y=rot90(Y,2); % Rotate Y for correct orientation

phi_in = 0;     % Incident azimuth angle [deg]
theta_in = 0;   % Incident elevation angle [deg]

%% Define parameters for RCS computation

pattern_resolution = 0.1; % Resolution of RCS [deg]

% Define ranges for reflected waves
phi_out_min=-90; 
phi_out_max=90; 
theta_out_min=-90; 
theta_out_max=90; 

% Number of modules along Y and Z axes
N=1; 
M=1; 

% Dimensions of modules along Y and Z axes [m]
a_dim=0.1; 
b_dim=0.1; 

% Horizontal and vertical slope angles [deg]
alpha=(10).*ones(length(M),length(N)); 
beta= 10.*ones(length(M),length(N)); 

% Compute reflected elevation and azimuth waves
theta_out = theta_out_min:pattern_resolution:theta_out_max; 
phi_out = phi_out_min:pattern_resolution:phi_out_max; 

% Compute wavelength [m]
lambda=c/freq; 

%% Calculate Gain and RCS

[sigma_sum_max_n, sigma_sum_mat_n, eta_a_n, eta_b_n] = fn_SF_HELIOS_array_eff(a_dim,b_dim,N,M,alpha,beta,phi_in,theta_in,phi_out,theta_out,lambda);
sigma_sum_mat_n=(10.*log10(sigma_sum_mat_n.^2));

% Extract valid RCS values
vector_RCS = sigma_sum_mat_n(:);
valid_values_RCS = vector_RCS(~isnan(vector_RCS) & vector_RCS ~= -inf);  

%% Load simulation data and perform kernel density estimation

file_name = ['file_name.csv']; % Change to real simulation file name
assert(exist(file_name, 'file') == 2, 'File does not exist')

% Load simulation data
[RCS_mat, phi_vec, theta_vec] = load_simulation_data(file_name);

% Extract valid simulation data values
vector_SIM = RCS_mat(:);
valid_values_SIM = vector_SIM(~isnan(vector_SIM) & vector_SIM ~= -inf);  

% Sort RCS and simulation data values
sorted_values_RCS = sort(valid_values_RCS);
sorted_values_SIM = sort(valid_values_SIM);

% Perform kernel density estimation
[f_RCS, x_RCS] = ksdensity(sorted_values_RCS);
[f_sim, x_sim] = ksdensity(sorted_values_SIM); 

% Calculate cumulative distribution functions (CDF)
cdf_values_RCS = cumsum(f_RCS) / sum(f_RCS); 
cdf_values_SIM = cumsum(f_sim) / sum(f_sim); 

% Plot empirical cumulative distribution functions (ECDF)
plot(x_RCS, cdf_values_RCS, 'r--',x_sim, cdf_values_SIM, 'b--','LineWidth',1.5);
xlabel('Gain [dB]');
ylabel('ECDF');
hold on;
grid on;
xlim([-100 60]);
ylim([0 1]);
