% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear the command window, workspace, and close all figures
clc; clear all; close all;

% Add paths to helper functions
addpath('UMi_functions\') % load helper_function
addpath('UMa_functions\') % load helper_function
addpath('FSPL_function\') % load helper_function

%% Constants
c = physconst('lightspeed'); % Speed of light in meters per second
D_SR = 10:0.1:100;           % 2D distance in meters; Source to reflector
D_RD = 10:0.1:100;           % 2D distance in meters; reflector to Destination
D_SD = 10:0.1:100;           % 2D distance in meters; Source to Destination
freq = 28;                    % Frequency in GHz
lambda = c / freq;            % Wavelength [m/s]
h_BS_UMi = 10;                % Antenna height for BS in meters
h_BS_UMa = 25;                % Antenna height for BS in meters
h_Ref = 5;                    % Height of reflector in meters
h_UT = 1.5;                   % Antenna height for UT in meters
h_E = 1;                      % Effective environment height
Ref_location = 50;

%% FSPL Equations: General and Radar
% Calculate Free-Space Path Loss (FSPL) for LOS and Radar scenarios
[PL_LOS_dB, PL_Radar_dB] = Fn_FSPL(freq, D_SR, D_RD, D_SD, c, Ref_location);

% Plot FSPL results
figure;
plot(D_SD, PL_LOS_dB, 'r', D_SD, PL_Radar_dB, 'b--', 'LineWidth', 2);
xlabel('Distance [m]');
ylabel('Path Loss [dB]');
legend('FSPL LOS', 'FSPL NLOS', 'FSPL Radar Eq');

%% 3GPP UMi: General and Radar
% Calculate UMi (Urban Macrocell) path loss for LOS, NLOS, and Radar scenarios
[PL_UMi_LOS, PL_UMi_NLOS, PL_UMi_Radar] = UMi_scenario(freq, D_SR, D_RD, D_SD, h_BS_UMi, h_Ref, h_UT, h_E, c, Ref_location);

% Plot UMi results
figure;
plot(D_SD, PL_UMi_LOS, 'r', D_SD, PL_UMi_NLOS, 'g', D_SD, PL_UMi_Radar, 'b--', 'LineWidth', 2);
xlabel('Distance [m]');
ylabel('Path Loss [dB]');
legend('UMi LOS', 'UMi NLOS', 'UMi Radar');

%% 3GPP UMa: General and Radar Eq
% Calculate UMa (Urban Macrocell) path loss for LOS, NLOS, and Radar scenarios
[PL_UMa_LOS, PL_UMa_NLOS, PL_UMa_Radar] = UMa_scenario(freq, D_SR, D_RD, D_SD, h_BS_UMa, h_Ref, h_UT, h_E, c, Ref_location);

% Plot UMa results
figure;
plot(D_SD, PL_UMa_LOS, 'r', D_SD, PL_UMa_NLOS, 'g', D_SD, PL_UMa_Radar, 'b--', 'LineWidth', 2);
xlabel('Distance [m]');
ylabel('Path Loss [dB]');
legend('UMa LOS', 'UMa NLOS', 'UMa Radar');
