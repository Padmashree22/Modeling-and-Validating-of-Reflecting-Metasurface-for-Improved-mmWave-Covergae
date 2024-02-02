% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, close all figures, and clear workspace variables
clc;
close all;
clear all;

% Add path to UMi_functions folder
addpath("UMi_functions\");

%% Variables and parameters

% Speed of light
c = physconst('lightspeed');

% 2D distance in meters; Source to Reflector
D_SR = 10:1:100;

% 2D distance in meters; Reflector to Destination
D_RD = 10:1:100;

% 2D distance in meters; Source to Destination
D_SD = 10:1:100;

% Operating frequency in GHz
freq = 28;

% Wavelength
lambda = c / freq;

% Heights of different components
h_BS_UMi = 10;  % Height of base station in UMi environment
h_BS_UMa = 25;  % Height of base station in UMa environment
h_Ref = 5;      % Height of reflector
h_UT = 1.5;     % Height of user terminal
h_E = 1;        % Height of environment (ground)

% Initialize an array to store PL_Radar_dB values for different Ref_location
PL_Radar_dB_values = zeros(size(D_SR));
shift_values = zeros(size(D_SR));

% Loop through different reflector locations
for Ref_location = 11:100
    
    % Call UMi_scenario function to calculate path loss for different cases
    [PL_UMi_LOS, PL_UMi_NLOS, PL_UMi_Radar, delta_cases] = UMi_scenario(freq, D_SR, D_RD, D_SD, h_BS_UMi, h_Ref, h_UT, h_E, c, Ref_location);
    
    % Calculate shift in PL_Radar_dB values for each reflector location
    shift = PL_UMi_Radar(D_SD == Ref_location+1) - PL_UMi_Radar(D_SD == Ref_location);
    
    % Store the PL_Radar_dB values for each reflector location
    PL_Radar_dB_values = [PL_Radar_dB_values; PL_UMi_Radar];
    shift_values = [shift_values, shift];
end

% Remove the initial zeros from the array
PL_Radar_dB_values = PL_Radar_dB_values(2:end, :);
shift_values = shift_values(:, 91:end);

% Plot the shift values against reflector distance
plot(11:100, shift_values, 'LineWidth', 2);
hold on;
xlabel('Reflector distance [m]');
ylabel('Shift [dB]');
