% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, workspace, and close all figures
clc;
clear all;
close all;

% Add path to the folder containing the FSPL function
addpath("FSPL_function\");

%% Variables and parameters
c = physconst('lightspeed'); % Speed of light
D_SR = 1:1:100;    % 2D distance in meters: Source to reflector
D_RD = 1:1:100;    % 2D distance in meters: reflector to Destination
D_SD = 1:1:100;     % 2D distance in meters: Source to Destination
freq = 28;          % Frequency in GHz
lambda = c / freq;  % Wavelength

% Heights of different entities
h_BS_UMi = 10; % Height of the Base Station in UMi (Urban Micro) environment
h_BS_UMa = 25; % Height of the Base Station in UMa (Urban Macro) environment
h_Ref = 5;     % Height of the reflector
h_UT = 1.5;    % Height of the User Terminal
h_E = 1;       % Height of the Environment

% Initialize arrays to store path loss and shift values
PL_Radar_dB_values = zeros(size(D_SR)); % Array to store PL_Radar_dB values for different Ref_location
shift_values = zeros(size(D_SR));       % Array to store shift values

% Loop over different reflector locations
for Ref_location = 1:100
    % Call the FSPL function to calculate path loss values
    [PL_LOS_dB, PL_NLOS_dB, PL_Radar_dB] = Fn_FSPL(freq, D_SR, D_RD, D_SD, c, Ref_location);
    
    % Calculate shift in PL_Radar_dB for each Ref_location
    shift = PL_Radar_dB(D_SD == Ref_location+1) - PL_Radar_dB(D_SD == Ref_location);
    
    % Store the PL_Radar_dB values for each Ref_location
    PL_Radar_dB_values = [PL_Radar_dB_values; PL_Radar_dB];
    shift_values = [shift_values, shift];
end

% Remove the initial zeros from the array
PL_Radar_dB_values = PL_Radar_dB_values(2:end, :);
shift_values = shift_values(:, 100:end);

% Plot the shift values against reflector distance
plot(1:100, shift_values, 'LineWidth', 2);
xlabel('Reflector distance [m]');
ylabel('Shift [dB]');
