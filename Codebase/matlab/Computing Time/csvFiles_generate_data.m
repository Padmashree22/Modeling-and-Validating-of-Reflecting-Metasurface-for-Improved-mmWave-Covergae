% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, close all figures, and clear all variables
clc;
close all;
clear all;

% Add path to helper functions folder
addpath('Analytical_Model_functions\') % load helper_functions

%% Parameters
c = physconst('Light'); % speed of light [m/s]
f = 28e9; % frequency [Hz]

%% Define parameters for RCS pattern
pattern_resolution = 2; % resolution of RCS [deg]
phi_out_min = -90; % Minimum reflected azimuth wave [deg]
phi_out_max = 90; % Maximum reflected azimuth wave [deg]
theta_out_min = -90; % Minimum reflected elevation wave [deg]
theta_out_max = 90; % Maximum reflected elevation wave [deg]

%% Dimensions and modules configuration
N = 1:4; % Modules along Y axis 
M = 1:4; % Modules along Z axis
a = 0.4; % Dimensions along Y axis [m]
b = 0.4; % Dimensions along Z axis [m]

%% Define angles for reflected waves and calculate wavelength
phi_out = phi_out_min:pattern_resolution:phi_out_max; % Reflected azimuth wave [deg]
theta_out = theta_out_min:pattern_resolution:theta_out_max; % Reflected elevation wave [deg]
lambda = c/f; % wavelength [m]

% Upper limit of RCS
RCS_max = 10*log10(4*pi*(a*b/lambda)^2); % [dB]
RCS_max = ceil(RCS_max/10)*10;
RCS_min = RCS_max - 50; % [dB]

%% File handling
folderPaths = "path_to_current_folder";  % path to the existing file
fileList = dir(fullfile(folderPaths, '*.csv'));
fileName = "file_name.csv";     % file name
filePath = fullfile(folderPaths, fileName); % Construct the full file path
data = csvread(filePath, 1, 0); % Read the CSV file
l = size(data,1);

% Preallocate variables
sigma_sum_max = zeros(l, 1);
Analytical_time = zeros(l, 1);

% Define the folder and filename for the new CSV file
outputFolder = 'path_to_new_folder';  % Replace with your desired folder path
outputFileName = 'output_File_Name.csv';  % Replace with your desired file name

% Extract the simulation data from the existing csv file and calculate analytical
% calculations and save them in a new csv file.
for i = 1:l
    alpha = data(i,1)*ones(length(M),length(N));
    beta = data(i,2)*ones(length(M),length(N));
    phi_in = data(i,3);
    theta_in = data(i,4) - 90;

    % Perform analytical calculations using the helper function
    tic;
    [sigma_sum_max(i), ~, ~, ~] = fn_SF_HELIOS_array_eff(a, b, N, M, alpha, beta, phi_in, theta_in, phi_out, theta_out, lambda);
    Analytical_time(i) = toc;

    % Write results to a new CSV file
    if i == 1
        newData = table(sigma_sum_max(i), Analytical_time(i), 'VariableNames', {'sigma_sum_max', 'Analytical_time'});
        writetable(newData, fullfile(outputFolder, outputFileName));
    else
        newData = table(sigma_sum_max(i), Analytical_time(i), 'VariableNames', {'sigma_sum_max', 'Analytical_time'});
        writetable(newData, fullfile(outputFolder, outputFileName),'WriteMode','Append',...
        'WriteVariableNames',false,'WriteRowNames',true);
    end
end
