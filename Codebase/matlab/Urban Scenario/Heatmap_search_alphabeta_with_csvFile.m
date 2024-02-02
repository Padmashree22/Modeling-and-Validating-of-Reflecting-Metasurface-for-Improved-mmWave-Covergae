% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear the command window, workspace, and close all figures
clc;
clear all;
close all;

% Add the path to the folder containing helper functions
addpath('Helper_function\');

% Define the output folder path and file name
outputFolder = 'output_folder_path';  % Replace with your desired folder path
outputFileName = 'outputFileName.csv';  % Replace with your desired file name

% Construct the full file path
filePath = fullfile(outputFolder, outputFileName);

% Read data from the CSV file into the 'data' matrix
data = csvread(filePath, 0, 0);

% Extract Mean_Power, alpha, and beta from the data
Mean_Power = data(:, 4);
alpha = unique(data(:, 2));
beta = unique(data(:, 3));

% Reshape Mean_Power into a 2D matrix based on unique alpha and beta values
Mean_Power = reshape(Mean_Power, length(alpha), length(beta));

% Find the maximum mean power and its corresponding indices
[max_mean_power, I] = max(Mean_Power(:));

% Display the corresponding alpha and beta values for the maximum mean power
disp(['alpha= ', num2str(data(I, 2)), ' and beta= ', num2str(data(I, 3))]);

% Create an image plot of Mean_Power with alpha on the x-axis and beta on the y-axis
imagesc(alpha, beta, Mean_Power);

% Set colormap to 'parula'
colormap(parula);

% Display colorbar
colorbar;

% Label x-axis and y-axis
xlabel('Alpha');  
ylabel('Beta');
