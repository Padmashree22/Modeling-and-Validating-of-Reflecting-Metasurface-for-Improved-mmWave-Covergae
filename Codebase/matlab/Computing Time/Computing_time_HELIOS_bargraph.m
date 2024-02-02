% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear the command window, close all figures, and clear workspace
clc;
close all;
clear all;

% Add paths to helper functions and plotting functions
addpath('Analytical_Model_functions\')
addpath('Plotting_functions\')

% Specify folder paths for the computing time CSV files
folder1 = 'path_to_computing_time_analytical_csv_files'; % Replace with the folder path
folder2 = 'path_to_computing_time_128_core_simulation_csv_files'; % Replace with the folder path
folder3 = 'path_to_computing_time_4_core_simulation_csv_files'; % Replace with the folder path

% File names for analytical model, 128-core simulation, and 4-core simulation
RCS_files = {'file1.csv', 'file2.csv'};
Sim128_files = {'file3.csv', 'file4.csv'};
Sim4_files = {'file5.csv', 'file6.csv'};

% Get the number of files for each set
RCS_num_files = numel(RCS_files);
Sim128_num_files = numel(Sim128_files);
Sim4_num_files = numel(Sim4_files);

% Initialize vectors to store mean values for each set of files
RCS_means = zeros(RCS_num_files, 1);
Sim128_means = zeros(Sim128_num_files, 1);
Sim4_means = zeros(Sim4_num_files, 1);

% Calculate mean values for analytical model files
for i = 1:RCS_num_files
    data = readtable(fullfile(folder1, RCS_files{i}));
    column_data = table2array(data);
    RCS_means(i) = mean(column_data(:, end));
end

% Calculate mean values for 128-core simulation files
for i = 1:Sim128_num_files
    data = readtable(fullfile(folder2, Sim128_files{i}));
    column_data = table2array(data);
    Sim128_means(i) = mean(column_data(:, end));
end

% Calculate mean values for 4-core simulation files
for i = 1:Sim4_num_files
    data = readtable(fullfile(folder3, Sim4_files{i}));
    column_data = table2array(data);
    Sim4_means(i) = mean(column_data(:, end));
end

% Group names for the bar plot
groups = {'1x1 HELIOS','2x2 HELIOS','4x4 HELIOS'}; 

% Arrange mean values for bar plotting
Values_1=[RCS_means(2) RCS_means(1)];
Values_2=[Sim128_means(2) Sim128_means(1)];
Values_3=[Sim4_means(2) Sim4_means(1)];

% Combine all mean values for bar plotting
data_all=[Values_1; Values_2; Values_3];

% Create a bar plot
bar(data_all);
set(gca, 'XTickLabel', groups);
xlabel('Processors');
ylabel('Computing Time [s] ');
legend('2° Resolution', '1° Resolution');
colormap('parula'); % Change the color map if desired
grid on;
