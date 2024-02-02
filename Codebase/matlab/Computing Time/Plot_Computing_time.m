% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, workspace, and close all figures
clc;
clear all;
close all;

% Specify the folder path and file names
folder = 'path_to_the_folder'; % Replace with the folder path
files = {'file1.csv', 'file2.csv', 'file3.csv'};

% Get the number of files
num_files = numel(files);

% Initialize an array to store means for each file
means = zeros(num_files, 1);

% Reads all the files and calculates the mean for the last column
for i = 1:num_files
    % Read data from the current file
    data = readtable(fullfile(folder, files{i}));
    
    % Convert the table data to an array
    column_data = table2array(data);
    
    % Calculate the mean of the last column for the first 54 rows
    means(i) = mean(column_data(1:54, end));
end

% Specify groups for the bar graph
groups = {'a=b=1cm', 'a=b=10cm', 'a=b=20cm'}; % Assign proper groups

% Create a vector with mean values for each group
values_Sim = [means(1); means(2); means(3)];

% Create a bar graph
figure;
bar(values_Sim);
set(gca, 'XTickLabel', groups);
xlabel('Reflector size a,b [m]');
ylabel('Time measurements [s]');
legend('2° Resolution', '1° Resolution', '0.5° Resolution');
colormap('parula'); % Change the color map if desired
grid on;
set(gcf, 'Position', [100, 100, 800, 500]); % Set figure position
