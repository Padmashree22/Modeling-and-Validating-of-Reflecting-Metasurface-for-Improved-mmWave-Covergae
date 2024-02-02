% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 


% Clear Command Window, Close All Figures, and Clear Workspace
clc;
close all;
clear all;

% Define the folder path where CSV files are located
folderPaths = 'folderPaths';

% Initialize loop counter
n = 1;

% Run the loop while n is less than or equal to 1
while n <= 1
    % List all CSV files in the specified folder
    fileList = dir(fullfile(folderPaths, '*.csv'));
    
    % Initialize cell array to store frames
    frames = cell(1, numel(fileList));
    
    % Loop through each CSV file
    for i = 1:numel(fileList)
        % Extract file name and full file path
        fileName = fileList(i).name;
        filePath = fullfile(folderPaths, fileName);
        
        % Check if the file exists
        assert(exist(filePath, 'file') == 2, 'File does not exist');
        
        % Load simulation data from the CSV file
        [RCS_mat, phi_vec, theta_vec] = load_simulation_data(filePath); % -> code appended
        
        % Calculate maximum RCS and mean RCS for each file
        RCS_Max(i) = max(RCS_mat(:));
        Mean_Sim(i) = mean(RCS_mat(:));
    end

    % Upper limit of RCS
    RCS_max = max(RCS_mat, [], 'all'); % dB
    RCS_max = ceil(RCS_max./10).*10;
    RCS_min = RCS_max - 50; % dB

    % Plotting: Maximum Gain vs. Theta_in
    figure;
    theta_in = 0:10:90;
    plot(theta_in, RCS_Max, 'LineWidth', 2);
    xlabel('Theta_{in} vary');
    ylabel('Maximum Gain [dB]')

    % Plotting: Heatmap Simulation
    figure;
    plot_heatmap_simulation(RCS_mat, phi_vec, theta_vec, ...
        'Gain [dB]', 'Phi [deg]', 'Theta [deg]', ...
        RCS_min, RCS_max); % -> code appended

    % Plotting: 2D Pattern Slice
    figure;
    theta_slice = 0;
    k = find(theta_vec == theta_slice);
    plot_pattern_slice(RCS_mat(:, k), phi_vec, ...
        'Gain [dB]', 'Azimuth angle at receiver [deg]', ...
        RCS_min, RCS_max, ...
        -90, 90, 'b');
    
    % Increment loop counter
    n = n + 1;
end
