function [RCS_mat, phi_vec, theta_vec] = load_simulation_data(file_name)
    %% Specify *csv file options
    opts = delimitedTextImportOptions("NumVariables", 3);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["Phideg", "Thetadeg", "dBRCSTotal"];
    opts.VariableTypes = ["double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    %% Load data from the specified CSV file
    sim_data = readtable(file_name, opts);
    
    %% Reshape data for further analysis
    phi_sim = sim_data{:, 1};            % Phi values in degrees
    theta_sim = sim_data{:, 2};          % Theta values in degrees
    RCS_sim = sim_data{:, 3};            % RCS values in dB
    
    % Unique phi and theta values
    phi_vec = unique(phi_sim);
    theta_vec = 90 - unique(theta_sim);  % Change coordinate system definition
    
    % Reshape RCS data into a matrix
    RCS_mat = reshape(RCS_sim, length(phi_vec), length(theta_vec));
    
    %% Find maximum RCS and corresponding angles
    [sigma_max, I] = max(RCS_sim);
    phi_max = phi_sim(I);
    theta_max = 90 - theta_sim(I); % Change coordinate system definition
    
    % Display maximum gain information
    disp(['Max. Gain (dB): ', num2str(round(sigma_max, 2))])
    disp(['At [phi, theta] (deg): [', ...
         num2str(round(phi_max, 2)), ', ', ...
         num2str(round(theta_max, 2)), ']'])
end
