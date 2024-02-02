% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, close all figures, and clear workspace
clc; 
close all; 
clear all;

% Add path to helper functions
addpath('Helper_function')

% Load data files
load('LOS_scenario1_grid','LOS_grid')
load('NLOS_region')

%% Define simulation parameters

% Communication parameters
freq = 28e9;  % Frequency in Hz (28 GHz)
c = physconst('Light');  % Speed of light

% Transmitter and reflector positions
tx_position = [309.9, 170.1, 10.0];  
reflector_center = [191.5, 170, 5];  
reflector_size = [3, 3];  

% Effective environment height
h_E = 1;          

% Transmitter and receiver hardware parameters
P_TX = 20.3;    % Transmit power in dBm
G_TX = 19.7;    % Transmit antenna gain in dBi
G_RX = 13.5;    % Receive antenna gain in dBi
R_Sensitivity = 83.5;  % Receiver sensitivity in dBm

% Grid for receiver
receiver_x = 120:1:240;
receiver_y = 120:1:240;
receiver_z = 1.5;
[X, Y] = meshgrid(receiver_x, receiver_y);
Y = rot90(Y,2);

% Incident angles
phi_in = 80.17;    % Incident azimuth angle [deg]
theta_in = -2.46;   % Incident elevation angle [deg]

%% Define pattern resolution and angles

% Pattern resolution
pattern_resolution = 0.1; 

% Reflected wave angles
phi_out_min = -90; 
phi_out_max = 90; 
theta_out_min = -90; 
theta_out_max = 90;

% Array parameters
N = 1:32; 
M = 1; 
a_dim = 3/32; 
b_dim = 3;

% Angle parameters
theta_out = theta_out_min:pattern_resolution:theta_out_max; 
phi_out = phi_out_min:pattern_resolution:phi_out_max; 
lambda = c/freq; 

%% Compute distances

% Distances between transmitter, receiver, and reflector
distance_tx_rx = sqrt((X - tx_position(1)).^2 + (Y - tx_position(2)).^2 + (receiver_z - tx_position(3)).^2);
distance_tx_ref = sqrt((reflector_center(1) - tx_position(1)).^2 + (reflector_center(2) - tx_position(2)).^2 + (reflector_center(3) - tx_position(3)).^2);
distance_ref_rx = sqrt((X - reflector_center(1)).^2 + (Y - reflector_center(2)).^2 + (receiver_z - reflector_center(3)).^2);

% Compute angles for reflection
phi_rx = round(asind((X - reflector_center(1))./distance_ref_rx), 1);
theta_rx = round(asind((reflector_center(2) - Y)./distance_ref_rx), 1);

a1 = abs(round((phi_out_max + phi_rx)./pattern_resolution));
a2 = abs(round((theta_out_max + theta_rx)./pattern_resolution));
a1(isnan(a1)) = 0;
a2(isnan(a2)) = 0;

var = 2.5;

% Calculate path loss
[PL_LOS_dB, PL_NLOS_dB, PL_Radar_dB] = Fn_FSPL(freq.'*1e-9, (distance_tx_ref), (distance_ref_rx),(distance_tx_rx), c, distance_tx_ref);
NLOS_region_case_sens1 = NLOS_region_case_sens(1:49,52:70);

% Output folder and file names
outputFolder = 'Output_folder_path';  
outputFileName = 'file_name.csv';  

i = 1;

% Loop over alpha and beta values
for alpha_var = -6.5
    for beta_var = -4.6
        tStart = tic; 

        %% Narrow beams
        % beta= (beta_var).*ones(length(M),length(N)); % Vert. Slope angle [deg]
        % alpha=(alpha_var).*ones(length(M),length(N)); % Hor. Slope angle [deg]

        %% Broad beams
        alpha = [(alpha_var)*(ones(1,2)) (alpha_var+var)*(ones(1,4)) (alpha_var+2*var)*(ones(1,5)) (alpha_var+3*var)*(ones(1,6)) (alpha_var+4*var)*(ones(1,6)) (alpha_var+5*var)*(ones(1,5)) (alpha_var+6*var)*(ones(1,4))];
        alpha = flip(alpha);
        beta = [(beta_var)*(ones(1,2)) (beta_var+var)*(ones(1,4)) (beta_var+2*var)*(ones(1,5)) (beta_var+3*var)*(ones(1,6)) (beta_var+4*var)*(ones(1,6)) (beta_var+5*var)*(ones(1,5)) (beta_var+6*var)*(ones(1,4))];
        beta = flip(beta);

        %% Gain_RCS
        [sigma_sum_max(i), sigma_sum_mat, eta_a, eta_b] = fn_SF_HELIOS_array_eff(a_dim,b_dim,N,M,alpha,beta,phi_in,theta_in,phi_out,theta_out,lambda);
        sigma_sum_mat = (10.*log10(sigma_sum_mat.^2));

        %% Different IRS models instead of HELIOS
        % [sigma_sum_mat] = RIS_model1(a_dim,b_dim,N,M,phi_in,theta_in,phi_out,theta_out,lambda);
        % [sigma_sum_mat] =RIS_model2(a_dim,b_dim,N,M,phi_in,theta_in,phi_out,theta_out,lambda);
        % [sigma_sum_mat] = RIS_model3(a_dim,b_dim,N,M,phi_in,theta_in,phi_out,theta_out,lambda);
        %     theta_rx=-25.75;
        %     phi_rx=-2.5;
        % [sigma_sum_mat]= RIS_model1_dev(a_dim,b_dim,N,M,phi_in,theta_in,phi_out,theta_out,theta_rx,phi_rx,lambda);

        %%
        j = 1;
        for f = 1:1:49
            j = 1;
            for g = 52:1:70
                PL2(f,j) = PL_Radar_dB(f,g);
                gain_rcs2(f,j) = sigma_sum_mat(a1(f,g)+1,a2(f,g)+1);
                P_RX2(f,j) = P_TX + G_TX + G_RX + gain_rcs2(f,j) - PL2(f,j);       
                Combined_P_RX(f,j) = P_RX2(f,j);
                if reflector_center(1) == 195 && NLOS_region_case_sens(f,g) == 1
                    gain_rcs2(f,j) = -1000; 
                    Combined_P_RX(f,j) = -1000; 
                end
                j = j+1;
            end
        end
        tEnd = toc(tStart);
        Mean_gain = mean(gain_rcs2(:));
        Mean_prx2 = mean(Combined_P_RX(:));
        
        % Write results to a CSV file
        newData = table(i, alpha_var, beta_var, Mean_gain, Mean_prx2, tEnd, 'VariableNames', {'Number','Alpha','Beta','Mean Gain','Mean Prx2','Run Time'});
        writetable(newData, fullfile(outputFolder, outputFileName),'WriteMode','Append',...
        'WriteVariableNames',false,'WriteRowNames',true); 
        i = i+1;
        clear vars alpha beta Mean_power Mean_gain Mean_prx2
    end
end
