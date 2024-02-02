% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, close all figures, and clear workspace
clc;
close all;
clear all;

% Add path to helper functions and load necessary data
addpath('Helper_function') % Load helper function
load('LOS_scenario1_grid','LOS_grid') % Load LOS grid data
load('NLOS_region') % Load NLOS region data

%% Define simulation parameters

% Communication frequency
freq = 28e9;  % Frequency in Hz (28 GHz)

% Speed of light
c = physconst('Light');

% Transmitter position [x, y, z] in meters 
tx_position = [309.9, 170.1, 10.0];

% Reflector footprint [x, y, z] in meters
reflector_center = [191.5, 170, 5];
reflector_size = [3, 3];  % [length, width] in meters

% Effective environment height
h_E = 1;

% Transmitter parameters
P_TX = 20.3; % dBm
G_TX = 19.7; % dBi
G_RX = 13.5; % dBi
R_Sensitivity = 83.5; % dBm

% Grid for receiver positions
receiver_x = 120:1:240;
receiver_y = 120:1:240;
receiver_z = 1.5;
[X, Y] = meshgrid(receiver_x, receiver_y);
Y = rot90(Y,2);

% Incident angles
phi_in = 80.17; % Incident azimuth angle [deg]
theta_in = -2.46; % Incident elevation angle [deg]

%% Pattern Resolution

% Resolution of RCS [deg]
pattern_resolution = 0.1;

% Define range for reflected azimuth and elevation waves
phi_out_min = -90; % Minimum reflected azimuth wave [deg] 
phi_out_max = 90;  % Maximum reflected azimuth wave [deg]
theta_out_min = -90; % Minimum reflected elevation wave [deg]
theta_out_max = 90; % Maximum reflected elevation wave [deg]

%% IRS Module Configuration

% Number of modules along Y and Z axes
N = 1:32; % Modules along Y axis 
M = 1;    % Modules along Z axis

% Dimensions along Y and Z axes [m]
a_dim = 3/32; % Dimensions along Y axis [m]
b_dim = 3;    % Dimensions along Z axis [m]

% Slope angle variation
alpha_var = 10.4;
beta_var = 2.2;
var = 2.5;

%% (Uncomment the following lines if you are calculating narrow beam)
% alpha=(alpha_var).*ones(length(M),length(N)); % Hor. Slope angle [deg]
% beta= (beta_var).*ones(length(M),length(N)); % Vert. Slope angle [deg]

%% Broad beam calculations
alpha=[(alpha_var)*(ones(1,2)) (alpha_var+var)*(ones(1,4)) (alpha_var+2*var)*(ones(1,5)) (alpha_var+3*var)*(ones(1,6)) (alpha_var+4*var)*(ones(1,6)) (alpha_var+5*var)*(ones(1,5)) (alpha_var+6*var)*(ones(1,4))];
alpha=flip(alpha);

beta=[(beta_var)*(ones(1,2)) (beta_var+var)*(ones(1,4)) (beta_var+2*var)*(ones(1,5)) (beta_var+3*var)*(ones(1,6)) (beta_var+4*var)*(ones(1,6)) (beta_var+5*var)*(ones(1,5)) (beta_var+6*var)*(ones(1,4))];
beta=flip(beta);
%% Reflected Waves Calculation

% Define range for reflected elevation and azimuth waves
theta_out = theta_out_min:pattern_resolution:theta_out_max;
phi_out = phi_out_min:pattern_resolution:phi_out_max;

% Wavelength [m]
lambda = c / freq;

%% Distances Calculation

% Compute distances
distance_tx_rx = sqrt((X - tx_position(1)).^2 + (Y - tx_position(2)).^2 + (receiver_z - tx_position(3)).^2);
distance_tx_ref = sqrt((reflector_center(1) - tx_position(1)).^2 + (reflector_center(2) - tx_position(2)).^2 + (reflector_center(3) - tx_position(3)).^2);
distance_ref_rx = sqrt((X - reflector_center(1)).^2 + (Y - reflector_center(2)).^2 + (receiver_z - reflector_center(3)).^2);

% Far-field distance
FF_dist = max(((a_dim / 10).^2), ((b_dim / 10).^2)) / (lambda);

%% Gain_RCS Calculation

% Call the function to compute effective RCS
[sigma_sum_max, sigma_sum_mat, eta_a, eta_b] = fn_SF_HELIOS_array_eff(a_dim, b_dim, N, M, alpha, beta, phi_in, theta_in, phi_out, theta_out, lambda);
sigma_sum_mat = (10 .* log10(sigma_sum_mat.^2));

% Uncomment and replace with the desired IRS models
% [sigma_sum_mat] = IRS_model1(a_dim, b_dim, N, M, phi_in, theta_in, phi_out, theta_out, lambda);
% [sigma_sum_mat] = IRS_model2(a_dim, b_dim, N, M, phi_in, theta_in, phi_out, theta_out, lambda);
% [sigma_sum_mat] = IRS_model3(a_dim, b_dim, N, M, phi_in, theta_in, phi_out, theta_out, lambda);

% Define reflected angles based on receiver position and reflector center
phi_rx = round(asind((X - reflector_center(1)) ./ distance_ref_rx), 1);
theta_rx = round(asind((reflector_center(2) - Y) ./ distance_ref_rx), 1);

a1 = abs(round((phi_out_max + phi_rx) ./ pattern_resolution));
a2 = abs(round((theta_out_max + theta_rx) ./ pattern_resolution));
a1(isnan(a1)) = 0;
a2(isnan(a2)) = 0;

% Initialize LOS_grid based on receiver position
for f = 1:length(receiver_y)
    for g = 1:length(receiver_x)
        if receiver_y(f) <= 190 && receiver_y(f) >= 170
            LOS_grid(f, g) = 1;
        elseif (receiver_y(f) >= 190 && (receiver_x(g) < 170 || receiver_x(g) > 190)) || (receiver_y(f) <= 170 && (receiver_x(g) <= 170 || receiver_x(g) > 190))
            LOS_grid(f, g) = 2;
        end
    end
end

% Calculate Path Loss and Power Received for different cases
[PL_LOS_dB, PL_NLOS_dB, PL_Radar_dB] = Fn_FSPL(freq.' * 1e-9, distance_tx_ref, distance_ref_rx, distance_tx_rx, c, distance_tx_ref);

% Initialize arrays for power and gain calculations
P_RX1 = zeros(size(X));
P_RX2 = zeros(size(X));
PL1 = zeros(size(X));
PL2 = zeros(size(X));
gain_rcs1 = zeros(size(X));
gain_rcs2 = zeros(size(X));
delta_power = zeros(size(X));

% Compute power and gain for LOS and NLOS cases
for f = 1:length(receiver_y)
    for g = 1:length(receiver_x)
        if (receiver_y(f) <= 190 && receiver_y(f) >= 170) || (receiver_x(g) < 190 && receiver_x(g) > 170)
            if LOS_grid(f, g) == 1  % LOS case
                PL1(f, g) = PL_LOS_dB(f, g);
                P_RX1(f, g) = P_TX + G_TX + G_RX - PL1(f, g);
            else
                PL1(f, g) = -Inf; 
                P_RX1(f, g) = -Inf;
                gain_rcs1(f, g) =  -Inf;
            end
        else
            gain_rcs1(f, g) = NaN;
            PL1(f, g) = NaN; 
            P_RX1(f, g) = NaN;
        end
    end
end

for f = 1:length(receiver_y)
    for g = 1:length(receiver_x)
        if (receiver_y(f) <= 190 && receiver_y(f) >= 170) || (receiver_x(g) < 190 && receiver_x(g) > 170)
            if (receiver_y(f) <= 190 && receiver_y(f) >= 170) || (receiver_x(g) < 190 && receiver_x(g) > 170 && receiver_y(f) <= 190)  % LOS
                PL2(f, g) = PL_Radar_dB(f, g);
                gain_rcs2(f, g) = sigma_sum_mat(a1(f, g) + 1, a2(f, g) + 1);
                P_RX2(f, g) = P_TX + G_TX + G_RX + gain_rcs2(f, g) - PL2(f, g);
            else  
                gain_rcs2(f, g) = -Inf;
                PL2(f, g) = -Inf; 
                P_RX2(f, g) = -Inf;
            end
        else
            gain_rcs2(f, g) = NaN;
            PL2(f, g) = NaN; 
            P_RX2(f, g) = NaN;
        end
    end
end

% Combine power from both cases
Combined_P_RX = max(P_RX1, P_RX2);

% Check conditions and adjust power values
for f = 1:length(receiver_y)
    for g = 1:length(receiver_x)
        if (receiver_y(f) <= 190 && receiver_y(f) >= 170) || (receiver_x(g) < 190 && receiver_x(g) > 170)
            if (receiver_y(f) <= 190 && receiver_y(f) >= 170) || (receiver_x(g) < 190 && receiver_x(g) > 170 && receiver_y(f) <= 190)  % LOS
                if LOS_grid(f, g) == 1
                    delta_power(f, g) = P_RX2(f, g) - P_RX1(f, g);
                elseif LOS_grid(f, g) == 0
                    delta_power(f, g) = P_RX2(f, g) + 100;
                end
                if distance_ref_rx(f, g) <= FF_dist
                    delta_power(f, g) = -Inf;
                    gain_rcs2(f, g) = -Inf;
                    Combined_P_RX(f, g) = -Inf;
                end
            else
                delta_power(f, g) = -Inf;
            end
        else
            delta_power(f, g) = NaN;
        end
        % Check condition and adjust power values based on reflector center
        if reflector_center(1) == 195 && NLOS_region_case_sens(f, g) == 1
            delta_power(f, g) = -Inf;
            gain_rcs2(f, g) = -Inf; 
            Combined_P_RX(f, g) = -Inf; 
        end
    end
end

% Plot results
figure;
subplot(1,2,1);
plot_heatmap_analytical(gain_rcs2, receiver_x, receiver_y, ...
              'Gain [dB]', 'X-position [m]', 'Y-position [m]', ...
              -30, 60) 
subplot(1,2,2);
plot_heatmap_analytical(Combined_P_RX, receiver_x, receiver_y, ...
              'Maximum Power Received [dBm]', 'X-position [m]', 'Y-position [m]', ...
              -130, -10) 
