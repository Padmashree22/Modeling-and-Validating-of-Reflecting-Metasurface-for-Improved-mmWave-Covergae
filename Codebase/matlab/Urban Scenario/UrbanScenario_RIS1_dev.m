% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear command window, close all figures, and clear workspace
clc;
close all;
clear all;

% Add path to helper function
addpath('Helper_function')

% Load necessary data
load('LOS_scenario1_grid','LOS_grid')
load('NLOS_region')

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

% Transmit power, transmit antenna gain, receive antenna gain, receiver sensitivity
P_TX = 20.3; % dBm
G_TX = 19.7; % dBi
G_RX = 13.5; % dBi
R_Sensitivity = 83.5; % dBm

% Grid for receivers
receiver_x = 120:1:240;
receiver_y = 120:1:240;
receiver_z = 1.5;
[X, Y] = meshgrid(receiver_x, receiver_y);
Y = rot90(Y, 2);

% Incident angles
phi_in = 80.17; % Incident azimuth angle [deg]
theta_in = -2.46; % Incident elevation angle [deg]

%% Pattern resolution for RCS calculation

pattern_resolution = 0.1; % Resolution of RCS [deg]

% Define ranges for reflected azimuth and elevation waves
phi_out_min = -90; % Minimum reflected azimuth wave [deg] 
phi_out_max = 90; % Maximum reflected azimuth wave [deg]
theta_out_min = -90; % Minimum reflected elevation wave [deg]
theta_out_max = 90; % Maximum reflected elevation wave [deg]

%% Module parameters

% Number of modules along Y axis and Z axis
N = 1:32;
M = 1;

% Dimensions along Y and Z axis [m]
a_dim = 3 / 32;
b_dim = 3;

% Ranges for reflected elevation and azimuth waves
theta_out = theta_out_min:pattern_resolution:theta_out_max;
phi_out = phi_out_min:pattern_resolution:phi_out_max;

% Wavelength
lambda = c / freq; % Wavelength [m]

%% Distances

% Tx-Rx distance
distance_tx_rx = sqrt((X - tx_position(1)).^2 + (Y - tx_position(2)).^2 + (receiver_z - tx_position(3)).^2);

% Tx-Reflector distance
distance_tx_ref = sqrt((reflector_center(1) - tx_position(1)).^2 + (reflector_center(2) - tx_position(2)).^2 + (reflector_center(3) - tx_position(3)).^2);

% Rx-Reflector distance
distance_ref_rx = sqrt((X - reflector_center(1)).^2 + (Y - reflector_center(2)).^2 + (receiver_z - reflector_center(3)).^2);

% Free space path loss
FF_dist = max(((a_dim / 10).^2), ((b_dim / 10).^2)) / (lambda);

% Calculate the Radar Cross Section (RCS) using a helper function for
% narrow beam case of IRS model 1
theta_rx1=-25.75;
phi_rx1=-2.5;
[sigma_sum_mat]= IRS_model1_dev(a_dim,b_dim,N,M,phi_in,theta_in,phi_out,theta_out,theta_rx1,phi_rx1,lambda);

% Calculate angles for reflected waves
phi_rx = round(asind((X - reflector_center(1)) ./ distance_ref_rx), 1);
theta_rx = round(asind((reflector_center(2) - Y) ./ distance_ref_rx), 1);

% Initialize variables
a1 = abs(round((phi_out_max + phi_rx) ./ pattern_resolution));
a2 = abs(round((theta_out_max + theta_rx) ./ pattern_resolution));
a1(isnan(a1)) = 0;
a2(isnan(a2)) = 0;

% Initialize LOS grid based on receiver position
LOS_grid = zeros(size(X));
for f = 1:length(receiver_y)
    for g = 1:length(receiver_x)
        if receiver_y(f) <= 190 && receiver_y(f) >= 170
            LOS_grid(f, g) = 1;
        elseif (receiver_y(f) >= 190 && (receiver_x(g) < 170 || receiver_x(g) > 190)) || (receiver_y(f) <= 170 && (receiver_x(g) <= 170 || receiver_x(g) > 190))
            LOS_grid(f, g) = 2;
        end
    end
end

% Calculate Free Space Path Loss for LOS, NLOS, and Radar scenarios
[PL_LOS_dB, PL_NLOS_dB, PL_Radar_dB] = Fn_FSPL(freq.' * 1e-9, (distance_tx_ref), (distance_ref_rx), (distance_tx_rx), c, distance_tx_ref);

% Initialize variables for LOS and Radar scenarios
PL1 = zeros(size(X));
P_RX1 = zeros(size(X));
gain_rcs1 = zeros(size(X));

PL2 = zeros(size(X));
P_RX2 = zeros(size(X));
gain_rcs2 = zeros(size(X));

% Calculate received power for LOS and Radar scenarios
for f = 1:length(receiver_y)
    for g = 1:length(receiver_x)
        if (receiver_y(f) <= 190 && receiver_y(f) >= 170) || (receiver_x(g) < 190 && receiver_x(g) > 170)
            if LOS_grid(f, g) == 1
                % LOS scenario
                PL1(f, g) = PL_LOS_dB(f, g);
                P_RX1(f, g) = P_TX + G_TX + G_RX - PL1(f, g);
            else
                % NLOS scenario
                PL1(f, g) = -Inf;
                P_RX1(f, g) = -Inf;
                gain_rcs1(f, g) = -Inf;
            end
        else
            gain_rcs1(f, g) = NaN;
            PL1(f, g) = NaN;
            P_RX1(f, g) = NaN;
        end
    end
end

% Calculate received power for Radar scenario
for f = 1:length(receiver_y)
    for g = 1:length(receiver_x)
        if (receiver_y(f) <= 190 && receiver_y(f) >= 170) || (receiver_x(g) < 190 && receiver_x(g) > 170)
            if (receiver_y(f) <= 190 && receiver_y(f) >= 170) || (receiver_x(g) < 190 && receiver_x(g) > 170 && receiver_y(f) <= 190)
                % LOS scenario
                PL2(f, g) = PL_Radar_dB(f, g);
                gain_rcs2(f, g) = sigma_sum_mat(a1(f, g) + 1, a2(f, g) + 1);
                P_RX2(f, g) = P_TX + G_TX + G_RX + gain_rcs2(f, g) - PL2(f, g);
            else
                % NLOS scenario
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

% Combine received powers from LOS and Radar scenarios
Combined_P_RX = max(P_RX1, P_RX2);

% Calculate delta power and handle special cases
delta_power = zeros(size(X));
for f = 1:length(receiver_y)
    for g = 1:length(receiver_x)
        if (receiver_y(f) <= 190 && receiver_y(f) >= 170) || (receiver_x(g) < 190 && receiver_x(g) > 170)
            if (receiver_y(f) <= 190 && receiver_y(f) >= 170) || (receiver_x(g) < 190 && receiver_x(g) > 170 && receiver_y(f) <= 190)
                % LOS scenario
                if LOS_grid(f, g) == 1
                    delta_power(f, g) = P_RX2(f, g) - P_RX1(f, g);
                elseif LOS_grid(f, g) == 0
                    delta_power(f, g) = P_RX2(f, g) + 100;
                end
                % Check for near field condition
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
        % Check for specific condition related to reflector center
        if reflector_center(1) == 195 && NLOS_region_case_sens(f, g) == 1
            delta_power(f, g) = -Inf;
            gain_rcs2(f, g) = -Inf;
            Combined_P_RX(f, g) = -Inf;
        end
    end
end

% Plotting
figure;
subplot(1, 2, 1);
plot_heatmap_analytical(gain_rcs2, receiver_x, receiver_y, ...
    'Gain [dB]', 'X-position [m]', 'Y-position [m]', ...
    -30, 60)

subplot(1, 2, 2);
plot_heatmap_analytical(Combined_P_RX, receiver_x, receiver_y, ...
    'Maximum Power Received [dBm]', 'X-position [m]', 'Y-position [m]', ...
    -130, -10)