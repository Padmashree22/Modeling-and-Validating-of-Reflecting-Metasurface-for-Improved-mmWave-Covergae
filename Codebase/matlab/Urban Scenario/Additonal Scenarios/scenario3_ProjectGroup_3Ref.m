clc; close all;
clear all;

addpath('Helper_function') % load helper_function
%load('NLOS_scenario2_grid','NLOS_grid')

%% Define simulation parameters
freq = 28e9;  % Frequency in Hz (28 GHz)
c = physconst('Light');  % Speed of light
tx_position = [5.5, 8.1, 1.2];  % Transmitter position [x, y, z] in meters 
reflector1_center = [1.5, 8.1, 1.2];  % Reflector footprint [x, y, z] in meters
reflector2_center = [1.5, 8.1, 2.2];  % Reflector footprint [x, y, z] in meters
reflector3_center = [0, 7.1, 1.2];  % Reflector footprint [x, y, z] in meters
reflector_size = [3, 3];  % [length, width] in meters

P_TX= 20.3; % dBm
G_TX= 19.7; % dBi
G_RX= 13.5; % dBi

% Grid receiver
receiver_x= linspace(0, 3, 7);
receiver_y= linspace(0, 5.1, 18);
receiver_z=1.2; 
[X, Y] = meshgrid(receiver_x, receiver_y);
%Y=rot90(Y,2);

phi_out = linspace(-55, -35, 91);
theta_out = -1:1:1;

%%
pattern_resolution = 0.1; % resolution of RCS [deg]
phi_out_min=-10; % Minimum reflected azimuth wave [deg] 
phi_out_max=10; % Maximum reflected azimuth wave [deg]
theta_out_min=-1; % Minimum reflected elevation wave [deg]
theta_out_max=1; % Maximum reflected elevation wave [deg]
% phi_out_min=-90; % Minimum reflected azimuth wave [deg] 
% phi_out_max=90; % Maximum reflected azimuth wave [deg]
% theta_out_min=-90; % Minimum reflected elevation wave [deg]
% theta_out_max=90; % Maximum reflected elevation wave [deg]

%%
N=1:3; % Modules along Y axis 
M=1:3; % Modules along Z axis
a_dim=0.1; % Dimensions along Y axis [m]
b_dim=0.1; % Dimensions along Z axis [m]
alpha=(10).*ones(length(M),length(N)); % Hor. Slope angle [deg]
beta= 0.*ones(length(M),length(N)); % Vert. Slope angle [deg]

%%
theta_out = theta_out_min:pattern_resolution:theta_out_max; % Reflected elevation wave [deg]
phi_out = phi_out_min:pattern_resolution:phi_out_max; % Reflected azimuth wave [deg]
lambda=c/freq; % wavelength [m]

%% Distances
distance_tx_rx=abs(sqrt((X-tx_position(1)).^2+(Y-tx_position(2)).^2+(receiver_z-tx_position(3)).^2)); % Tx-Rx distance
distance_tx_ref=abs(sqrt((reflector_center(1)-tx_position(1)).^2+(reflector_center(2)-tx_position(2)).^2+(reflector_center(3)-tx_position(3)).^2)); % Tx-Reflector distance
distance_ref_rx=abs(sqrt((X-reflector_center(1)).^2+(Y-reflector_center(2)).^2+(receiver_z-reflector_center(3)).^2)); % Rx-Reflector distance
%% Gain_RCS
[sigma_sum_max, sigma_sum_mat, eta_a, eta_b] = fn_SF_HELIOS_array_eff(a_dim,b_dim,N,M,alpha,beta,phi_in,theta_in,phi_out,theta_out,lambda);
sigma_sum_mat=(10*log10(sigma_sum_mat.^2));

phi_rx= round(asind((reflector_center(1)-X)./distance_ref_rx),1);
% theta_rx= round(asind((reflector_center(2)-Y)./distance_ref_rx),1);
theta_rx= zeros(length(receiver_y),length(receiver_x));
%NLOS_grid= ones(length(receiver_y),length(receiver_x));
a1=abs(round((phi_out_max+phi_rx)./pattern_resolution));
 a2=abs(round((theta_out_max+theta_rx)./pattern_resolution));
% a2=zeros(length(receiver_y),length(receiver_x));
a1(isnan(a1))=0;
a2(isnan(a2))=0;

[h,i]=size(sigma_sum_mat);
for f=1:length(receiver_y)
    for g=1:length(receiver_x)
        if NLOS_grid(f,g)==1 && a2(f,g)+1 <= h && a1(f,g)+1 <= i
           % PL(f,g)= -295 + (20*log10(distance_tx_ref)) + (20*log10(distance_ref_rx(f,g))) + 2*(20*log10(freq)); % 2 times FSPL
            PL(f,g)= -136.5 + (20*log10(distance_tx_ref)) + (20*log10(distance_ref_rx(f,g))) + (20*log10(freq)); %Radar Equation
            gain_rcs(f,g)=sigma_sum_mat(a2(f,g)+1,a1(f,g)+1);
            P_RX(f,g)= P_TX + G_TX + G_RX + gain_rcs(f,g) - PL(f,g);
        else
            %  PL(f,g)= -147.5 + (20*log10(distance_tx_rx(f,g))) + (20*log10(freq));
            % P_RX(f,g)= P_TX + G_TX + G_RX - PL(f,g);
            gain_rcs(f,g)=NaN;
            PL(f,g)= NaN; 
            P_RX(f,g)= NaN;
        end
    end
end
 
    subplot(1,3,1);
   plot_heatmap_analytical(gain_rcs, receiver_x, receiver_y, ...
                 'Gain RCS [dB]', 'X-position [m]', 'Y-position [m]', ...
                 min(min(gain_rcs)), max(max(gain_rcs))) 
    subplot(1,3,2);
   plot_heatmap_analytical(PL, receiver_x, receiver_y, ...
                 'PL Radar [dB]', 'X-position [m]', 'Y-position [m]', ...
                  min(min(PL)), max(max(PL))) 
   subplot(1,3,3);
   plot_heatmap_analytical(P_RX, receiver_x, receiver_y, ...
                 'P_{RX} [dB]', 'X-position [m]', 'Y-position [m]', ...
                  min(min(P_RX)), max(max(P_RX))) 