% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

% Clear the command window, close all figures, and clear all variables
clc; 
close all; 
clear all;

% Add the path to the helper functions
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
reflector_center = [195, 170, 5];  
reflector_size = [3, 3];  % [length, width] in meters

% Effective environment height
h_E=1;          

% Transmitter power, transmitter and receiver gains, and receiver sensitivity
P_TX= 20.3; % dBm
G_TX= 19.7; % dBi
G_RX= 13.5; % dBi
R_Sensitivity=83.5; %dBm

% Grid for the receiver
receiver_x= 120:1:240;
receiver_y= 120:1:240;
receiver_z=1.5;
[X, Y] = meshgrid(receiver_x, receiver_y);
Y=rot90(Y,2);

% Incident azimuth and elevation angles
phi_in = 80.17; % Incident azimuth angle [deg]
theta_in = -2.46; % Incident elevation angle [deg]

%% Pattern Resolution and Variation Parameters

% Resolution of RCS [deg]
pattern_resolution = 0.1; 

% Ranges for reflected azimuth and elevation waves [deg]
phi_out_min=-90;
phi_out_max=90;
theta_out_min=-90;
theta_out_max=90;

% Module and dimension parameters
N=1:32; % Modules along Y axis 
M=1; % Modules along Z axis
a_dim=3/32; % Dimensions along Y axis [m]
b_dim=3; % Dimensions along Z axis [m]

% Slope angle variation
alpha_var=17.6;
beta_var=0;
var=2.5;

%% (Uncomment the following lines if you are calculating narrow beam)
% alpha=(alpha_var).*ones(length(M),length(N)); % Hor. Slope angle [deg]
% beta= (beta_var).*ones(length(M),length(N)); % Vert. Slope angle [deg]

%% Broad beam calculations
alpha=[(alpha_var)*(ones(1,2)) (alpha_var+var)*(ones(1,4)) (alpha_var+2*var)*(ones(1,5)) (alpha_var+3*var)*(ones(1,6)) (alpha_var+4*var)*(ones(1,6)) (alpha_var+5*var)*(ones(1,5)) (alpha_var+6*var)*(ones(1,4))];
alpha=flip(alpha);

beta=[(beta_var)*(ones(1,2)) (beta_var+var)*(ones(1,4)) (beta_var+2*var)*(ones(1,5)) (beta_var+3*var)*(ones(1,6)) (beta_var+4*var)*(ones(1,6)) (beta_var+5*var)*(ones(1,5)) (beta_var+6*var)*(ones(1,4))];
beta=flip(beta);

%% Reflected elevation and azimuth wave angles
theta_out = theta_out_min:pattern_resolution:theta_out_max; % Reflected elevation wave [deg]
phi_out = phi_out_min:pattern_resolution:phi_out_max; % Reflected azimuth wave [deg]
lambda=c/freq; % wavelength [m]

%% Distances
% Compute distances
distance_tx_rx=(sqrt((X-tx_position(1)).^2+(Y-tx_position(2)).^2+(receiver_z-tx_position(3)).^2)); % Tx-Rx distance
distance_tx_ref=(sqrt((reflector_center(1)-tx_position(1)).^2+(reflector_center(2)-tx_position(2)).^2+(reflector_center(3)-tx_position(3)).^2)); % Tx-Reflector distance
distance_ref_rx=(sqrt((X-reflector_center(1)).^2+(Y-reflector_center(2)).^2+(receiver_z-reflector_center(3)).^2)); % Rx-Reflector distance

FF_dist= (max(((a_dim/9).^2),((b_dim/9).^2)))./(lambda);

%% Gain_RCS
[sigma_sum_max, sigma_sum_mat, eta_a, eta_b] = fn_SF_HELIOS_array_eff(a_dim,b_dim,N,M,alpha,beta,phi_in,theta_in,phi_out,theta_out,lambda);
sigma_sum_mat=(10.*log10(sigma_sum_mat.^2));

%% Radar Cross-Section Calculation using HELIOS Array
% (Uncomment the following lines if you want to use different IRS models)

% [sigma_sum_mat] = IRS_model1(a_dim,b_dim,N,M,phi_in,theta_in,phi_out,theta_out,lambda);
% [sigma_sum_mat] = IRS_model2(a_dim,b_dim,N,M,phi_in,theta_in,phi_out,theta_out,lambda);
% [sigma_sum_mat] = IRS_model3(a_dim,b_dim,N,M,phi_in,theta_in,phi_out,theta_out,lambda);

%     theta_rx=-25.75;
%     phi_rx=-2.5;
% [sigma_sum_mat]= IRS_model1_dev(a_dim,b_dim,N,M,phi_in,theta_in,phi_out,theta_out,theta_rx,phi_rx,lambda);

%% Calculate Incident Angles at the Receiver
phi_rx= round(asind((X-reflector_center(1))./distance_ref_rx),1);
theta_rx= round(asind((reflector_center(2)-Y)./distance_ref_rx),1);

a1=abs(round((phi_out_max+phi_rx)./pattern_resolution));
a2=abs(round((theta_out_max+theta_rx)./pattern_resolution));
a1(isnan(a1))=0;
a2(isnan(a2))=0;

% Update LOS_grid based on receiver positions
for f=1:length(receiver_y)
    for g=1:length(receiver_x)
        if receiver_y(f)<=190 && receiver_y(f)>=170
            LOS_grid(f,g)=1;
        elseif (receiver_y(f)>=190 && (receiver_x(g)<170 || receiver_x(g)>190)) || (receiver_y(f)<=170 && (receiver_x(g)<=170 || receiver_x(g)>190))
            LOS_grid(f,g)=2;
        end
    end
end

%% Path Loss Calculation
[PL_LOS_dB, PL_NLOS_dB, PL_Radar_dB] = UMi_scenario(freq.'*1e-9,(distance_tx_ref), (distance_ref_rx),(distance_tx_rx), tx_position(3), reflector_center(3), receiver_z, h_E, c, distance_tx_ref);

% Initialize arrays for path loss and received power
PL1 = zeros(size(X));
PL2 = zeros(size(X));
P_RX1 = zeros(size(X));
P_RX2 = zeros(size(X));

% Update PL1 and P_RX1 based on LOS_grid
for f=1:length(receiver_y)
    for g=1:length(receiver_x)
        if (receiver_y(f)<=190 && receiver_y(f)>=170) || (receiver_x(g)<190 && receiver_x(g)>170) 
            if LOS_grid(f,g)==1      %% LOS
                PL1(f,g)= PL_LOS_dB(f,g);
            elseif LOS_grid(f,g)==0
                PL1(f,g)= PL_NLOS_dB(f,g); 
            else
                PL1(f,g)= -Inf; 
            end
        else
            PL1(f,g)= NaN; 
        end
        P_RX1(f,g)= P_TX + G_TX + G_RX - PL1(f,g) ;
    end
end

% Update PL2 and P_RX2 based on conditions
for f=1:length(receiver_y)
    for g=1:length(receiver_x)
       if (receiver_y(f)<=190 && receiver_y(f)>=170) || (receiver_x(g)<190 && receiver_x(g)>170) 
         if (receiver_x(g)<190 && receiver_x(g)>170 && receiver_y(f)<170)        %% LOS
            PL2(f,g)= max(PL_Radar_dB(f,g),PL_NLOS_dB(f,g));
            gain_rcs2(f,g)=sigma_sum_mat(a1(f,g)+1,a2(f,g)+1);
         else  
            gain_rcs2(f,g)= -Inf;
            PL2(f,g)=-Inf; 
         end
        else
            gain_rcs2(f,g)= NaN;
             PL2(f,g)= NaN; 
       end
       P_RX2(f,g)= P_TX + G_TX + G_RX + gain_rcs2(f,g) - PL2(f,g) ; 
    end
end

% Combine the received power from both paths
Combined_P_RX= max(P_RX1,P_RX2);

% Apply additional conditions (near field region) and uncovered area due to reflector placement at 195 and update gain and power arrays
for f=1:length(receiver_y)
    for g=1:length(receiver_x)
       if (receiver_y(f)<=190 && receiver_y(f)>=170) || (receiver_x(g)<190 && receiver_x(g)>170) 
         if (receiver_y(f)<=190 && receiver_y(f)>=170) || (receiver_x(g)<190 && receiver_x(g)>170 && receiver_y(f)<=190) 
            if distance_ref_rx(f,g) <= FF_dist
                gain_rcs2(f,g)=-Inf;
                Combined_P_RX(f,g)=-Inf;
            end
         end
        if reflector_center(1)==195 && NLOS_region_case_sens(f,g)==1
            gain_rcs2(f,g)=-Inf; 
            Combined_P_RX(f,g)=-Inf; 
        end
       end
    end
end

% Plot the results
figure;
subplot(1,2,1);
plot_heatmap_analytical(gain_rcs2, receiver_x, receiver_y, ...
             'Gain [dB]', 'X-position [m]', 'Y-position [m]', ...
             -30, 60) 
subplot(1,2,2);
plot_heatmap_analytical(Combined_P_RX, receiver_x, receiver_y, ...
             'Maximum Power Received [dBm]', 'X-position [m]', 'Y-position [m]', ...
             -130,-10)

