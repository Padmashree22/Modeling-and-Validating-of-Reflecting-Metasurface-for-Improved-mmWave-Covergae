% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [PL_UMa_LOS, PL_UMa_NLOS, PL_UMa_Radar] = UMa_scenario(freq, D_SR, D_RD, D_SD, h_BS_UMa, h_Ref, h_UT, h_E, c, Ref_location)
    % UMa_scenario - UMa (Urban Macro) scenario path loss model
    %
    %   [PL_UMa_LOS, PL_UMa_NLOS, PL_UMa_Radar] = UMa_scenario(freq, D_SR, D_RD, D_SD, h_BS_UMa, h_Ref, h_UT, h_E, c, Ref_location)
    %   calculates the path loss components for the UMa scenario based on
    %   the given parameters.
    %
    % Input:
    %   freq - Frequency of operation (in Hertz)
    %   D_SR - Distance from Base Station (BS) to Reflector (in meters)
    %   D_RD - Distance from Reflector to User Terminal (UT) (in meters)
    %   D_SD - Distance from BS to UT (in meters)
    %   h_BS_UMa - Height of the BS antenna (in meters)
    %   h_Ref - Height of the Reflector (in meters)
    %   h_UT - Height of the UT antenna (in meters)
    %   h_E - Height of the environment (in meters)
    %   c - Speed of light (in meters per second)
    %   Ref_location - Location of the Reflector (in meters)
    %
    % Output:
    %   PL_UMa_LOS - Path loss for Line-of-Sight (LOS) scenario from Tx-Rx
    %   PL_UMa_NLOS - Path loss for Non-Line-of-Sight (NLOS) scenario from Tx-Rx
    %   PL_UMa_Radar - Path loss for Radar scenario based on LOS or NLOS conditions
    %
    % Note: This function relies on the helper functions UMa_LOS and UMa_NLOS.

    % LOS calculations from TX-Reflector and Reflector-Rx
    [PL_UMa_LOS1, ~, ~] = UMa_LOS(freq, h_BS_UMa, h_Ref, D_SR, h_E);
    [PL_UMa_LOS2, ~, ~] = UMa_LOS(freq, h_Ref, h_UT, D_RD, h_E);
    PL_UMa_LOS_Total = PL_UMa_LOS1 + PL_UMa_LOS2;

    % LOS calculations from Tx-Rx and NLOS calculations from Tx-Rx
    [PL_UMa_LOS, d_3D_LOS] = UMa_LOS(freq, h_BS_UMa, h_UT, D_SD, h_E);
    PL_UMa_NLOS = UMa_NLOS(d_3D_LOS, freq, h_UT, PL_UMa_LOS);

    % Radar scenario calculations
    [a, b] = size(D_SD);
    delta_UMa = 57.15 - 20 * log10(c) + 20 * log10(freq);
    PL_UMa_Radar1 = PL_UMa_LOS_Total + delta_UMa;

    % Radar path loss based on LOS or NLOS conditions
    PL_UMa_Radar = zeros(a, b);
    for a_idx = 1:a
        for b_idx = 1:b
            if D_SD(a_idx, b_idx) <= Ref_location
                PL_UMa_Radar(a_idx, b_idx) = PL_UMa_LOS(a_idx, b_idx); % LOS from Tx-Rx
            else
                PL_UMa_Radar(a_idx, b_idx) = PL_UMa_Radar1(a_idx, b_idx);
            end
        end
    end
end
