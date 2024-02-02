% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 
function [PL_UMi_LOS, PL_UMi_NLOS, PL_UMi_Radar] = UMi_scenario(freq, D_SR, D_RD, D_SD, h_BS_UMi, h_Ref, h_UT, h_E, c, Ref_location)
    
    % UMi_scenario - UMi (Urban Macrocell) scenario path loss model
    %
    %   [PL_UMi_LOS, PL_UMi_NLOS, PL_UMi_Radar] = UMi_scenario(freq, D_SR, D_RD, D_SD, h_BS_UMi, h_Ref, h_UT, h_E, c, Ref_location)
    %   calculates path loss values for UMi scenario with LOS, NLOS, and Radar conditions.
    %
    % Inputs:
    %   - freq: Frequency of operation (in Hz)
    %   - D_SR: Distance from the Base Station (BS) to Reflector (in meters)
    %   - D_RD: Distance from Reflector to User Terminal (UT) (in meters)
    %   - D_SD: Distance from BS to UT through direct path (in meters)
    %   - h_BS_UMi: Height of the Base Station (BS) in UMi scenario (in meters)
    %   - h_Ref: Height of the Reflector (in meters)
    %   - h_UT: Height of the User Terminal (UT) (in meters)
    %   - h_E: Height of the environment (in meters)
    %   - c: Speed of light (in meters/second)
    %   - Ref_location: Location of the reflector (in meters)
    %
    % Outputs:
    %   - PL_UMi_LOS: Path loss in LOS (Line-of-Sight) condition
    %   - PL_UMi_NLOS: Path loss in NLOS (Non-Line-of-Sight) condition
    %   - PL_UMi_Radar: Path loss in Radar condition
    %
    
    % LOS calculations from TX-Reflector
    [PL_UMi_LOS1, d_3D1, PL_case1] = UMi_LOS(freq, h_BS_UMi, h_Ref, D_SR, h_E);
    
    % LOS calculations from Reflector-Rx
    [PL_UMi_LOS2, d_3D2, PL_case2] = UMi_LOS(freq, h_Ref, h_UT, D_RD, h_E);
    
    % Total LOS path loss
    PL_UMi_LOS_Total = PL_UMi_LOS1 + PL_UMi_LOS2;
    
    % LOS calculations from Tx-Rx
    [PL_UMi_LOS, d_3D_LOS] = UMi_LOS(freq, h_BS_UMi, h_UT, D_SD, h_E);
    
    % NLOS calculations from Tx-Rx
    [PL_UMi_NLOS] = UMi_NLOS(d_3D_LOS, freq, h_UT, PL_UMi_LOS);
    
    % Calculate delta for LOS_Radar
    delta_UMi = 57.15 - 20 * log10(c) + 20 * log10(freq);
    
    % LOS_Radar with delta
    PL_UMi_Radar1 = PL_UMi_LOS_Total + delta_UMi;
    
    [i, j] = size(D_SD);
    
    % Radar path loss calculation
    for a = 1:i
        for b = 1:j
            if D_SD(a, b) <= Ref_location
                % LOS from Tx-Rx
                PL_UMi_Radar(a, b) = PL_UMi_LOS(a, b);
            else
                % Use LOS_Radar1
                PL_UMi_Radar(a, b) = PL_UMi_Radar1(a, b);
            end
        end
    end
end
