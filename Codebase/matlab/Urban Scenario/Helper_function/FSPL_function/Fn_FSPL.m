% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [PL_LOS_dB, PL_Radar_dB] = Fn_FSPL(freq, D_SR, D_RD, D_SD, c, Ref_location)
    % Fn_FSPL - Free Space Path Loss calculation
    %
    %   [PL_LOS_dB, PL_Radar_dB] = Fn_FSPL(freq, D_SR, D_RD, D_SD, c, Ref_location)
    %   calculates the Free Space Path Loss for Line-of-Sight (LOS) and Radar
    %   scenarios based on the input parameters.
    %
    % Input:
    %   freq - Frequency (in Hz)
    %   D_SR - Distance from the radar to the source (in meters)
    %   D_RD - Distance from the radar to the destination (in meters)
    %   D_SD - Distance from the source to the destination (in meters)
    %   c - Speed of light (in meters per second)
    %   Ref_location - Reference location (in meters)
    %
    % Output:
    %   PL_LOS_dB - Free Space Path Loss for Line-of-Sight (in dB)
    %   PL_Radar_dB - Free Space Path Loss for Radar scenario (in dB)
    %

    % Calculate Free Space Path Loss for Line-of-Sight (LOS) from Tx to Rx
    PL_LOS_dB = -147.5 + 20 * log10(D_SD) + 20 * log10(freq.' * 1e9);

    % Calculate Free Space Path Loss for Radar scenario
    PL_Radar_dB1 = -136.56 + 20 * log10(D_SR) + 20 * log10(D_RD) + 20 * log10(freq.' * 1e9);

    % Determine the size of D_SD
    [i, j] = size(D_SD);

    % Initialize PL_Radar_dB matrix
    PL_Radar_dB = zeros(i, j);

    % Iterate through each element of D_SD
    for a = 1:i
        for b = 1:j
            % Check if D_SD is less than or equal to Ref_location
            if D_SD(a, b) <= Ref_location
                % Calculate Free Space Path Loss for Radar scenario with LOS
                PL_Radar_dB(a, b) = -147.5 + 20 * log10(D_SD(a, b)) + 20 * log10(freq.' * 1e9);
            else
                % Use the pre-calculated PL_Radar_dB1 value
                PL_Radar_dB(a, b) = PL_Radar_dB1(a, b);
            end
        end
    end
end
