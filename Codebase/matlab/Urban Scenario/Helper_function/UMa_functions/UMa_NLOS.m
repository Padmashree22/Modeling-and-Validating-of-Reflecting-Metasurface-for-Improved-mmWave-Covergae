% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [PL_NLOS] = UMa_NLOS(d_3D, freq, h_UT, PL_LOS)
    % UMa_NLOS - Calculate Non-Line-of-Sight (NLOS) Path Loss in UMa model
    %
    %   [PL_NLOS] = UMa_NLOS(d_3D, freq, h_UT, PL_LOS) calculates the
    %   Non-Line-of-Sight (NLOS) path loss using the UMa (Urban Macro)
    %   model based on the specified inputs.
    %
    % Inputs:
    %   d_3D  - 3D distance between transmitter and receiver (in meters)
    %   freq  - Frequency of the communication (in Hertz)
    %   h_UT  - Height of the user terminal (in meters)
    %   PL_LOS - Path Loss in Line-of-Sight (LOS) condition
    %
    % Output:
    %   PL_NLOS - Computed Non-Line-of-Sight (NLOS) path loss
    %
    % Reference:
    %   The calculation is based on the UMa model as specified in:
    %   https://www.etsi.org/deliver/etsi_tr/138900_138999/138901/16.01.00_60/tr_138901v160100p.pdf
    %

    % Calculate the NLOS path loss using the UMa model formula
    PL_q_NLOS = 13.54 + 39.08 * log10(d_3D) + 20 * log10(freq) - 0.6 * (h_UT - 1.5);

    % Get the size of the distance matrix
    [a, b] = size(d_3D);

    % Initialize PL_NLOS matrix
    PL_NLOS = zeros(a, b);

    % Loop through each element of the distance matrix
    for x = 1:a
        for y = 1:b
            % Calculate NLOS path loss as the maximum of PL_LOS and PL_q_NLOS
            PL_NLOS(x, y) = max(PL_LOS(x, y), PL_q_NLOS(x, y));
        end
    end
end
