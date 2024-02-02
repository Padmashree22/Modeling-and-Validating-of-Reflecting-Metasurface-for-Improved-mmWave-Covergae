% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [PL_NLOS] = UMi_NLOS(d_3D, freq, h_UT, PL_LOS)
    % UMi_NLOS - Urban Macrocell NLOS path loss calculation
    %
    %   PL_NLOS = UMi_NLOS(d_3D, freq, h_UT, PL_LOS) calculates the
    %   non-line-of-sight (NLOS) path loss for an Urban Macrocell (UMi)
    %   environment based on the given input parameters.
    %
    % Input:
    %   d_3D - 3D distance between the transmitter and receiver
    %   freq - Frequency of the signal in Hertz
    %   h_UT - Height of the user terminal (receiver) in meters
    %   PL_LOS - Line-of-sight (LOS) path loss
    %
    % Output:
    %   PL_NLOS - Non-line-of-sight (NLOS) path loss
    %
    % Reference:
    %   The formula is based on the reference:
    %   https://www.etsi.org/deliver/etsi_tr/138900_138999/138901/16.01.00_60/tr_138901v160100p.pdf
    %
    
    % Calculate the NLOS path loss using the UMi NLOS model
    PL_q_NLOS = 22.4 + 35.3 * log10(d_3D) + 21.3 * log10(freq) - 0.3 * (h_UT - 1.5);
    
    % Get the dimensions of d_3D
    [a, b] = size(d_3D);
    
    % Initialize PL_NLOS
    PL_NLOS = zeros(a, b);
    
    % Iterate through the elements of d_3D to calculate PL_NLOS
    for x = 1:a
        for y = 1:b
            % Calculate PL_NLOS for each element
            PL_NLOS(x, y) = max(PL_LOS(x, y), PL_q_NLOS(x, y));
        end
    end
end
