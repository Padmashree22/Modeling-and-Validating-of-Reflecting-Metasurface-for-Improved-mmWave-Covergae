% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [PL_LOS, d_3D, PL_case] = UMi_LOS(freq, h_BS, h_UT, d_2D_in, h_E)
    % UMi_LOS - UMi (Urban Macro) LOS (Line-of-Sight) path loss model
    %
    % [PL_LOS, d_3D, PL_case] = UMi_LOS(freq, h_BS, h_UT, d_2D_in, h_E) 
    % calculates the path loss for Line-of-Sight (LOS) conditions in a UMi
    % environment based on the given parameters.
    %
    % Input:
    %   freq - Frequency of operation in GHz
    %   h_BS - Height of the base station antenna in meters
    %   h_UT - Height of the user terminal antenna in meters
    %   d_2D_in - 2D distance matrix between the base station and user terminal in meters
    %   h_E - Height of the environment (Earth) in meters
    %
    % Output:
    %   PL_LOS - Path loss for Line-of-Sight conditions
    %   d_3D - 3D distance matrix between the base station and user terminal in meters
    %   PL_case - Path loss case (1 or 2) indicating the LOS path loss model used
    %
    % Reference: https://www.etsi.org/deliver/etsi_tr/138900_138999/138901/16.01.00_60/tr_138901v160100p.pdf
    
    % Calculate the height difference between the base station and environment
    h_d_BS = h_BS - h_E;
    
    % Propagation velocity in free space (m/s)
    c = physconst('lightspeed');
    
    % Calculate the height difference between the user terminal and environment
    h_d_UT = h_UT - h_E;
    
    % Breakpoint distance in meters
    d_d_BP = (4 * h_d_BS * h_d_UT * freq * 1e9) / c;

    % Initialize matrices to store results
    [i, j] = size(d_2D_in);
    PL_case = zeros(i, j);
    PL_LOS = zeros(i, j);
    d_3D = zeros(i, j);

    % Loop through each element of the 2D distance matrix
    for x = 1:i
        for y = 1:j
            % Extract 2D distance
            d_2D = d_2D_in(x, y);
            
            % Calculate 3D distance between Tx and Rx
            d_3D(x, y) = sqrt((h_BS - h_UT)^2 + d_2D^2);

            % Pathloss model 1
            PL_1(x, y) = 32.4 + 21 * log10(d_3D(x, y)) + 20 * log10(freq);

            % Pathloss model 2
            PL_2(x, y) = 32.4 + 40 * log10(d_3D(x, y)) + 20 * log10(freq) - 9.5 * log10(d_d_BP^2 + (h_BS - h_UT)^2);

            % Determine the LOS path loss model based on the 2D distance
            if d_2D <= d_d_BP && d_2D >= 10
                PL_LOS(x, y) = PL_1(x, y);
                PL_case(x, y) = 1;
            end

            if d_2D <= 5000 && d_2D >= d_d_BP
                PL_LOS(x, y) = PL_2(x, y);
                PL_case(x, y) = 2;
            end
        end
    end
end
