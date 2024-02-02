% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [PL_LOS, d_3D, PL_case] = UMa_LOS(freq, h_BS, h_UT, d_2D_in, h_E)
    % UMa_LOS - UMa (Urban Macrocell) LOS (Line-of-Sight) path loss model
    % Calculates the path loss for the LOS scenario based on the UMa model.
    %
    % Input:
    %   freq - Frequency of operation in GHz
    %   h_BS - Height of the base station antenna in meters
    %   h_UT - Height of the user terminal antenna in meters
    %   d_2D_in - 2D distance between BS and UT in meters (matrix)
    %   h_E - Height of the environment (ground) in meters
    %
    % Output:
    %   PL_LOS - Calculated path loss for the LOS scenario
    %   d_3D - 3D distance between BS and UT in meters
    %   PL_case - Path loss case (1 for Pathloss 1, 2 for Pathloss 2)
    %
    % Reference: https://www.etsi.org/deliver/etsi_tr/138900_138999/138901/16.01.00_60/tr_138901v160100p.pdf
    %

    h_d_BS = h_BS - h_E;
    c = physconst('lightspeed');  % Propagation velocity in free space m/s
    h_d_UT = h_UT - h_E;
    d_d_BP = (4 * h_d_BS * h_d_UT * freq * 10^9) / c;   % Breakpoint distance in meters

    [i, j] = size(d_2D_in);
    PL_case = zeros(i, j);
    PL_LOS = zeros(i, j);
    d_3D = zeros(i, j);

    for x = 1:i
        for y = 1:j
            d_2D = d_2D_in(x, y);

            % Calculate 3D distance
            d_3D(x, y) = sqrt((h_BS - h_UT).^2 + (d_2D).^2);

            % Calculate Pathloss 1
            PL_1(x, y) = 28 + 22 * log10(d_3D(x, y)) + 20 * log10(freq);

            % Calculate Pathloss 2
            PL_2(x, y) = 28 + 40 * log10(d_3D(x, y)) + 20 * log10(freq) - 9 * log10((d_d_BP).^2 + (h_BS - h_UT).^2);

            % Determine Path Loss based on conditions
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
