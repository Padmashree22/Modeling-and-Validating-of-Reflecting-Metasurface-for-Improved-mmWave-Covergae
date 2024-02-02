% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 


function [tx_gamma] = gamma_tx(theta_k_in, phi_k_in)
    % gamma_tx - Calculate the transmission coefficient alpha for given angles
    %
    %   [tx_gamma] = gamma_tx(theta_k_in, phi_k_in) calculates the
    %   transmission coefficient alpha for a set of incident angles
    %   (theta_k_in) and azimuthal angles (phi_k_in) in a wireless
    %   communication scenario.
    %
    % Input:
    %   theta_k_in - Vector of incident angles in degrees
    %   phi_k_in   - Vector of azimuthal angles in degrees
    %
    % Output:
    %   tx_gamma   - Matrix containing the calculated transmission
    %                coefficients for each combination of incident and
    %                azimuthal angles
    
    % Constants
    f = 28e9; % Frequency of operation in Hz
    c = 3e8;  % Speed of light in m/s
    
    k = (2 * pi * f) / c; % Wavenumber
    
    % Initialize variables
    y = 1;
    for theta_k = theta_k_in
        x = 1;
        for phi_k = phi_k_in
            % Calculate components of the wavenumber vector
            kx = k * sind(theta_k) * cosd(phi_k);
            ky = k * sind(theta_k) * sind(phi_k);
            kz = k * cosd(theta_k);
            
            % Calculate the magnitude of the wavenumber vector
            K = sqrt(kx^2 + ky^2 + kz^2);
            
            % Calculate the transmission angles in spherical coordinates
            tx_theta(x, y) = acosd(kz / K);
            tx_phi(x, y) = atand(ky / kx);
            
            % Calculate the transmission coefficient alpha
            tx_gamma(x, y) = sqrt(tx_theta(x, y)^2 + tx_phi(x, y)^2);
            
            x = x + 1;
        end
        y = y + 1;
    end
end
