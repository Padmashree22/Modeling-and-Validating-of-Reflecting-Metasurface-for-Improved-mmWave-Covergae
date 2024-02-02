% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [rx_gamma] = gamma_rx(theta_r_in, phi_r_in)
    % gamma_rx - Calculate gamma_rx for given theta_r and phi_r
    %
    %   [rx_gamma] = gamma_rx(theta_r_in, phi_r_in) calculates the
    %   gamma_rx values for given theta_r and phi_r. The function uses
    %   parameters such as frequency (f), speed of light (c), and
    %   calculates rx, ry, rz, R, rx_theta, rx_phi, and finally, rx_gamma.
    %
    % Input:
    %   theta_r_in - Input vector of elevation angles (in degrees)
    %   phi_r_in - Input vector of azimuth angles (in degrees)
    %
    % Output:
    %   rx_gamma - Matrix containing calculated gamma_rx values

    % Constants
    f = 28e9;  % Frequency in Hertz
    c = 3e8;   % Speed of light in meters per second
    r = (2 * pi * f) / c;  % Wavelength
    
    % Loop over elevation angles
    a = 1;
    for theta_r = theta_r_in
        % Initialize the inner loop index
        b = 1;
        
        % Loop over azimuth angles
        for phi_r = phi_r_in
            % Calculate Cartesian coordinates (rx, ry, rz)
            rx = r * sind(theta_r) * cosd(phi_r);
            ry = r * sind(theta_r) * sind(phi_r);
            rz = r * cosd(theta_r);
            
            % Calculate distance from the transmitter (R)
            R = sqrt(rx^2 + ry^2 + rz^2);
            
            % Calculate elevation angle (rx_theta) in degrees
            rx_theta(a, b) = acosd(rz / R);
            
            % Calculate azimuth angle (rx_phi) in degrees
            rx_phi(a, b) = atand(ry / rx);
            
            % Calculate gamma_rx
            rx_gamma(a, b) = sqrt(rx_theta(a, b)^2 + rx_phi(a, b)^2);
            
            % Increment inner loop index
            b = b + 1;
        end
        
        % Increment outer loop index
        a = a + 1;
    end
end
