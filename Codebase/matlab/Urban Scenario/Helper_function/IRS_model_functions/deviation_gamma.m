% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [dev_gamma, dev_theta, dev_phi] = deviation_gamma(theta_dev_in, phi_dev_in)
    % deviation_gamma - Compute deviation angles (theta, phi, gamma)
    %                   for given theta_dev and phi_dev values
    %
    %   [dev_gamma, dev_theta, dev_phi] = deviation_gamma(theta_dev_in, phi_dev_in)
    %   calculates the deviation angles (theta, phi, gamma) for given
    %   theta_dev and phi_dev values.
    %
    % Input:
    %   theta_dev_in - Vector of input theta_dev values (in degrees)
    %   phi_dev_in   - Vector of input phi_dev values (in degrees)
    %
    % Output:
    %   dev_gamma    - Matrix of calculated gamma deviation angles (in degrees)
    %   dev_theta    - Matrix of calculated theta deviation angles (in degrees)
    %   dev_phi      - Matrix of calculated phi deviation angles (in degrees)

    % Constants
    f = 28e9;  % Frequency (in Hz)
    c = 3e8;   % Speed of light (in m/s)
    k = (2 * pi * f) / c;

    % Initialize indices for the output matrices
    x = 1;
    y = 1;

    % Loop over theta_dev values
    for theta_dev = theta_dev_in
        % Reset y index for each theta_dev value
        y = 1;

        % Loop over phi_dev values
        for phi_dev = phi_dev_in
            % Calculate directional cosines
            dx = k * sind(theta_dev) * cosd(phi_dev);
            dy = k * sind(theta_dev) * sind(phi_dev);
            dz = k * cosd(theta_dev);
            K = sqrt(dx^2 + dy^2 + dz^2);

            % Calculate deviation angles
            dev_theta(x, y) = acosd(dz / K);
            dev_phi(x, y) = atand(dy / dx);

            % Calculate gamma angle
            dev_gamma(x, y) = sqrt(dev_theta(x, y)^2 + dev_phi(x, y)^2);

            % Increment y index
            y = y + 1;
        end

        % Increment x index
        x = x + 1;
    end
end

