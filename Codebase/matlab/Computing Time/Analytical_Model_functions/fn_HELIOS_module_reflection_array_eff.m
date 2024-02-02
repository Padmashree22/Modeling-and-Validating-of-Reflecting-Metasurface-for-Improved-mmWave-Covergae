% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function sigma_mat = fn_HELIOS_module_reflection_array_eff(a_eff, b_eff, phi_in, theta_in, phi_out, theta_out, alpha, beta, p, q, lambda) 
    %   sigma_mat = fn_HELIOS_module_reflection_array_eff(a_eff, b_eff, phi_in, theta_in, phi_out, theta_out, alpha, beta, p, q, lambda)
    %   calculates the RCS for the HELIOS module
    %
    % Input:
    %   a_eff, b_eff - Effective dimensions of the module
    %   phi_in, theta_in - Incident azimuthal and zenith angles in degrees
    %   phi_out, theta_out - Reflected azimuthal and zenith angles in degrees
      %   alpha, beta - Horizontal and vertical slope angles
    %   p, q - Beamwidth correction factor for the adapted HELIOS module formula
    %   lambda - Wavelength of the incident light
    %
    % Output:
    %   sigma_mat - Gain Matrix of reflection array 
    %
    
    % Convert phi_out to phi_out_range for later use
    phi_out_range = phi_out;

    % Wavenumber and shape factor calculations
    k = 2 * pi / lambda;
    shape_factor = 4 * pi * ((a_eff * b_eff) / lambda)^2;

    % Impact of dimensions
    phi_in = deg2rad(phi_in - alpha) * (-1); % rad
    theta_in = deg2rad(-theta_in - beta) * (-1); % rad
    phi_out = deg2rad(phi_out - alpha); % rad
    theta_out = deg2rad(theta_out + (beta))'; % rad 

    % Adapted formula for HELIOS module
    X = k * b_eff * (sin(theta_out + 2 * deg2rad(beta)) .* cos(phi_out) - sin(theta_in) .* cos(phi_in) .* ones(length(theta_out), length(phi_out))) .* (p); 
    Y = k * a_eff * (cos(theta_out) .* sin(phi_out) - cos(theta_in) .* sin(phi_in) .* ones(length(theta_out), length(phi_out))) .* (q);
    Z = (cos(theta_in) .* cos(theta_out) .* cos(phi_out)).^2 ;

    %% Narrowing reflection to not reflect into the building
    sigma_mat = (shape_factor * Z .* (si_eff(X) .* si_eff(Y)).^2);
    sigma_mat(:, abs(phi_out_range) > 90) = 0; 
end
