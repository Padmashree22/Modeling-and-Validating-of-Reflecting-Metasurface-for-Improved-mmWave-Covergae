% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [IRS_gain] = IRS_model3(a, b, N, M, phi_in, theta_in, phi_out, theta_out, lambda)
    % IRS_model3 - Compute the gain of an IRS (Intelligent Reflecting Surface) based on a given model
    %
    %   [IRS_gain] = IRS_model3(a, b, N, M, phi_in, theta_in, phi_out, theta_out, lambda)
    %
    %   This function calculates the gain of an IRS based on a given model.
    %   The model takes into account parameters such as reflection coefficient,
    %   dimensions of the IRS surface, incident and outgoing angles, and wavelength.
    %
    % Input:
    %   a, b - Dimensions of the IRS surface
    %   N, M - Number of reflecting elements in the horizontal and vertical directions
    %   phi_in, theta_in - Incident angles in azimuth and elevation (degrees)
    %   phi_out, theta_out - Outgoing angles in azimuth and elevation (degrees)
    %   lambda - Wavelength of the signal
    %
    % Output:
    %   IRS_gain - Gain of the IRS (in dB)
    %
    % Parameters:
    %   R_coeff - Reflection coefficient

    % Reflection coefficient
    R_coeff = 0.9;

    % Select the last element for N and M
    N = N(end);
    M = M(end);

    % Calculate angles of departure and arrival
    [tx_gamma] = gamma_tx(theta_in, phi_in);
    [rx_gamma] = gamma_rx(theta_out, phi_out);

    % Calculate transmit and receive coefficients
    F_t = (cosd(tx_gamma ./ sqrt(2))).^3;
    F_r = (cosd(rx_gamma ./ sqrt(2))).^3;

    % Calculate IRS gain using the provided model
    IRS_gain = ((R_coeff .* a .* b .* M .* N) ./ lambda).^2 .* (F_t .* F_r);

    % Convert gain to dB
    IRS_gain = 10 .* log10(abs(IRS_gain));
end
