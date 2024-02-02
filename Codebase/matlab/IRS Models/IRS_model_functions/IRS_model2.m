% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [IRS_gain] = IRS_model2(a, b, N, M, phi_in, theta_in, phi_out, theta_out, lambda)
    % IRS_model2 - Calculate the gain of an IRS (Intelligent Reflecting Surface) channel model
    %
    %   [IRS_gain] = IRS_model2(a, b, N, M, phi_in, theta_in, phi_out, theta_out, lambda)
    %   computes the channel gain of an IRS with given parameters.
    %
    % Input:
    %   a, b - Dimensions of the IRS element
    %   N, M - Number of IRS elements in the horizontal and vertical directions
    %   phi_in, theta_in - Incident angle of the signal at the IRS (in degrees)
    %   phi_out, theta_out - Exit angle of the reflected signal from the IRS (in degrees)
    %   lambda - Wavelength of the signal
    %
    % Output:
    %   IRS_gain - Calculated channel gain in dB
    %
    % Parameters:
    %   R_coeff - Reflection coefficient of the IRS elements
    %   tx_gamma - Incident angle gain factor
    %   rx_gamma - Exit angle gain factor
    %

    % Reflection coefficient of the IRS elements
    R_coeff = 0.9;

    % Extract the last element from N and M if they are vectors
    N = N(end);
    M = M(end);

    % Calculate the incident angle gain factor
    [tx_gamma] = gamma_tx(theta_in, phi_in);

    % Calculate the exit angle gain factor
    [rx_gamma] = gamma_rx(theta_out, phi_out);

    % Calculate the IRS gain using the simplified model
    IRS_gain = ((R_coeff .* a .* b .* M .* N) ./ lambda).^2 .* (cosd(tx_gamma ./ sqrt(2)) .* cosd(rx_gamma ./ sqrt(2)));

    % Convert gain to dB
    IRS_gain = 10 .* log10(abs(IRS_gain));
end
