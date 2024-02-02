% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 


function [IRS_gain] = IRS_model1_dev(a, b, N, M, phi_in, theta_in, phi_dev, theta_dev, theta_rx, phi_rx, lambda)
    % IRS_model1_dev - Simulates the gain of an IRS (Intelligent Reflecting Surface) in a wireless communication system
    %
    %   [IRS_gain] = IRS_model1_dev(a, b, N, M, phi_in, theta_in, phi_dev, theta_dev, theta_rx, phi_rx, lambda)
    %   calculates the gain of an IRS based on specified parameters.
    %
    % Input:
    %   a, b - Dimensions of the IRS element
    %   N, M - Number of elements in the horizontal and vertical directions, respectively
    %   phi_in, theta_in - Incident angle of the signal
    %   phi_dev, theta_dev - Deviation angle of the IRS elements
    %   theta_rx, phi_rx - Receiver angle
    %   lambda - Wavelength of the signal
    %
    % Output:
    %   IRS_gain - Gain of the IRS in dB
    %
    % Constants:
    %   R_coeff - Reflection coefficient of the IRS
    %
    % Functions:
    %   gamma_tx - Computes the transmission angle based on incident angles
    %   gamma_rx - Computes the reflection angle based on receiver angles
    %   deviation_gamma - Computes the deviation angle for IRS elements
    %   si_eff - Computes the effective sinc function

    % Reflection coefficient
    R_coeff = 0.9;

    % Extracting the last element if N and M are vectors
    N = N(end);
    M = M(end);

    % Compute transmission and reflection angles
    [tx_gamma] = gamma_tx(theta_in, phi_in);
    [rx_gamma] = gamma_rx(theta_rx, phi_rx);

    % Compute deviation angle for IRS elements
    [gamma_dev] = deviation_gamma(theta_dev, phi_dev);

    % Compute p_var parameter
    p_var = pi .* b .* (sind(rx_gamma ./ sqrt(2)) - sind(gamma_dev ./ sqrt(2))) ./ lambda;

    % Compute the second term using the effective sinc function
    second_term = si_eff(p_var).^2;

    % Compute IRS gain
    IRS_gain = ((R_coeff).^2 .* (a .* b .* M .* N).^2 .* cosd(tx_gamma ./ sqrt(2)).^2 .* second_term) ./ (lambda).^2;

    % Convert gain to dB
    IRS_gain = 10 * log10(IRS_gain);
end
