% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [IRS_gain] = IRS_model1(a, b, N, M, phi_in, theta_in, phi_out, theta_out, lambda)
    % IRS_model1 - Simulates the gain of an IRS (Intelligent Reflecting Surface) model
    %
    %   [IRS_gain] = IRS_model1(a, b, N, M, phi_in, theta_in, phi_out, theta_out, lambda)
    %   calculates the gain of an IRS model based on the given parameters.
    %
    % Input:
    %   a, b - Dimensions of the IRS surface (meters)
    %   N, M - Number of IRS elements in the horizontal and vertical directions
    %   phi_in, theta_in - Incident angle of the incoming signal (degrees)
    %   phi_out, theta_out - Angle of the reflected signal (degrees)
    %   lambda - Wavelength of the signal (meters)
    %
    % Output:
    %   IRS_gain - Matrix representing the gain of the IRS for different
    %              output angles (theta_out, phi_out)
    %
    
    % Add helper functions path
    addpath('Helper_function'); 

    % Reflection coefficient
    R_coeff = 0.9;

    % Extract the last element of N and M (assuming they are arrays)
    N = N(end);
    M = M(end);

    % Calculate the gamma of the transmitter
    [tx_gamma] = gamma_tx(theta_in, phi_in);

    % Calculate IRS gain based on the provided formula
    IRS_gain = ((R_coeff).^2 .* (a .* b .* M .* N).^2 .* cosd(tx_gamma / sqrt(2)).^2) ./ (lambda).^2;

    % Extend the calculated gain matrix to match the dimensions of theta_out and phi_out
    IRS_gain = IRS_gain .* ones(length(theta_out), length(phi_out));

    % Convert gain to dB scale
    IRS_gain = 10 .* log10(abs(IRS_gain));
end
