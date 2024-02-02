% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 


function [a_eff, b_eff, eta_a, eta_b] = fn_ShadowFactor_eff(a, b, N, M, alpha, beta, phi_in, theta_in)
    % fn_ShadowFactor_eff - Compute effective shadow factors for given parameters
    %
    %   [a_eff, b_eff, eta_a, eta_b] = fn_ShadowFactor_eff(a, b, N, M, alpha, beta, phi_in, theta_in)
    %   calculates the effective shadow factors for arrays of HELIOS module.
    %
    % Input:
    %   a, b - Parameters defining the shape of the shadowing function
    %   N, M - Reflector array number
    %   alpha - Horizontal slope angle
    %   beta - Vertical slope angle
    %   phi_in - Azimuth incidence angle
    %   theta_in -  Elevation incidence angle
    %
    % Output:
    %   a_eff, b_eff - Effective parameters for shadowing function
    %   eta_a, eta_b - Effective shadow factors
    %
    
    % Verify if the dimensions of slope angles match the array size
    [x, y] = size(alpha);
    if x ~= length(M) || y ~= length(N)
        error('Size of alpha should be equal to NxM');
    end
    
    % Verify if the dimensions of slope angles match the array size
    [x, y] = size(beta);
    if x ~= length(M) || y ~= length(N)
        error('Size of beta should be equal to NxM');
    end
    
    % Initialize output variables
    a_eff = zeros(length(M), length(N));
    b_eff = zeros(length(M), length(N));
    eta_a = zeros(length(M), length(N));
    eta_b = zeros(length(M), length(N));
    
    % Loop through the reflector arrays
    for x = 1:length(M)
        for y = 1:length(N)
            % Call the function to compute shadow factors for specific alpha and beta
            [a_eff(x, y), b_eff(x, y), eta_a(x, y), eta_b(x, y)] = fn_ShadowFactor_alphabeta(a, b, M(x), N(y), alpha(x, y), beta(x, y), phi_in, theta_in);
        end
    end
end
