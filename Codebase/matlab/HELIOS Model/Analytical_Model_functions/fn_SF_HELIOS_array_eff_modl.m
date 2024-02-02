% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 


function [sigma_sum_modl_max, sigma_sum_mat, eta_a, eta_b] = fn_SF_HELIOS_array_eff_modl(a, b, N, M, alpha, beta, phi_in, theta_in, phi_out, theta_out, lambda)
    % Function to compute effective shadow factors for HELIOS array
    %
    % Input:
    %   a, b - Dimensions of the array
    %   N, M - Arrays specifying sub-array sizes
    %   alpha, beta - Horizontal and vertical slope angles
    %   phi_in, theta_in - Incident angles
    %   phi_out, theta_out - Reflected angles
    %   lambda - Wavelength
    %
    % Output:
    %   sigma_sum_modl_max - Maximum reflection for each sub-array size
    %   sigma_sum_mat - Reflection gain matrix
    %   eta_a, eta_b - Effective shadow factors for the array

    for M1 = M
        for N1 = N
            % Sub-array indices
            N_sub = 1:N1;
            M_sub = 1:M1;

            % Extract sub-arrays
            alpha_sub = alpha(1:M1, 1:N1);
            beta_sub = beta(1:M1, 1:N1);

            % Calculate effective shadow factors
            [a_eff, b_eff, eta_a, eta_b] = fn_ShadowFactor_eff(a, b, N_sub, M_sub, alpha_sub, beta_sub, phi_in, theta_in);

            % Initialize reflection losses matrix
            sigma_mat = NaN * ones(length(theta_out), length(phi_out));
            r = 1;

            % Compute reflection losses for each sub-array element
            for p = M_sub
                for q = N_sub
                    sigma_mat(:, :, r) = fn_HELIOS_module_reflection_array_eff(a_eff(p, q), b_eff(p, q), phi_in, theta_in, phi_out, theta_out, alpha_sub(p, q), beta_sub(p, q), (p * q) / 2, (p * q) / 2, lambda);
                    r = r + 1;
                end
            end

            % Reshape and compute square root of reflection losses matrix
            sigma_mat2 = reshape(sqrt(sigma_mat), length(theta_out), length(phi_out), length(M_sub), length(N_sub));

            % Initialize maximum reflection losses matrix
            sigma_sum_max = zeros(length(M_sub), length(N_sub));

            % Compute maximum reflection losses for each sub-array size
            for r = 1:length(M_sub)
                for s = 1:length(N_sub)
                    sigma_sum_mat = sum(sigma_mat2(:, :, 1:M_sub(r), 1:N_sub(s)), [3, 4]);
                    sigma_sum_max(r, s) = 10 * log10(max(sigma_sum_mat(:).^2));
                end
            end

            % Save results for the current sub-array size
            sigma_sum_modl(M1, N1, :, :) = sigma_sum_mat(:, :, end, end);
            sigma_sum_modl_max(M1, N1) = sigma_sum_max(end, end);
        end
    end
end
