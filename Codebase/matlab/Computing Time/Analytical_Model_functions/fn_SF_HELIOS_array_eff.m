% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [sigma_sum_max, sigma_sum_mat, eta_a, eta_b] = fn_SF_HELIOS_array_eff(a,b,N,M,alpha,beta,phi_in,theta_in,phi_out,theta_out,lambda) 
    % Function to calculate the effective shadowing factor and overall reflector gain for HELIOS array
    
    % Verify if the given dimensions of slope angles match the array size
    [x,y]=size(alpha);
    if x~=length(M) || y~=length(N)
        error('Size of alpha should be equal to NxM')
    end
    [x,y]=size(beta);
    if x~=length(M) || y~=length(N)
        error('Size of beta should be equal to NxM')
    end
    
    % Preallocate matrices for efficiency
    a_eff = NaN(length(M), length(N)); % Effective a
    b_eff = NaN(length(M), length(N)); % Effective b
    eta_a = NaN(length(M), length(N)); % Shadowing factor for a
    eta_b = NaN(length(M), length(N)); % Shadowing factor for b
    
    % Calculating the shadowing factor for each element in the array
    for x=1:length(M)
        for y=1:length(N)
            [a_eff(x,y), b_eff(x,y), eta_a(x,y), eta_b(x,y)] = fn_ShadowFactor_alphabeta(a,b,M(x),N(y),alpha(x,y),beta(x,y),phi_in,theta_in);       
        end
    end
    
    % Calculating overall reflector gain
    sigma_mat = NaN(length(theta_out), length(phi_out), length(M) * length(N)); % Preallocate matrix
    r = 1;
    for p=1:length(M)
        for q=1:length(N)
            sigma_mat(:,:,r) = fn_HELIOS_module_reflection_array_eff(a_eff(p,q), b_eff(p,q), phi_in, theta_in, phi_out, theta_out, alpha(p,q), beta(p,q), M(p)/2, N(q)/2, lambda); 
            r = r + 1;
        end
    end
    
    % Calculate sigma_sum_mat for all elements in one operation
    sigma_sum_mat = sum(sqrt(sigma_mat), [3, 4]);
    
    % Compute sigma_sum_max (maximum of the squared sum of sigma)
    sigma_sum_max = 10 * log10(max(sigma_sum_mat(:).^2));
    
    % Display maximum gain and corresponding [phi, theta] values
    [sigma_max, I] = max(sigma_sum_mat(:));
    [I_row, I_col] = ind2sub(size(sigma_sum_mat), I);
    phi_max = phi_out(I_col);
    theta_max = theta_out(I_row);
    display(['Max. Gain (dB): ', num2str(round(sigma_sum_max,2))])
    display(['At [phi, theta] (deg): [', num2str(round(phi_max,2)), ', ', num2str(round(theta_max,2)), ']'])
end