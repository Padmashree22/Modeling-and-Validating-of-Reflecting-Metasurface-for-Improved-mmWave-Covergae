function y = si_eff(x)
    % si_eff - Compute the effective sinc function for an input vector x
    %
    %   y = si_eff(x) calculates the values of the effective sinc function
    %   for each element in the input vector x. The effective sinc
    %   function is defined as sin(x)/x, with a special handling for the
    %   case where x is equal to zero (avoids division by zero).
    %
    % Input:
    %   x - Input vector
    %
    % Output:
    %   y - Vector containing the computed effective sinc values
    %
    % Notes:
    %   - The function handles the case where x is equal to zero by
    %     explicitly setting the corresponding y value to 1.

    % Initialize the output vector with ones
    y = ones(size(x));

    % Compute the effective sinc values for non-zero elements
    y(x ~= 0) = sin(x(x ~= 0)) ./ x(x ~= 0);
end