% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function a = sinc_f(x)
    % sinc_f - Compute the sinc function for an input vector x
    %
    %   a = sinc_f(x) calculates the values of the sinc function for each
    %   element in the input vector x. The sinc function is defined as
    %   sin(x)/x. This function assumes that x is in radians.
    %
    % Input:
    %   x - Input vector
    %
    % Output:
    %   a - Vector containing the computed sinc values

    % Compute the sinc values using sin(x)/x formula
    a = sin(x) ./ x;
end
