% Author: Padmashree Reddy Malur Vemana
% Masters in Automation and Robotics
% Affiliation: TU Dortmund, Germany 

function [a_eff, b_eff, eta_a, eta_b] = fn_ShadowFactor_alphabeta(a, b, x, y, alpha, beta, phi_in, theta_in)
    % fn_ShadowFactor_alphabeta - Compute shadowing factors along alpha and beta axes
    %
    % Inputs:
    %   a, b - Dimensions of the reflector module
    %   x, y - Array size
    %   alpha, beta - Horizontal and vertical slope angles
    %   phi_in, theta_in - Incident angles
    %
    % Outputs:
    %   a_eff, b_eff - Effective dimensions after considering shadowing
    %   eta_a, eta_b - Shadowing factors along the alpha and beta axes
    %

    % Calculations in radians
    alpha = deg2rad(alpha);
    beta = deg2rad(beta);
    phi_in = deg2rad(phi_in);
    theta_in = deg2rad(theta_in);

    %% Shadowing factor eta along Y-axis involving alpha, phi_in, and a
    [l, m] = size(alpha);
    if alpha(x, y) >= 0
        % For positive alpha
        if y < m
            % Calculation if the calculating element is not the first element
            a1 = a .* sqrt(1 + tan(alpha(x, y + 1)).^2);

            % Check for the previous element and handle different cases
            if y ~= 1 && (rad2deg(alpha(x, y - 1)) < 0 || rad2deg(alpha(x, y - 1)) > 0)
                if rad2deg(phi_in) > 0
                    eta_a = 0;
                else
                    % Handle cases based on the sign of the previous element
                    if rad2deg(alpha(x, y - 1)) < 0
                        % for the previous element being negative
                        if rad2deg(phi_in) <= -90 - rad2deg(alpha(x, y - 1))
                            eta_a = 1;
                        else
                            eta_a = 0;
                        end
                    elseif rad2deg(alpha(x, y - 1)) > 0
                        % for the previous element being positive
                        if rad2deg(phi_in) <= -90 + rad2deg(alpha(x, y - 1))
                            eta_a = 1;
                        else
                            eta_a = 0;
                        end
                    end
                end
            end

            % Handle cases based on the sign of the next element
            if rad2deg(alpha(x, y + 1)) < 0
                % for the next element being negative
                if rad2deg(phi_in) > 0
                    if rad2deg(phi_in) >= 90 + rad2deg(alpha(x, y + 1))
                        SF_a1(x, y) = (a .* sqrt(1 + tan(-alpha(x, y + 1)).^2) .* sin((phi_in - alpha(x, y + 1) - deg2rad(90))))./(cos((-phi_in + alpha(x, y))));
                        eta_a = SF_a1(x, y) ./ a1;
                    else
                        eta_a = 0;
                    end
                else
                    if rad2deg(phi_in) >= -90 + rad2deg(alpha(x, y))
                        eta_a = 0;
                    else
                        eta_a = 1;
                    end
                end
            elseif rad2deg(alpha(x, y + 1)) > 0
                % for the next element being positive
                if rad2deg(phi_in) > 0
                    SF_a1(x, y) = (a .* tan(alpha(x, y + 1)) .* sin(phi_in))./(cos(-phi_in + alpha(x, y)));
                    eta_a = SF_a1(x, y) ./ a1;
                else
                    if rad2deg(phi_in) < -90 + rad2deg(alpha(x, y))
                        eta_a = 1;
                    else
                        SF_a2(x, y) = (a .* tan(alpha(x, y + 1)) .* sin(-phi_in + (2 .* alpha(x, y + 1))))./(cos(-phi_in + alpha(x, y) - 2 .* alpha(x, y + 1)));
                        eta_a = SF_a2(x, y) ./ a1;
                    end
                end
            end
        else
            % Calculation if the calculating element is the first or the last element
            a1 = a .* sqrt(1 + tan(alpha(x, y)).^2);

            % Handle cases based on the sign of phi_in
            if rad2deg(phi_in) > 0
                eta_a = 0;
            else
                if rad2deg(phi_in) <= -90 + rad2deg(alpha(x, y))
                    eta_a = 1;
                else
                    eta_a = 0;
                end
            end

            % Check for the previous element and handle different cases
            if y ~= 1 && (rad2deg(alpha(x, y - 1)) < 0 || rad2deg(alpha(x, y - 1)) > 0)
                if rad2deg(phi_in) > 0
                    eta_a = 0;
                else
                    % Handle cases based on the sign of the previous element
                    if rad2deg(alpha(x, y - 1)) < 0
                        % for the previous element being negative
                        if rad2deg(phi_in) <= -90 - rad2deg(alpha(x, y - 1))
                            eta_a = 1;
                        else
                            eta_a = 0;
                        end
                    elseif rad2deg(alpha(x, y - 1)) > 0
                        % for the previous element being positive
                        if rad2deg(phi_in) <= -90 + rad2deg(alpha(x, y - 1))
                            eta_a = 1;
                        else
                            eta_a = 0;
                        end
                    end
                end
            elseif y == 1
                eta_a = 0;
            end
        end
    elseif alpha(x, y) <= 0
        % for negative alpha
        if y == 1
            % Calculation if the calculating element is the first or the last element
            a1 = a .* sqrt(1 + tan(alpha(x, y)).^2);

            % Handle cases based on the next element and theta_in
            if y ~= m && rad2deg(alpha(x, y + 1)) <= 0
                % Calculation if the calculating element is the first element and the next element is negative
                if rad2deg(theta_in) > 0
                    eta_a = 0;
                else
                    if rad2deg(beta(x + 1, y)) <= 0
                        if rad2deg(theta_in) > -90 - rad2deg(beta(x + 1, y))
                            eta_a = 0;
                        else
                            eta_a = 1;
                        end
                    elseif rad2deg(beta(x + 1, y)) > 0
                        if rad2deg(theta_in) > -90 + rad2deg(beta(x + 1, y))
                            eta_a = 0;
                        else
                            eta_a = 1;
                        end
                    end
                end
            elseif y == m
                eta_a = 0;
            end
        else
            % for calculating element not being the first element and the previous element is either positive or negative
            a1 = a .* sqrt(1 + tan(alpha(x, y - 1)).^2);

            % Handle cases based on the sign of the previous element
            if rad2deg(alpha(x, y - 1)) > 0
                % the previous element is positive
                if rad2deg(phi_in) > 0
                    if rad2deg(phi_in) >= 90 + rad2deg(alpha(x, y))
                        eta_a = 1;
                    else
                        eta_a = 0;
                    end
                else
                    if rad2deg(phi_in) <= -90 + rad2deg(alpha(x, y - 1))
                        SF_a1(x, y) = (a .* sqrt(1 + tan(alpha(x, y - 1)).^2) .* sin((-phi_in + alpha(x, y - 1) - deg2rad(90))))./(cos((-phi_in + alpha(x, y))));
                        if SF_a1(x, y) > a1
                            SF_a1(x, y) = a1;
                        elseif SF_a1(x, y) < 0
                            SF_a1(x, y) = 0;
                        end
                        eta_a = SF_a1(x, y) ./ a1;
                    else
                        eta_a = 0;
                    end
                end
            elseif rad2deg(alpha(x, y - 1)) < 0
                % the previous element is negative
                if rad2deg(phi_in) > 0
                    if rad2deg(phi_in) >= 90 + rad2deg(alpha(x, y))
                        eta_a = 1;
                    else
                        SF_a1(x, y) = (a .* tan(alpha(x, y - 1)) .* sin(-phi_in + (2 .* alpha(x, y - 1))))./(cos(-phi_in + alpha(x, y) - 2 .* alpha(x, y - 1)));
                        if SF_a1(x, y) > a1
                            SF_a1(x, y) = a1;
                        elseif SF_a1(x, y) < 0
                            SF_a1(x, y) = 0;
                        end
                        eta_a = SF_a1(x, y) ./ a1;
                    end
                else
                    SF_a2(x, y) = (a .* tan(alpha(x, y - 1)) .* sin(phi_in))./(cos(-phi_in + alpha(x, y)));
                    eta_a = SF_a2(x, y) ./ a1;
                end
            end
        end
    end

    a_eff = (1 - eta_a) .* a1;

    %% Shadowing factor eta along Z-axis involving beta, theta_in, and b
    if beta(x, y) >= 0
        if x == 1
            b1 = b .* sqrt(1 + tan(beta(x, y)).^2);
            if x ~= l && (rad2deg(beta(x + 1, y)) <= 0 || rad2deg(beta(x + 1, y)) > 0)
                if rad2deg(theta_in) > 0
                    eta_b = 0;
                else
                    if rad2deg(beta(x + 1, y)) <= 0
                        if rad2deg(theta_in) > -90 - rad2deg(beta(x + 1, y))
                            eta_b = 0;
                        else
                            eta_b = 1;
                        end
                    elseif rad2deg(beta(x + 1, y)) > 0
                        if rad2deg(theta_in) > -90 + rad2deg(beta(x + 1, y))
                            eta_b = 0;
                        else
                            eta_b = 1;
                        end
                    end
                end
            elseif x == l
                eta_b = 0;
            end
        else
            b1 = b .* sqrt(1 + tan(beta(x - 1, y)).^2);
            if rad2deg(beta(x - 1, y)) >= 0
                if rad2deg(theta_in) > 0
                    SF_b1(x,y)= (b.*tan(beta(x-1,y)).*sin(theta_in))./(cos(-theta_in+beta(x,y)));    % 1
                     if SF_b1(x,y) > b1 
                            SF_b1(x,y) =b1;
                     elseif SF_b1(x,y) < 0
                                SF_b1(x,y) =0;
                     end
                    eta_b=SF_b1(x,y)./b1;
                else
                    if rad2deg(theta_in) <= -90+rad2deg(beta(x,y))   % 3
                        eta_b=1;
                    else
                        SF_b2(x,y)=(b.*tan(beta(x-1,y)).*sin(-theta_in+(2.*beta(x-1,y))))./(cos(-theta_in+beta(x-1,y)-2.*beta(x,y)));
                        if SF_b2(x,y) > b1 
                                SF_b2(x,y) =b1;
                         elseif SF_b2(x,y) < 0
                                    SF_b2(x,y) =0;
                        end
                        eta_b=SF_b2(x,y)./b1;
                    end
                 end
             elseif rad2deg(beta(x-1,y)) < 0
                 if rad2deg(theta_in) > 0 % case 5- reflector 1- Eq(9 and 11)
                    if rad2deg(theta_in) >= 90+rad2deg(beta(x-1,y))     % 9
                        SF_b1(x,y)= (b.*sqrt(1+(tan((-beta(x-1,y)))).^2).*sin((theta_in-beta(x-1,y)-deg2rad(90))))./(cos((-theta_in+beta(x,y))));
                         if SF_b1(x,y) > b1 
                                SF_b1(x,y) =b1;
                         elseif SF_b1(x,y) < 0
                                    SF_b1(x,y) =0;
                         end
                        eta_b=SF_b1(x,y)./b1;
                    else 
                        eta_b=0;
                    end
                 else
                    if rad2deg(theta_in) <= -90+rad2deg(beta(x,y))  % 11
                        eta_b=1;
                    else
                        eta_b=0;
                    end
                 end
             end
        end
    elseif beta(x,y) <0
        if x==1
            b1=b.*sqrt(1+(tan((beta(x,y)))).^2);
            if x~=l && rad2deg(beta(x+1,y)) < 0 %6,8
                if rad2deg(theta_in) > 0  %Reflector 1 (Eq. 6 and 8)
                    if rad2deg(theta_in) >= 90+rad2deg(beta(x,y))
                        eta_b=1;
                    else
                        SF_b1(x,y)=(b.*tan(beta(x+1,y)).*sin(-theta_in+(2.*beta(x,y))))./(cos(-theta_in + beta(x+1,y)-2.* beta(x,y)));
                         if SF_b1(x,y) > b1 
                                SF_b1(x,y) =b1;
                         elseif SF_b1(x,y) < 0
                                    SF_b1(x,y) =0;
                         end
                       eta_b=SF_b1(x,y)./b1;
                    end
                else
                    SF_b2(x,y)=(b.*tan((beta(x+1,y))).*sin(theta_in))./(cos(-theta_in+beta(x,y)));
                    eta_b=SF_b2(x,y)./b1;
                end
            end
            if x~=l && rad2deg(beta(x+1,y)) > 0 % Eq. 10,12
                if rad2deg(theta_in) > 0
                    if rad2deg(theta_in) >= 90 + rad2deg(beta(x,y))
                        eta_b=1;
                    else 
                        eta_b=0;
                    end
                else
                    if rad2deg(theta_in) <= -90 + rad2deg(beta(x+1,y)) 
                        SF_b1(x,y)=(b.*sqrt(1+(tan((-beta(x+1,y)))).^2).*sin((theta_in-beta(x+1,y)-deg2rad(90))))./(cos((-theta_in+beta(x,y))));
                        eta_b=SF_b1(x,y)./b1;
                    else
                        eta_b=0;
                    end
                end
            end
            if x==l 
                eta_b=0;
            end
        else
            b1=b.*sqrt(1+(tan((beta(x-1,y)))).^2);
            if rad2deg(beta(x-1,y)) < 0 %5,7
                if rad2deg(theta_in) > 0  %Reflector 1 (Eq. 6 and 8)
                    if rad2deg(theta_in) >= 90+rad2deg(beta(x-1,y))
                        eta_b=1;
                    else
                        eta_b=0;
                    end
                else
                    eta_b=0;
                end
            end
            if rad2deg(beta(x-1,y)) > 0 %13,15
                if rad2deg(theta_in) > 0
                    if rad2deg(theta_in) >= 90 + rad2deg(beta(x,y))
                        eta_b=1;
                    else 
                        eta_b=0;
                    end
                else
                        eta_b=0;
                end
            end
            if rad2deg(theta_in) <= 0 && (rad2deg(beta(x-1,y)) < 0 || rad2deg(beta(x-1,y)) > 0) %Eq. 15,7
                 eta_b=0;       
            end
            if x~=l 
                b1=b.*sqrt(1+(tan((beta(x+1,y)))).^2);
                if rad2deg(beta(x+1,y)) > 0 % Eq. 10,12
                    if rad2deg(theta_in) > 0
                        if rad2deg(theta_in) >= 90 + rad2deg(beta(x,y))
                            eta_b=1;
                        else 
                            eta_b=0;
                        end
                    else
                        if rad2deg(theta_in) <= -90 + rad2deg(beta(x+1,y)) 
                            SF_b1(x,y)=(b.*sqrt(1+(tan((beta(x+1,y)))).^2).*sin((theta_in+beta(x+1,y)-deg2rad(90))))./(cos((-theta_in-beta(x,y))));
                            if SF_b1(x,y) > b1 
                                SF_b1(x,y) =b1;
                            elseif SF_b1(x,y) < 0
                                SF_b1(x,y) =0;
                            end
                            eta_b=SF_b1(x,y)./b1;
                        else
                            eta_b=0;
                        end
                    end
                end
                if rad2deg(beta(x+1,y)) < 0 %6,8
                    if rad2deg(theta_in) > 0  %Reflector 1 (Eq. 6 and 8)
                        if rad2deg(theta_in) >= 90+rad2deg(beta(x,y))
                            eta_b=1;
                        else
                            SF_b1(x,y)=(b.*tan(beta(x+1,y)).*sin(-theta_in+(2.*beta(x,y))))./(cos(-theta_in + beta(x+1,y)-2.* beta(x,y)));
                             if SF_b1(x,y) > b1 
                                    SF_b1(x,y) =b1;
                             elseif SF_b1(x,y) < 0
                                        SF_b1(x,y) =0;
                             end
                           eta_b=SF_b1(x,y)./b1;
                        end
                    else
                        SF_b2(x,y)=(b.*tan((beta(x+1,y))).*sin(theta_in))./(cos(-theta_in+beta(x,y)));
                        eta_b=SF_b2(x,y)./b1;
                    end
                end
            end
        end
    end
    b_eff=(1-eta_b).*b1;
end