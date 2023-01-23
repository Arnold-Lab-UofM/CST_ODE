%% lhs_probiotic_ode_gLV.m
% ODE file that is compatible with the Kirschner group Global Sensitivity 
% Analysis Code.
%
% The manuscript code is a three species model, but this function could
% generalize to any number of species. See Guide for more information.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Feb, 19, 2021
% Updated: Jan 20, 2023 (added comments)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function dy = lhs_ode_gLV(t, y, params)

    % Defines number of equations based on y
    neq = size(y,1);
    dy = zeros(neq,1);

    % Reformat parameters to easily be able to loop through
    p_grow = params(1:neq);
    p_int = reshape(params(neq+1:end),[neq neq]);

    % Generate ODEs
    % Core equation: dyi = (p_growi + p_intii*yi + pintij*yj +... pintij*y(j)) * y(i)
    for i  = 1:length(y)
        temp(i) = p_grow(i) + p_int(i,i)*y(i); % growth rate + self-interaction
        for j =  1:length(y)
            if i ~= j
                temp(i) = temp(i) + p_int(i,j)*y(j); % add affect of all other species into equation
            end
        end
        dy(i) = temp(i)*y(i); % muliply by the abundance of species to finish equation
    end

end