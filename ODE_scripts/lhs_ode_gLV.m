%% lhs_probiotic_ode_gLV.m
% ODE file that is compatible with the Kirschner group GSUA.
%
% Example: 3 Species
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Feb, 19, 2021
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function dy = lhs_ode_gLV(t, y, params)

    neq = size(y,1);
    dy = zeros(neq,1);

    p_grow = params(1:neq);

    p_int = reshape(params(neq+1:end),[neq neq]);

    % Generate ODEs
    for i  = 1:length(y)
        temp(i) = p_grow(i) + p_int(i,i)*y(i);
        for j =  1:length(y)
            if i ~= j
                temp(i) = temp(i) + p_int(i,j)*y(j);
            end
        end
        dy(i) = temp(i)*y(i);
    end

end