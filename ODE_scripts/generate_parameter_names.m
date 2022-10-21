%% [nm_out] = generate_parameter_names(sp_names)
%
% INPUT:
%   * sp_names: names of species
%
% OUTPUT: 
%   * nm_out: parameter names (here for ABX model)
%
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2020
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [nm_out] = generate_parameter_names(sp_names)

    pnames = {'k_{grow}','k_{kill}','k_{uptake}','k_{met}','EC50','K'};
    nm_out = {};
    for i = 1:length(pnames)
        for j = 1:length(sp_names)
             temp = strcat(pnames(i),'-',sp_names{j});
             nm_out{end+1} = temp;
        end
    end

end