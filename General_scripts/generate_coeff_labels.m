%% generate_coeff_labels.m
% GOAL: Write out names of interaction coeff

function [nm_out] = generate_coeff_labels(pname,sp_names)

nm_out = {};
for i = 1:length(sp_names)
    for j = 1:length(sp_names)
         temp = strcat(pname,'_{',sp_names{i},'->',sp_names{j},'}');
         nm_out{end+1} = temp;
    end
end

end