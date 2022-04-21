function outputLabels = lhs_ode_default_output_labels_new(outputCount)

outputLabels = {};

for i = 1:outputCount
    outputLabels{end+1} = sprintf('O%d', i);
end

end
