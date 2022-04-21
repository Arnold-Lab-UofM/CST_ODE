%% HC_probiotic.m

% GOAL: visualize a heatmap that describes the relationship between model
% input and response metrics

% INPUTS
% - Xblock
% - Yblock
% - xnames (feature names)
% - classes (leave empty if response metric is not categorical, otherwise
%   enter names of the categories)
% - normalize (true if zscore, false if no normalize is wanted)
% - cmap_hm (colormap for heatmap, usually we use myredblue
% - cmap_bar (if categorical, enter your colors in a # classes x 3 digit
%   RGB matrix, if numerical call a default colormap by using @colormap)
% - cluster_by ('all', 'ROW', 'COLUMN') see clustergram on mathworks for
%   help
% - cluster_type ('euclidean', 'correlation', 'spearman') see clustergram
%   on mathworks for other options

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% November 13th, 2020
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%

function HC_probiotic(Xblock,Yblock,xnames,classes,normalize,cmap_hm,cmap_bar,...
    cluster_by, cluster_type)
    
    if normalize == true
        normdata = zscore(Xblock)';
    else
        normdata = Xblock';
    end

    if isempty(classes)
        % Define Colors for continuous variable
        Y = Yblock;
        
        cmap = cmap_bar(length(Y));
        color_define = cell(size(Y));
        [B,I] = sort(Y);
        for i = 1:length(Y)
            color_define(I(i)) = {cmap(i,:)};
        end
    else
        num_class = size(Yblock,2); % Number of classes
        Y = zeros(size(Yblock,1),1);
        for cat_num = 1:num_class
            Y = Y + Yblock(:,cat_num).*cat_num;
        end

        if size(cmap_bar,1) ~= num_class
            disp(['Error: cmap_bar must include the same number as there are classes'])
        end

        color_define = cell(size(Y));
        for k = 1:num_class
            color_define(Y == k) = {cmap_bar(k,:)}; % Integrated
        end
    end

    s = struct('Labels',cellstr(num2str(Y)),'Colors',color_define);

    % Clustergram
    c1 = clustergram(normdata, 'RowLabels', cellstr(xnames), ...
        'ColumnLabels',Y,...
        'Colormap',cmap_hm,'Cluster',cluster_by,'RowPDist',cluster_type,'ColumnPDist',...
        cluster_type);


    c1.LabelsWithMarkers = true;
    c1.ColumnLabelsColor = s;
end