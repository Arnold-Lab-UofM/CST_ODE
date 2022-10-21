 %% [SS_type, SS_typenms, SS_map] = bin_to_SS(run_mat)
%
% Enter run_mat, which is an p x 3 matrix of species relative abundances in
% the order: BV LI LB
%
% Returns converted "SS_type" in 16S -> model steady-state format
%
%

function [SS_type, SS_typenms, SS_map] = bin_to_SS(run_mat,plotFlag,SSref)
    SSnms = {'[Only BV]', '[Only Li]', '[Only Lc]', ...
        '[Bv and Li]', '[Bv and Lc]', '[Li and Lc]', ...
        '[Bv, Li and Lc]','NaN'};
    
    if nargin < 3 % default reference
        SSref = [1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0
            0.5 0.5 0;
            0.5 0 0.5;
            0.0 0.5 0.5;
            0.33 0.33 0.33];
    end

    SS_type = ones(size(run_mat,1),1);
    SS_typenms = cell(size(run_mat,1),1);
    SS_map = zeros(size(run_mat,1),size(SSref,1));

    for i = 1:size(run_mat,1)
        tmp = run_mat(i,:);
        if sum(isnan(tmp)) > 0
            id = 8;
        else

            dif = sum((SSref - repmat(tmp,size(SSref,1),1)).^2,2);

            [val,id] = min(dif);
        end

        SS_type(i) = id;
        SS_typenms(i) = SSnms(id);
        SS_map(i,id) = 1;
    end
    
    SScols = [145 25 25;
    201 201 201;
    19 235 213;
    224 170 150;
    80 14 201;
    114 157 194;
    247 205 64]./255;
    
    if plotFlag

        imagesc(SS_type')
        ylabel('SS Type')
        colormap(SScols)
        caxis([1 7])
        cm = colorbar('southoutside');
        cm.TickLabels = SSnms(1:end-1);
        cm.FontSize = 7;
    end
    
end