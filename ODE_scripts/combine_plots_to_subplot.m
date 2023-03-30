%% Compile Menses Trajectories for Supplement



function FF = combine_plots_to_subplot(fig_list,idx,nrows,ncols)
    % Load saved figures

    %%
    FF = figure;
    tcl=tiledlayout(nrows,ncols);
    for i = 1:length(fig_list)
        if class(fig_list) == "cell"
            f=hgload(fig_list{i});
        else
            f = fig_list(i);
        end
        ax =findobj(f,'type','axe');
        if isempty(ax)
            ax = f.CurrentAxes;
            axP = cell2mat({ax.Position}');
            sub_order = reorder_axes(axP);
            for k = 1:length(ax)
                tidx = idx{i};
                ax(k).FontName = 'Arial';
                ax(k).FontSize= 7;
                ax(k).Parent=tcl;
                if length(ax) > 1
                    ax(k).Layout.Tile=tidx(sub_order(k));
                else
                    ax(k).Layout.Tile=[tidx(1),tidx(end)];
                    ax(k).Layout.TileSpan = [1 tidx(end)-tidx(1)+1];
                end
            end
        else
            axP = cell2mat({ax.Position}');
            sub_order = reorder_axes(axP);
            for k = 1:length(ax)
                tidx = idx{i};
                ax(k).XColor = 'k';
                ax(k).YColor = 'k';
                ax(k).LineWidth = 1;
                ax(k).FontName = 'Arial';
                ax(k).FontSize= 7;
                ax(k).Parent=tcl;
                if length(ax) > 1
                    ax(k).Layout.Tile=tidx(sub_order(k));
                else
                    ax(k).Layout.Tile=[tidx(1),tidx(end)];
                    ax(k).Layout.TileSpan = [1 tidx(end)-tidx(1)+1];
                end
            end
        end
        close(f)
    end
    
end

%%

function indx_list = reorder_axes(axP)
    row_values = axP(:,2);
    col_values = axP(:,1);
    [col_idx,~] = sort(unique(col_values),'ascend');
    [row_idx,~] = sort(unique(row_values),'descend');
    
    indx_list = [];
    for i = 1:length(row_idx)
        for j = 1:length(col_idx)
            indx = find(col_idx(j) == col_values & row_idx(i) == row_values);
            indx_list = [indx_list,indx];
        end
    end
end
