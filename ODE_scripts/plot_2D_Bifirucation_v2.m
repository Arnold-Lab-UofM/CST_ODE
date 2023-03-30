
function [p1_range,p2_range,all_valSSmap] = plot_2D_Bifirucation_v2(alldata_out,p1min,p1max,p2min,p2max,pnum,p1,p2,param_names)

    p1_range = linspace(p1min,p1max,pnum);
    p2_range = linspace(p2min,p2max,pnum);

    SS_names = {'1SS: [Li] CST-III';
        '1SS: [oLB] CST-I/II/V';
        '1SS: [NO] CST-IV';
        '2SS: [NO] CST-IV or [oLB] CST-I/II/V';
        '2SS: [Li] CST-III or [oLB] CST-I/II/V';
        '2SS: [Li] CST-III or [Li] CST-III';
        '2SS: [NO] CST-IV or [Li] CST-III';
        '3SS: [NO] CST-IV or [Li] CST-III or [oLB] CST-I/II/V';
        '2SS: [oLB] CST-I/II/V or [oLB] CST-I/II/V';
        '2SS: [NO] CST-IV or [NO] CST-IV'};

    n = pnum;
    L = length(alldata_out);
    reshape_dataout = cell(L,n,n);
    all_valSSmap = NaN(L,n,n);
    for i = 1:L
        reshape_dataout(i,:,:) = alldata_out{i};

        data_out = alldata_out{i};
        val_SSmap = NaN(n,n);
        for k = 1:size(data_out,1)
            for j = 1:size(data_out,2)
                run_mat = cell2mat(data_out{k,j});
                nmf = get_VALENCIA_class(run_mat);
                id = find(contains(SS_names,nmf));
                if isempty(id)
                    id = 11;
                end
                val_SSmap(k,j) = id;
            end
        end
        all_valSSmap(i,:,:) = val_SSmap;
    end

    % Get Proportions of each
    datain = all_valSSmap;
    for k = 1:n
        for j = 1:n
            x = tabulate(datain(:,k,j));
            [vd,id] = maxk(x(:,3),2);
            mostfreq = x(id(1),1);
            mostfreqval = vd(1);

            all_mostfreq(k,j) = mostfreq;
            all_mostfreqval(k,j) = mostfreqval;

            tmp = cell2mat([reshape_dataout{:,k,j}]');
            if ~isempty(tmp)
                relBV = tmp(1)./sum(tmp,2);
            else
                relBV = NaN;
            end

            all_medianBV(k,j) = nanmedian(relBV);
        end
    end

    colors = [147	149	152;
        77	190	236;
        175	30	0;
        107	68	197;
        155	168	253;
        38	38	38;
        237	181	211;
        255	242	204;
        48	84	150;
        99	0	0;
        255	255	255]./255;

        figure(1)
        ax(1) = subplot(1,4,[1 3]); % oriented so that 2nd plot looks like legend bar
        [X,Y] = meshgrid(p2_range,p1_range); % get parameter 1 on y-axis, parameter 2 on x-axis
        surface(X,Y,all_mostfreq,'EdgeAlpha',0)
        colormap(ax(1),colors);
        ax = gca;
        ax.XColor = 'k';
        ax.LineWidth = 1;
        ax.FontName = 'Arial';
        caxis([1,11])
        ylabel([param_names{p1}])
        xlabel([param_names{p2}])
         hold on
        yline(0,'LineWidth',2)
        xline(0,'LineWidth',2)
        xlim([min(p2_range) max(p2_range)])
        ylim([min(p1_range) max(p1_range)])
        title('Most Common SS Config')

        % Save the SS that appeared on plot
        obs_ss = unique(all_mostfreq);
        chk_nan = obs_ss == 11;

        if sum(chk_nan) ~= 0
            obs_ss = obs_ss(~chk_nan);
            obs_ss_nms = vertcat(string(SS_names(obs_ss)),"NaN")';
            sum_table = array2table([obs_ss;21],'RowNames',obs_ss_nms,'VariableNames',{'mat_index'})
        else
            obs_ss_nms = horzcat(SS_names(obs_ss))';
            sum_table = array2table([obs_ss],'RowNames',obs_ss_nms,'VariableNames',{'mat_index'})
        end

        % Plot legend
        ax(2) = subplot(1,4,4);
        imagesc(sum_table.mat_index)
        
        yticks(1:length(sum_table.mat_index))
        yticklabels(obs_ss_nms)
        colormap(ax(2),colors);
        caxis([1,11])
        xticks([0 1])
        xticklabels({'',''})
        xlabel('Legend')
        set(gcf,'Renderer','painters')


        %
        figure(3)
        ax(3) = subplot(2,4,[1 4]);
        surface(X,Y,all_mostfreqval,'EdgeAlpha',0.5)
        colormap(ax(3),parula);
        caxis([0 100]) % NaNs (no stable SS) as NaN/0
        ylabel([param_names{p1}])
        xlabel([param_names{p2}])
        colorbar
        xlim([min(p2_range) max(p2_range)])
        ylim([min(p1_range) max(p1_range)])
        hold on
        yline(0,'LineWidth',2)
        xline(0,'LineWidth',2)
        title('Frequency of Most Common SS Config')

        ax(4) = subplot(2,4,[5 8]);
        surface(X,Y,all_medianBV,'EdgeAlpha',0.5)
        colormap(ax(4),redbluecmap);
        caxis([0 1]) % NaNs (no stable SS) as NaN/0
        ylabel([param_names{p1}])
        xlabel([param_names{p2}])
        colorbar
        hold on
        yline(0,'LineWidth',2)
        xline(0,'LineWidth',2)
        xlim([min(p2_range) max(p2_range)])
        ylim([min(p1_range) max(p1_range)])
        title('Relative Abundance BV')
end

function plot_only_2D_Bifurcation(all_valSSmap,p1_range,p2_range,p1,p2,param_names)

    SS_names = {'1SS: [Li] CST-III';
        '1SS: [oLB] CST-I/II/V';
        '1SS: [NO] CST-IV';
        '2SS: [NO] CST-IV or [oLB] CST-I/II/V';
        '2SS: [Li] CST-III or [oLB] CST-I/II/V';
        '2SS: [Li] CST-III or [Li] CST-III';
        '2SS: [NO] CST-IV or [Li] CST-III';
        '3SS: [NO] CST-IV or [Li] CST-III or [oLB] CST-I/II/V';
        '2SS: [oLB] CST-I/II/V or [oLB] CST-I/II/V';
        '2SS: [NO] CST-IV or [NO] CST-IV'};
    n = size(all_valSSmap,2);
    datain = all_valSSmap;
    for k = 1:n
        for j = 1:n
            x = tabulate(datain(:,k,j));
            [vd,id] = maxk(x(:,3),2);
            mostfreq = x(id(1),1);
            mostfreqval = vd(1);

            all_mostfreq(k,j) = mostfreq;
            all_mostfreqval(k,j) = mostfreqval;

        end
    end

    colors = [147	149	152;
        77	190	236;
        175	30	0;
        107	68	197;
        155	168	253;
        38	38	38;
        237	181	211;
        255	242	204;
        48	84	150;
        99	0	0;
        255	255	255]./255;

        figure(1)
%         ax(1) = subplot(1,4,[1 3]); % oriented so that 2nd plot looks like legend bar
        [X,Y] = meshgrid(p2_range,p1_range); % get parameter 1 on y-axis, parameter 2 on x-axis
        surface(X,Y,all_mostfreq,'EdgeAlpha',0)
        colormap(colors);
    %     caxis([1 size(colors,1)]) % NaNs (no stable SS) as NaN/0
%         ax = gca;
%         ax.XColor = 'k';
%         ax.LineWidth = 1;
%         ax.FontName = 'Arial';
        caxis([1,11])
        ylabel([param_names{p1}])
        xlabel([param_names{p2}])
         hold on
        yline(0,'LineWidth',2)
        xline(0,'LineWidth',2)
        xlim([min(p2_range) max(p2_range)])
        ylim([min(p1_range) max(p1_range)])
        title('Most Common SS Config')

        % Save the SS that appeared on plot
        obs_ss = unique(all_mostfreq);
        chk_nan = obs_ss == 11;

        if sum(chk_nan) ~= 0
            obs_ss = obs_ss(~chk_nan);
            obs_ss_nms = vertcat(string(SS_names(obs_ss)),"NaN")';
            sum_table = array2table([obs_ss;21],'RowNames',obs_ss_nms,'VariableNames',{'mat_index'})
        else
            obs_ss_nms = horzcat(SS_names(obs_ss))';
            sum_table = array2table([obs_ss],'RowNames',obs_ss_nms,'VariableNames',{'mat_index'})
        end

        % Plot legend
%         ax(2) = subplot(1,4,4);
%         imagesc(sum_table.mat_index)
        
%         yticks(1:length(sum_table.mat_index))

end