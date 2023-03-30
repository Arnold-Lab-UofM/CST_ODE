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
        [X,Y] = meshgrid(p2_range,p1_range); % get parameter 1 on y-axis, parameter 2 on x-axis
        surface(X,Y,all_mostfreq,'EdgeAlpha',0)
        ax = gca;
        colormap(ax,colors);
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

end
