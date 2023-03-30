function plot_trajectories(select_plot,eval_points,sp_p,ep_p,time_post)

    sp_cols = [0.9290 0.6940 0.1250;
    0.5 0.5 0.5;
    0.3010 0.7450 0.9330];
    m0 = sp_p; mf = ep_p;
    spl_sz = size(select_plot,1);
    comb_avg = [squeeze(nanmean(select_plot,1))];
    comb_std = [squeeze(nanstd(select_plot,[],1))];
    
    comb_CI = tinv(0.99,spl_sz-1)*comb_std/sqrt(spl_sz);
    upper_CI = comb_avg + comb_CI;
    lower_CI = comb_avg - comb_CI;
    
    for i = 1:3
        curve1 = upper_CI(:,i)';
        curve2 = lower_CI(:,i)';
        inBetweenRegionX = [eval_points, fliplr(eval_points)];
        inBetweenRegionY = [curve1, fliplr(curve2)];
        fill(inBetweenRegionX, inBetweenRegionY, sp_cols(i,:),...
            'FaceAlpha',0.3,'EdgeColor',sp_cols(i,:));
        hold on
        p = plot(eval_points,comb_avg(:,i),'LineWidth',1.5,'Color',sp_cols(i,:));
    end
    xlim([0 ep_p+time_post])
    xlabel('Days')
    ylabel('Relative Abundance')
    set(gca,'fontsize',14)
    ax = gca;
    ax.LineWidth = 1;
    ax.XColor = 'k'; % Red
    ax.YColor = 'k'; % Blue
    
    
    
    patch([m0 mf mf m0], [0 0 1 1], [1 0 0], 'FaceAlpha', 0.1);
    ylim([0 1])
    xline(ep_p+30,':',{'1mo'})
    set(gca,'fontname','Arial') 
end