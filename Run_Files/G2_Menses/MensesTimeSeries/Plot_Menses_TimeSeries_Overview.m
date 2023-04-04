

%% 1. Load workspace

ws_name = '29-Mar-2023-4param-mod-0  0  0  0-for-7d-run.mat';

dir_names = dir;
file_names = {dir_names.name};
isdir = [dir_names.isdir];
fdr_names = file_names(isdir);
fdr_names = fdr_names(3:end);

for fd_id = 1:2
    loc_name = strcat(fdr_names{fd_id},'/',ws_name);
    th = 0.5;
    sp_idx = 1;
    figure(fd_id)
    plot_multi_panel(loc_name,th,sp_idx)
    figtit = strcat('z/',extractBefore(loc_name,'/'),'.fig');
    savefig(gcf,figtit)
    close all
end


%%

function plot_multi_panel(ws_name,th,sp_idx)
    load(ws_name)
    eval_points = [0:0.25:ep_p,ep_p+1:1:time_post+ep_p];

    param_names = {'k_{grow}-nAB',	'k_{grow}-Li',...
    	'k_{grow}-oLB',	'\alpha_{nAB->nAB}',...
    	'\alpha_{nAB->Li}',	'\alpha_{nAB->oLB}',...
    	'\alpha_{Li->nAB}',	'\alpha_{Li->Li}',...
    	'\alpha_{Li->oLB}',	'\alpha_{oLB->nAB}',...
    	'\alpha_{oLB->Li}'	'\alpha_{oLB->oLB}'};

    for dose_id = 1:size(newValueMat,1)
        sel_run_mat = all_run_mat(dose_id,:);
        A = cellfun(@(x) size(x,1),sel_run_mat,'UniformOutput',false);
        sz = max(cell2mat(A));
        
        select_plot = NaN(length(sel_run_mat),length(eval_points),3);
        
        EvaluationMenses = [ep_p,ep_p+30];
    
        Evaluation_Data = NaN(length(sel_run_mat),length(EvaluationMenses),3);
        for net_id = 1:length(sel_run_mat)
           tmp = sel_run_mat{net_id};
           tcol = tmp(:,1);
           ycol = tmp(:,2:end);
           for i = 1:length(eval_points)
                [~,idx] = min(abs(tcol - eval_points(i))); % finds closest time point
                if ~isempty(idx)
                    select_plot(net_id,i,:) = ycol(idx(1),:) ./ sum(ycol(idx(1),:),2);
                end
    
                if A{net_id} < sz
                    select_plot(net_id,i,:)  = [NaN NaN NaN];
                end
           end
        
           for j = 1:length(EvaluationMenses)
                [~,idx] = min(abs(tcol - EvaluationMenses(j))); % finds closest time point
                if ~isempty(idx)
                    Evaluation_Data(net_id,j,:) = ycol(idx(1),:) ./ sum(ycol(idx(1),:),2);
                end
           end
        end
        
        
        sw_idx = Evaluation_Data(:,1,sp_idx) > th;
        Counts = squeeze(sum(Evaluation_Data(:,:,sp_idx) > th));
    
        r = 4;
        subplot(4,4,r*dose_id - (r-1))
        plot_trajectories(select_plot(:,:,:),eval_points,sp_p,ep_p,time_post)
        str = 'Average Trajectories:';
        title(strcat(str, " ", extractBefore(ws_nm,'_')),'FontSize',12)
    
        subplot(4,4,r*dose_id - (r-2))
        plot_trajectories(select_plot(sw_idx,:,:),eval_points,sp_p,ep_p,time_post)
        str = 'day 0/day 30 (%.1f%%/%.1f%%)';
        tit = sprintf(str,round(Counts/size(select_plot,1)*100,1));
        title(strcat("Sensitive: ", tit),'FontSize',12)
        
        subplot(4,4,r*dose_id - (r-3))
        plot_trajectories(select_plot(~sw_idx,:,:),eval_points,sp_p,ep_p,time_post)
        tit = sprintf(str,100 - round(Counts/size(select_plot,1)*100,1));
        title(strcat("Resilient: ", tit),'FontSize',12)
        
        nets_switch = sel_nets(sw_idx,:); % updated for error + Bv start
        nets_reb = sel_nets(~sw_idx,:);
        sel_nets1 = nets_switch;
        sel_nets2 = nets_reb;
        offset = 0.05;
        alpha = 0.05;
        
        classes = {'Switch','Rebound'};
        subplot(4,4,r*dose_id - (r-4))
        if sum(sw_idx) > 2
            [~,~] = plot_Volcano(sel_nets1,sel_nets2,alpha,offset,param_names,classes);
        end
        str = string([param_names(pidx)']) + repmat(" + ",4,1) + ...
            string(newValueMat(dose_id,:)') + repmat("x ",4,1);
        title(str,'FontSize',10)
    end
    
    f = gcf;
    f.Position = [200 200 1600 1600];
end
%%
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
    xlim([0 ep_p+10])
    xlabel('Days')
    ylabel('Relative Abundance')
    set(gca,'fontsize',14)
    
    
    patch([m0 mf mf m0], [0 0 1 1], [1 0 0], 'FaceAlpha', 0.1);
    ylim([0 1])
%     xline(ep_p+30,':',{'1mo'})
end