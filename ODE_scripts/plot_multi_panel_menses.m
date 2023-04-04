function [Dose_Counts,total_runs] = plot_multi_panel_menses(ws_name,th,sp_idx,dose_ids,xlimit)
    load(ws_name)

    %%
    eval_points = [0:0.25:ep_p,ep_p+1:1:time_post+ep_p];

    param_names = {'k_{grow}-nAB',	'k_{grow}-Li',...
    	'k_{grow}-oLB',	'\alpha_{nAB->nAB}',...
    	'\alpha_{nAB->Li}',	'\alpha_{nAB->oLB}',...
    	'\alpha_{Li->nAB}',	'\alpha_{Li->Li}',...
    	'\alpha_{Li->oLB}',	'\alpha_{oLB->nAB}',...
    	'\alpha_{oLB->Li}'	'\alpha_{oLB->oLB}'};
    
    
    Nd = length(dose_ids);
    Dose_Counts = NaN(Nd,2);
    for d_id = 1:Nd
        dose_id = dose_ids(d_id);
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
        
        if sp_idx == 1
            sw_idx = 1 - Evaluation_Data(:,1,1) > th;
            Counts = squeeze(1-sum(Evaluation_Data(:,:,1) > th));
        else
            sw_idx = Evaluation_Data(:,1,1) > th;
            Counts = squeeze(sum(Evaluation_Data(:,:,1) > th));
        end

        Dose_Counts(d_id,:) = Counts;
        total_runs = length(sw_idx);
    
        r = Nd;
        subplot(r,4,r*d_id - (r-1))
        plot_trajectories(select_plot(:,:,:),eval_points,sp_p,ep_p,time_post)
        str = 'Average Trajectories:';
        title(strcat(str, " ", extractBefore(ws_nm,'_')),'FontSize',6)
        set(gca,'fontname','Arial') 
        xlim([0 ep_p+xlimit])
        ax = gca;
        ax.LineWidth = 1;
        ax.XColor = 'k'; % Red
        ax.YColor = 'k'; % Blue
        ax.FontSize = 7;
    
        subplot(r,4,r*d_id - (r-2))
        plot_trajectories(select_plot(sw_idx,:,:),eval_points,sp_p,ep_p,time_post)
        str = 'd0/d30 (%.1f%%/%.1f%%)';
        tit = sprintf(str,round(Counts/total_runs*100,1));
        title(strcat("Sensitive: ", tit),'FontSize',6)
        set(gca,'fontname','Arial') 
        xlim([0 ep_p+xlimit])
        ax = gca;
        ax.LineWidth = 1;
        ax.XColor = 'k'; % Red
        ax.YColor = 'k'; % Blue
        ax.FontSize = 7;
    
        
        subplot(r,4,r*d_id - (r-3))
        plot_trajectories(select_plot(~sw_idx,:,:),eval_points,sp_p,ep_p,time_post)
        tit = sprintf(str,100 - round(Counts/total_runs*100,1));
        title(strcat("Resilient: ", tit),'FontSize',12)
        set(gca,'fontname','Arial') 
        xlim([0 ep_p+xlimit])
        ax = gca;
        ax.LineWidth = 1;
        ax.XColor = 'k'; % Red
        ax.YColor = 'k'; % Blue
        ax.FontSize = 7;
    
        
        nets_switch = sel_nets(sw_idx,:); % updated for error + Bv start
        nets_reb = sel_nets(~sw_idx,:);
        sel_nets1 = nets_switch;
        sel_nets2 = nets_reb;
        offset = 0.075;
        alpha = 0.05;
        
        classes = {'Switch','Rebound'};
        subplot(r,4,r*d_id - (r-4))
        if sum(sw_idx) > 2
            [~,~] = plot_Volcano(sel_nets1,sel_nets2,alpha,offset,param_names,classes);
        end
%         str = string([param_names(pidx)']) + repmat(" + ",4,1) + ...
%             string(newValueMat(dose_id,:)') + repmat("x ",4,1);
%         title(str,'FontSize',10)
    end
    set(gca,'fontname','Arial') 
    ax = gca;
    ax.LineWidth = 1;
    ax.XColor = 'k'; % Red
    ax.YColor = 'k'; % Blue
    ax.FontSize = 7;

end