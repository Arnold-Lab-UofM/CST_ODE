function [f1,f2,f3,f4,f5,Dose_Counts] = plot_multi_panel_ABX_v2(ws_name,th,sp_idx,dose_ids,xlimit)
    load(ws_name)

    eval_points = [0:0.25:ep_p,ep_p+1:1:time_post+ep_p];

    param_names = {'k_{grow}-nAB',	'k_{grow}-Li',...
    	'k_{grow}-oLB',	'\alpha_{nAB->nAB}',...
    	'\alpha_{nAB->Li}',	'\alpha_{nAB->oLB}',...
    	'\alpha_{Li->nAB}',	'\alpha_{Li->Li}',...
    	'\alpha_{Li->oLB}',	'\alpha_{oLB->nAB}',...
    	'\alpha_{oLB->Li}'	'\alpha_{oLB->oLB}'};

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

        [~,col_id] = intersect(SS_names,unique(all_nm_CST));
    
    
    Nd = length(dose_ids);
    Dose_Counts = NaN(Nd,4);
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
        
        
        sw_idx2 = 1-Evaluation_Data(:,2,sp_idx) > th;
        sw_idx1 = 1-Evaluation_Data(:,1,sp_idx) > th;
        Counts = [sum(~sw_idx1 & ~sw_idx2), sum(sw_idx1 & ~sw_idx2), sum(sw_idx1 & sw_idx2), sum(~sw_idx1 & sw_idx2)]/length(sw_idx1);

        Dose_Counts(d_id,:) = Counts*length(sw_idx2);
            % ORDER:
                % Treatment Failure:
                    % No Response
                    % Recurrence
                % Treatment Success:
                    % Sustained Response
                    % Delayed Response
    
        r = Nd;
        f1 = figure;
        s1 = ~sw_idx1 & ~sw_idx2;
        str = 'No Response: ';
        Count = Counts(1);
        color_sub = colors(col_id,:);
        format_w_bar_EB(str,s1,Count,color_sub,all_nm_CST,xlimit,select_plot,eval_points,sp_p,ep_p,time_post)
 
        
        f2 = figure;
        str = "Transient Response: ";
        Count = Counts(2);
        color_sub = colors(col_id,:);
        s2 = sw_idx1 & ~sw_idx2;
        format_w_bar_EB(str,s2,Count,color_sub,all_nm_CST,xlimit,select_plot,eval_points,sp_p,ep_p,time_post)
    
 
        f3 = figure;
        str = "Sustained Response: ";
        Count = Counts(3);
        color_sub = colors(col_id,:);
        s3 = sw_idx1 & sw_idx2;
        format_w_bar_EB(str,s3,Count,color_sub,all_nm_CST,xlimit,select_plot,eval_points,sp_p,ep_p,time_post)

     
        f4 = figure;
        str = "Delayed Response: ";
        Count = Counts(4);
        color_sub = colors(col_id,:);
        s4 = ~sw_idx1 & sw_idx2;
        format_w_bar_EB(str,s4,Count,color_sub,all_nm_CST,xlimit,select_plot,eval_points,sp_p,ep_p,time_post)

        %%
        nets_switch = sel_nets(s3|s4,:); % updated for error + Bv start
        nets_reb = sel_nets(s1|s2,:);
        sel_nets1 = nets_switch;
        sel_nets2 = nets_reb;
        offset = 0.075;
        alpha = 0.01;
        
        classes = {'Switch','Rebound'};
        f5 = figure;
        if sum(sw_idx2) > 2
            [~,~] = plot_Volcano(sel_nets1,sel_nets2,alpha,offset,param_names,classes);
        end
        title('Success vs Failed Treatment','FontSize',6)
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

    %%


end

function counts_EB = get_EB_freq(all_nm_CST,r_id)
    un_names = unique(all_nm_CST);
    names_dat = all_nm_CST(r_id);
    counts_EB = NaN(1,length(un_names));
    for nm_c = 1:length(un_names)
        counts_EB(nm_c) = sum(un_names(nm_c) == names_dat);
    end
end

function format_w_bar_EB(resp_name,s,Count,color_sub,all_nm_CST,xlimit,select_plot,eval_points,sp_p,ep_p,time_post)
        subplot(5,1,[1,4])
        plot_trajectories(select_plot(s,:,:),eval_points,sp_p,ep_p,time_post)
        str = '%.1f%% Samples';
        tit = sprintf(str,Count*100);
        title(strcat(resp_name, tit),'FontSize',6)
        set(gca,'fontname','Arial') 
        xlim([0 ep_p+xlimit])
        ax = gca;
        ax.LineWidth = 1;
        ax.XColor = 'k'; 
        ax.YColor = 'k'; 
        ax.FontSize = 7;
        
        counts_EB = get_EB_freq(all_nm_CST,s);
        subplot(5,1,5)
        Y = counts_EB/sum(counts_EB)*100;
        b = barh(1,Y,'stacked');
        ax = gca;
        colororder(ax,color_sub)
        ylim([0.75 1.25])
        yticklabels("")
        ax.FontSize = 7;
        ax.XColor = 'k';
        text(cumsum(Y),ones(size(Y)),string(round(Y,1)) + "%",'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle','FontSize',ax.FontSize)
end