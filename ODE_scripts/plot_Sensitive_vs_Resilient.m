%% [NS,NR,PrismFormat,SummaryStats] = plot_Sensitive_vs_Resilient(fdr_nm)
%
% Input: Folder name for the output of "simulate_CST_EB_response.m"
%
% Output: Average trajectory plots and volcano plots for sensitive vs
% resistant parameter sets
%   * NS: Number of sensitive parameter sets (changed states at evaluation
%           point)
%   * NR: Number of resilient or tolerant parameter sets (did not change
%       states at evaluation point)

function [NS,NR,PrismFormat,SummaryStats] = plot_Sensitive_vs_Resilient(fdr_nm)

    d = dir(fdr_nm);
    nms = {d.name};
    ws_nm = nms(contains(nms,'.mat'));
    load(strcat(fdr_nm,'/',ws_nm{1}))
    
    if size(newValueMat,2) > 1
        nvm = num2str(newValueMat);
        nvm = string(nvm);
    else
        nvm = string(newValueMat);
    end
        

    %% Get idx of interest:
    [indx,~] = listdlg('PromptString',{'Pick a dose:'},...
        'ListString',nvm);
    dose_id = indx;
    %% ######### GENERATE RESULTS #########

    % ~~~~~~~~~~ GET EFFECT OF MENSES ~~~~~~~~~~
    sel_run_mat = all_run_mat(dose_id,:);
    
    prompt = {'Initial State:','Swith Treshold:','Evaluation Point:','Significance Threshold:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'Enter LBinit or NOinit','0.6','30','0.05'};
    a = inputdlg(prompt,dlgtitle,dims,definput);

    sw_th = str2double(a{2});
    ev_val = str2double(a{3});
    eval_point = ep_p + ev_val;
    alpha = str2double(a{4});
    
    A = cellfun(@(x) size(x,1),sel_run_mat,'UniformOutput',false);
    sz = max(cell2mat(A));

    redoBV = [];
    eval_BV = NaN(size(sel_run_mat));
    init_BV = NaN(size(sel_run_mat));
    for i = 1:length(sel_run_mat)
        tmp = sel_run_mat{i};
        if size(tmp,1) < sz
            tpts = NaN(sz,1);
            idx_eval = NaN; % error runs (solver did not finish)
            eval_BV(i) = NaN;
        else
            tpts = tmp(:,1);
            idx_eval = find(eval_point == tpts);
            eval_BV(i) = tmp(idx_eval(1),2);
            p0 =find(sp_p == tpts); 
            init_BV(i) = tmp(p0(1),2);
        end
    end
    
    if lower(a{1}) == 'lbinit'
        idx_wr = init_BV > sw_th; % CHECK IF NO BY MENSES START
        idx_sw = eval_BV > sw_th; % CHECK IF NO BY MENSES END
    elseif lower(a{1}) == 'noinit'
        idx_wr = init_BV < sw_th; % CHECK IF LB BY ABX START
        idx_sw = eval_BV < sw_th; % CHECK IF LB BY ABX END
    else
        disp('Warning: check your entry. Please enter LBinit or NOinit')
    end
        
    nan_idx = isnan(init_BV);
    Nr = (length(sel_run_mat)-sum(nan_idx) - sum(idx_wr)); % REMOVE NAN RUNS
    x = sum(idx_sw(~idx_wr))/Nr;

    nets_switch = sel_nets(idx_sw & ~idx_wr & ~nan_idx,:); % updated for error + Bv start
    nets_reb = sel_nets(~idx_sw & ~nan_idx & ~idx_wr,:);

    NS = size(nets_switch,1);
    NR = size(nets_reb,1);
    disp('*********')
    disp(strcat(num2str(NS/(NS+NR)*100)," % switch at ", num2str(ev_val), "d after"))
    disp(strcat(num2str(NS), " of ", num2str(NS+NR)))

    %% ~~~~~~~~~~ Rank Sum Tests ~~~~~~~~~~
    nets_switch = sel_nets(idx_sw & ~idx_wr,:); % updated for error + Bv start
    nets_reb = sel_nets(~idx_sw & ~nan_idx & ~idx_wr,:);
    sel_nets1 = nets_switch;
    sel_nets2 = nets_reb;
    offset = 0.05;

    classes = {'Switch','Rebound'};
    subplot(2,2,4)
    [PrismFormat,SummaryStats] = plot_Volcano(sel_nets1,sel_nets2,alpha,offset,param_names,classes);
     title(strcat(extractBefore(dirName,'_4'),'-dose-',num2str(dose_id)))

    % ~~~~~~~~~~ PLOT INDIV RUNS ~~~~~~~~~~

    clear BV_mat1 BV_mat2 tplot1 tplot2
    effect_idx = idx_sw & ~idx_wr & ~nan_idx;
    noeffect_idx = ~idx_sw & ~nan_idx & ~idx_wr;
    sw_mat1 = sel_run_mat(effect_idx);
    sw_mat2 = sel_run_mat(noeffect_idx);

    run_mat = sw_mat1;
    newV = newValueMat(dose_id,:);
    subplot(2,2,1)
    get_average_traj_plots(run_mat,newV,pidx,sp_p,ep_p,time_post,'Average for Sensitive')
    set(gca,'fontsize',14,'fontname','arial')
    set(gcf,'Renderer','painters')
    %
    run_mat = sw_mat2;
    newV = newValueMat(dose_id,:);
    subplot(2,2,2)
    %LHS_trace_visualize_rmError(run_mat,newV,dirName,pidx,sp_p,ep_p,time_post)
    get_average_traj_plots(run_mat,newV,pidx,sp_p,ep_p,time_post,'Average for Resilient')
    set(gca,'fontsize',14,'fontname','arial')
    set(gcf,'Renderer','painters')

    % ~~~~~~~~~~ Re-plot Average Trajectories ~~~~~~~~~~
    run_mat = all_run_mat(dose_id,~idx_wr);
    newV = newValueMat(dose_id,:);
    subplot(2,2,3)
    %LHS_trace_visualize_rmError(run_mat,newV,dirName,pidx,sp_p,ep_p,time_post)
    get_average_traj_plots(run_mat,newV,pidx,sp_p,ep_p,time_post,'Average for All')
    set(gca,'fontsize',14,'fontname','arial')
    set(gcf,'Renderer','painters')
   


end

function get_average_traj_plots(run_mat,newValueMat,pidx,sp_p,ep_p,time_post,fig_title)
    sp_cols = [0.9290 0.6940 0.1250;
        0.5 0.5 0.5;
        0.3010 0.7450 0.9330];
    sp_names = {'NO','LI','oLB'};
    [nm_out1] = generate_parameter_names(sp_names);
    [nm_out2] = generate_coeff_labels('\alpha',sp_names);
    param_names = horzcat(nm_out1{1:length(sp_names)},nm_out2);

    %% Get length of runs that are successful
    for siz_id = 1:length(run_mat)
        sizPa(siz_id) = size(run_mat{siz_id},1);
    end

    [siztsep,idxm] = max(sizPa);

    %% Loop through each one and reformat output
    for val_id = 1:size(run_mat,1)
        NN = size(run_mat,2);
        runVect = run_mat(val_id,:);
        new_val = newValueMat(val_id,:);
        BV = NaN(siztsep,size(run_mat,2));
        oLB = NaN(siztsep,size(run_mat,2));
        LI = NaN(siztsep,size(run_mat,2));
        AT = NaN(siztsep,size(run_mat,2));

       for net_id = 1:NN
            fprintf('%d ', net_id);
            % Extract run information for given parameter set
            run_dat = runVect{net_id};
            tplot = run_dat(:,1);
            yplot = run_dat(:,2:end);
            yplot_rel = yplot./sum(yplot,2); % plots relative abundance
            
            sizP = size(tplot,1);
            BV(1:sizP,net_id) = yplot(:,1);
            LI(1:sizP,net_id) = yplot(:,2);
            oLB(1:sizP,net_id) = yplot(:,3);
            AT(1:sizP,net_id) = tplot;
       end

        xe = ep_p + time_post;
        xs = sp_p - 24;
        if xs < 0
            xs = 0;
        end
        
        tmp = run_mat{idxm};
        terr = tmp(:,1);

        % Plot the average and std
        cf = figure(val_id);
        colororder(sp_cols)
        comb_avg = [nanmean(BV,2), nanmean(LI,2), nanmean(oLB,2)];
        comb_std = [nanstd(BV,[],2), nanstd(LI,[],2), nanstd(oLB,[],2)];
        
        spl_sz = size(run_mat,2);
        comb_CI = tinv(0.95,spl_sz-1)*comb_std/sqrt(spl_sz);
        upper_CI = comb_avg + comb_CI;
        lower_CI = comb_avg - comb_CI;
        
        p = plot(terr,comb_avg,'LineWidth',1.5);
        hold on
        
        np = brighten(sp_cols,0.8);
        
        for i = 1:3
            p2 = plot(terr,upper_CI(:,i),'LineWidth',1,'Color',np(i,:));
            p3 = plot(terr,lower_CI(:,i),'LineWidth',1,'Color',np(i,:));
        end
        
        q = xline(sp_p,'--',{'Start', param_names{pidx}});
        q2 = xline(ep_p,'--',{'Finish', param_names{pidx}});
        q3 = xline(ep_p+30,'--',{'1mo post'});
        legend(p,[sp_names])
        xlabel('Time (days)')
        ylabel('Abundance')
        set(gca,'fontsize',12)
        xlim([xs xe])
        title([fig_title])
        
    end

end