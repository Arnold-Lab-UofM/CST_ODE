%% plot_temp_pert_results(fdr_nm)
%
% Input: Folder name for the output of "simulate_CST_EB_response.m"
%
% Output: Average trajectory plots and volcano plots for sensitive vs
% resistant parameter sets

function plot_temp_pert_results(fdr_nm)

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
    
    prompt = {'Initial State:','Swith Treshold:','Evaluation Point:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'Enter LBinit or NOinit','0.6','30'};
    a = inputdlg(prompt,dlgtitle,dims,definput);

    sw_th = str2double(a{2});
    ev_val = str2double(a{3});
    eval_point = ep_p + ev_val;
    
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
        idx_wr = init_BV > sw_th; % CHECK IF BV BY MENSES START
        idx_sw = eval_BV > sw_th;
    elseif lower(a{1}) == 'noinit'
        idx_wr = init_BV < sw_th; % CHECK IF BV BY MENSES START
        idx_sw = eval_BV < sw_th;
    else
        disp('Warning: check your entry. Please enter LBinit or NOinit')
    end
        
    nan_idx = isnan(init_BV);
    Nr = (length(sel_run_mat)-sum(nan_idx) - sum(idx_wr)); % REMOVE NAN RUNS
    x = sum(idx_sw(~idx_wr))/Nr;

    nets_switch = run_nets(idx_sw & ~idx_wr & ~nan_idx,:); % updated for error + Bv start
    nets_reb = run_nets(~idx_sw & ~nan_idx & ~idx_wr,:);

    NS = size(nets_switch,1);
    NR = size(nets_reb,1);
    disp('*********')
    disp(strcat(num2str(NS/(NS+NR)*100)," % switch at ", num2str(ev_val), "d after"))
    disp(strcat(num2str(NS), " of ", num2str(NS+NR)))

    %% ~~~~~~~~~~ Rank Sum Tests ~~~~~~~~~~
    nets_switch = run_nets(idx_sw & ~idx_wr,:); % updated for error + Bv start
    nets_reb = run_nets(~idx_sw & ~nan_idx & ~idx_wr,:);
    sel_nets1 = nets_switch;
    sel_nets2 = nets_reb;

    classes = {'Switch','Rebound'};
    Ns =[size(sel_nets1,1) size(sel_nets2,1)];
    for k = 1:length(param_names)
        tmp1 = sel_nets1(:,k);
        tmp2 = sel_nets2(:,k);
    %     [h,p] = ttest2(tmp1,tmp2);
        p = ranksum(tmp1,tmp2);
        svp(k) = p;
        tmp = NaN(max(Ns),2);
        tmp(1:Ns(1),1) = tmp1;
        tmp(1:Ns(2),2) = tmp2;
        X(:,k) = nanmean(tmp)';
        ylabel(param_names(k))

        [mean1,mean2]=  mean_rank(tmp1,tmp2);

        DIF(k) = mean2 - mean1;
    end

    RS = rows2vars(array2table(X,'RowNames',classes,'VariableNames',...
        param_names));

    % ~~~~~~~~~~ VOLCANO ~~~~~~~~~~
    [FDR] = mafdr(svp,'BHFDR',true);

    offset = 0.5;
    idxLR = DIF < 0;
    idxsig = FDR < 0.05;

    idxR = ~idxLR & idxsig;
    idxL = idxLR & idxsig;

    subplot(2,2,4)
    plot(DIF(idxL),-log10(FDR(idxL)),'o','MarkerFaceColor',[59 110 178]/255,'MarkerEdgeColor','k',...
        'MarkerSize',10)
    hold on
    plot(DIF(idxR),-log10(FDR(idxR)),'o','MarkerFaceColor',[194 171 131]/255,'MarkerEdgeColor','k',...
        'MarkerSize',10)
    hold on
    plot(DIF(~idxsig),-log10(FDR(~idxsig)),'o','MarkerFaceColor',[200 200 200]/255,'MarkerEdgeColor','k',...
        'MarkerSize',10)
    yline(-log10(0.05),'LineStyle',':')
    xline(0,'LineStyle',':','LineStyle',':')
    text(DIF+offset,-log10(FDR)+offset,param_names)
    xlabel('Rank Difference')
    ylabel('-log10(q-value)')
    set(gca,'fontsize',14,'fontname','arial')
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
    LHS_trace_visualize_rmError(run_mat,newV,dirName,pidx,sp_p,ep_p,time_post)
    set(gca,'fontsize',14,'fontname','arial')
    set(gcf,'Renderer','painters')
    %
    run_mat = sw_mat2;
    newV = newValueMat(dose_id,:);
    subplot(2,2,2)
    LHS_trace_visualize_rmError(run_mat,newV,dirName,pidx,sp_p,ep_p,time_post)
    set(gca,'fontsize',14,'fontname','arial')
    set(gcf,'Renderer','painters')

    % ~~~~~~~~~~ Re-plot Average Trajectories ~~~~~~~~~~
    run_mat = all_run_mat(dose_id,~idx_wr);
    newV = newValueMat(dose_id,:);
    subplot(2,2,3)
    LHS_trace_visualize_rmError(run_mat,newV,dirName,pidx,sp_p,ep_p,time_post)
    set(gca,'fontsize',14,'fontname','arial')
    set(gcf,'Renderer','painters')
   


end