%% plot_Combination_Heatmap.m
% [results,total_runs,param_values1,param_values2] = plot_Combination_Heatmap(fdr_loc,ev_val)
% INPUTS:
%   * fdr_loc: folder name and location (path) of simulate_CST_EB_response
%   * ev_val: evaluation time point (post therapy cessation)
%
% OUTPUTS:
%   * results: MxN matrix of percent parameter sets that switched states 
%   * param_values1: M length vector of parameter changes for parameter 1
%   * param_values2: N length vector of parameter changes for parameter 2
%   * total_runs: total runs that were not NaN
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 22, 2023
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [results,total_runs,param_values1,param_values2] = plot_Combination_Heatmap(fdr_loc,ev_val,th)
    if isfolder(fdr_loc)
        listing = dir(fdr_loc);
        file_nms = {listing.name};
        ws_nms = file_nms(contains(file_nms,'.mat'));
    else
        msg = 'Folder not found, please check folder name and location';
        error(msg)
    end

    % ORDER OF VALUES
    load(strcat(fdr_loc,'/',ws_nms{1}))

    % RUN CALCULATION
    param_values1 = vectorCell{1};
    param_values2 = vectorCell{2};

    per_switch = NaN(length(param_values1),length(param_values2));
    total_runs = NaN(length(param_values1),length(param_values2));

    for i = 1:size(all_run_mat,1)
        arm = all_run_mat(i,:);

        val = newValueMat(i,:);

        v1_id = find(param_values1 == val(1));
        v2_id = find(param_values2 == val(2));

        [per_switch(v2_id,v1_id),~,T] = calc_runs_switch(arm,ep_p,ev_val,th);
        total_runs(v2_id,v1_id) = T;
    end

    % PLOT OUTPUT
    X = per_switch;
    F = flip(X,1); % reorganizes plot
    results = F*100;
    h = heatmap(results);
    h.XDisplayLabels = param_values1;
   xlabel(strcat("Addition: ", param_names{pidx(1)}))
    h.YDisplayLabels = flip(param_values2);
    ylabel(strcat("Addition: ", param_names{pidx(2)}))
    title(fdr_loc)
    caxis([0,100])
    colormap(redblue(100))

    total_runs = flip(total_runs,1);
end

%%
function [per_switch,S,T] = calc_runs_switch(arm,ep_p,ev_val,sw_th)
    eval_point = ep_p + ev_val;

    eval_BV = NaN(size(arm));
    for k = 1:length(arm)
        tmp = arm{k};
        tpts = tmp(:,1);
        idx_eval = find(eval_point == tpts);
        if isempty(idx_eval)
            eval_BV(k) = NaN;
        elseif length(idx_eval) > 1
            eval_BV(k) = tmp(idx_eval(1),2);
        else
            eval_BV(k) = tmp(idx_eval,2);
        end
    end

    idx_sw = eval_BV < sw_th;
    idx_nan = isnan(eval_BV);
    S = sum(idx_sw(~idx_nan));
    T = sum(~idx_nan);
    per_switch = S/T;
end