%% PULL RUN INFO FROM MATLAB WORKSPACES


function plot_dosedur(fldrnm,dose,duration,ev_point)
    listing = dir(fldrnm);
    file_nms = {listing.name};
    ws_nms = file_nms(contains(file_nms,'.mat'));

    %% ORDER OF VALUES
    Lin = duration;
    var_vals = dose;

    % RUN CALCULATION

    per_switch = NaN(length(Lin),length(var_vals));
    ev_val = ev_point;

    for i = 1:length(ws_nms)
        ex = ws_nms{i};
        % val = extractBetween(ex,'v','-1param');
        % dur = extractBetween(ex,'mod-','hr');
        load(strcat(fldrnm,'/',ex),'all_run_mat','ep_p','sp_p','newValueMat')
        val = newValueMat;
        dur = ep_p - sp_p;

        dur_id = find(Lin == dur);
        val_id = find(var_vals == val);

        per_switch(dur_id,val_id) = calc_runs_switch(all_run_mat,ep_p,ev_val,0.4);
    end

    %% PLOT OUTPUT
    X = per_switch;
    Y = flip(per_switch,1);
    F = flip(Y,2);

    % F = per_switch;

    h = heatmap(F*100);
    h.YDisplayLabels = flip(Lin);
    ylabel('Length of Pert. (days)')
    h.XDisplayLabels = flip(string(var_vals));
    xlabel('Percent Change from Baseline')
    % title(strcat(ss_type,": ", param_names(pidx)))
    title(fldrnm)
    caxis([0,100])
end

%%
function [per_switch] = calc_runs_switch(arm,ep_p,ev_val,sw_th)
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