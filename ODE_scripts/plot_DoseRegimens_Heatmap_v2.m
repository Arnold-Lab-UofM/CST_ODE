%% plot_DoseRegimens_Heatmap.m


function [per_switch,total_runs,h] = plot_DoseRegimens_Heatmap_v2(fldrnm,dose,duration,ev_val,th)
    
    if isfolder(fldrnm)
        listing = dir(fldrnm);
        file_nms = {listing.name};
        ws_nms = file_nms(contains(file_nms,'.mat'));
    else
        msg = 'Folder not found, please check folder name and location';
        error(msg)
    end

    %% ORDER OF VALUES

    % RUN CALCULATION
    per_switch = NaN(length(duration),length(dose));
    total_runs = NaN(length(duration),length(dose));
    for i = 1:length(ws_nms)
        ex = ws_nms{i};

        load(strcat(fldrnm,'/',ex),'all_run_mat','ep_p','sp_p','newValueMat')
        val = newValueMat;
        dur = ep_p - sp_p;

        dur_id = find(duration == dur);
        for val_id = 1:length(newValueMat)
            [per_switch(dur_id,val_id),T] = calc_runs_switch(all_run_mat(val_id,:),ep_p,ev_val,th);
            total_runs(dur_id,val_id) = T;
        end
    end

    %% PLOT OUTPUT
    Y = flip(per_switch,1);
    F = flip(Y,2);

    h = heatmap(F*100);
    h.YDisplayLabels = flip(duration);
    ylabel('Length of Pert. (days)')
    h.XDisplayLabels = flip(string(dose));
    xlabel('Percent Change from Baseline')
    title(fldrnm)
    caxis([0,100])
    colormap(redblue(100))
    
end

%%
function [per_switch,T] = calc_runs_switch(arm,ep_p,ev_val,sw_th)
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