%% PULL RUN INFO FROM MATLAB WORKSPACES




function [O,v1,v2] = plot_comboresults(fldrnm,ev_val)
    listing = dir(fldrnm);
    file_nms = {listing.name};
    ws_nms = file_nms(contains(file_nms,'.mat'));

    % ORDER OF VALUES
    load(strcat(fldrnm,'/',ws_nms{1}))

    % RUN CALCULATION
    v1 = vectorCell{1};
    v2 = vectorCell{2};

    per_switch = NaN(length(v1),length(v2));

    for i = 1:size(all_run_mat,1)
        arm = all_run_mat(i,:);

        val = newValueMat(i,:);

        v1_id = find(v1 == val(1));
        v2_id = find(v2 == val(2));

        [per_switch(v2_id,v1_id),S,T] = calc_runs_switch(arm,ep_p,ev_val,0.4);
    end

    % PLOT OUTPUT
    X = per_switch;
    Y = flip(X,1);

    F = Y;
    O = F*100;
    h = heatmap(O);
    h.XDisplayLabels = v1;
    ylabel(strcat("Fold Addition: ", param_names{pidx(1)}))
    h.YDisplayLabels = flip(v2);
    xlabel(strcat("Fold Addition: ", param_names{pidx(2)}))
    title(fldrnm)
    caxis([0,100])

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