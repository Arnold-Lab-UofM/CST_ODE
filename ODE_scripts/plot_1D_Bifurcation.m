%% plot_1D_Bifurcation(fdr_loc)
%
% Enter file location for output of running the SS_1D_Bifurcation.m code
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% v1: Jan 12, 2020
% v2: Jan 21, 2023 (cleaned to not have remnants of 2D nomencalture, and
%   added comments)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plotdat  = plot_1D_Bifurcation(fdr_loc)
    listing = dir(fdr_loc);
    nms = {listing.name};
    ws_nms = nms(contains(nms,'.mat'));

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

    %% Extract info from files
    tmp = load(strcat(fdr_loc,'/',ws_nms{1}),'pnum');
    pnum = tmp.pnum;
    all_dataout = cell(length(ws_nms),pnum);
    all_valSSmap = NaN(length(ws_nms),pnum);
    for i = 1:length(ws_nms)
        load(strcat(fdr_loc,'/',ws_nms{i}))
        all_dataout(i,:) = data_out;

        val_SSmap = size(data_out);
        for k = 1:size(data_out,1)
            run_mat = cell2mat(data_out{k});
            nmf = get_VALENCIA_class(run_mat);
            id = find(contains(SS_names,nmf));
            if isempty(id)
                id = 11;
            end
            val_SSmap(k) = id;
        end
        all_valSSmap(i,:) = val_SSmap;
    end

    % Format Data for plotting
    X = all_valSSmap;
    ns = size(X,1);

    ssnames = vertcat(SS_names,'NaN'); % add label for unstable = NaN
    uval = unique(X); % get unique EB types
    plotdat = [];
    for i = 1:length(uval)
        Xt = X == uval(i);
        plotdat = [plotdat; sum(Xt)/ns*100];
    end

    [~,idv] = sort(mean(plotdat,2),'descend'); % stack based on frequency
    sdat = plotdat(idv,:);

    nu = length(uval);
    % add rows:
    for i = 1:nu
        if i == 1
            tmp = zeros(size(sdat(1,:)));
        else
            tmp = sdat(1:i-1,:);
        end
        fdat = sum([sdat(i,:);tmp]);
        sumdat(i,:) = fdat;
    end

    % PLOT AREA
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

    selcolors = colors(uval,:); % pull colors for unique EB on plot
    scols = selcolors(idv,:); % re-order the colors by frequency

    % Create area plot
    for i = 1:nu
        fidx = nu+1-i;
        area(p_range,sumdat(fidx,:),'FaceColor',scols(fidx,:))
        hold on
    end
    ylabel('Percent LHS Sets')
    xlabel(strcat("Fold Addition ",param_names(pidx)))

    selnms = ssnames(uval); % legend only for unique labels on plot
    legend(flip(selnms(idv))) % in order of frequency
end