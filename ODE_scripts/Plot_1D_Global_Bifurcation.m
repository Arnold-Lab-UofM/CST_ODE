%%

%%


function Plot_1D_Global_Bifurcation(fdr_loc)
    listing = dir(fdr_loc);
    nms = {listing.name};
    ws_nms = nms(find(contains(nms,'.mat')));

    SS_names = {'1SS: [Li] CST-III';'1SS: [oLB] CST-I/II/V';'1SS: [NO] CST-IV';'2SS: [NO] CST-IV or [oLB] CST-I/II/V';'2SS: [Li] CST-III or [oLB] CST-I/II/V';'2SS: [Li] CST-III or [Li] CST-III';'2SS: [NO] CST-IV or [Li] CST-III';'3SS: [NO] CST-IV or [Li] CST-III or [oLB] CST-I/II/V';'2SS: [oLB] CST-I/II/V or [oLB] CST-I/II/V';'2SS: [NO] CST-IV or [NO] CST-IV'};

    %% Extract info from files
    n1 = 2;
    n2 = 30;
    all_SSmap = NaN(length(ws_nms),n1,n2);
    all_dataout = cell(length(ws_nms),n1,n2);
    all_valSSmap = NaN(length(ws_nms),n1,n2);
    for i = 1:length(ws_nms)
        load(strcat(fdr_loc,'/',ws_nms{i}))
        all_SSmap(i,:,:) = SS_map;
        all_dataout(i,:,:) = data_out;

        val_SSmap = size(data_out);
        for k = 1:size(data_out,1)
            for j = 1:size(data_out,2)
                run_mat = cell2mat(data_out{k,j});
                nmf = get_VALENCIA_class(run_mat);
                id = find(contains(SS_names,nmf));
                if isempty(id)
                    id = 11;
                end
                val_SSmap(k,j) = id;
            end
        end
        all_valSSmap(i,:,:) = val_SSmap;
    end



    % Get Proportions of each

    datain = all_valSSmap;
    for k = 1:n1
        for j = 1:n2
            x = tabulate(datain(:,k,j));
            [vd,id] = maxk(x(:,3),2);
            mostfreq = x(id(1),1);
            mostfreqval = vd(1);
            if mostfreq == 11
                mostfreq = x(id(2),1);
                mostfreqval = vd(2);
                disp('Warning: NaN is most common')
            end

            all_mostfreq(k,j) = mostfreq;
            all_mostfreqval(k,j) = mostfreqval;

            tmp = cell2mat([all_dataout{:,k,j}]');
            relBV = tmp(1)./sum(tmp,2);

            all_medianBV(k,j) = nanmedian(relBV);
        end
    end


        %% for 1D

    %     SS_names(3)
        X = squeeze(all_valSSmap(:,1,:));
        ns = size(X,1);

        ssnames = vertcat(SS_names,'NaN');

        uval = unique(X);

        plotdat = [];
        for i = 1:length(uval)

            Xt = X == uval(i);
    %         plot(p2_range,sum(Xt)/ns*100,'linewidth',2)
    %         hold on
    %         xlim([-3,0])

            plotdat = [plotdat; sum(Xt)/ns*100];
        end

        %%
    [v,idv] = sort(mean(plotdat,2),'descend');
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

    selcolors = colors(uval,:);
    scols = selcolors(idv,:);
    for i = 1:nu
        fidx = nu+1-i;
        area(p2_range,sumdat(fidx,:),'FaceColor',scols(fidx,:))
        hold on
    end
    xline(-2.64)
    xlim([-3,0])
    ylim([0,100])
    ylabel('Percent LHS Sets')
    xlabel(strcat("Fold Addition ",param_names(p2)))

    selnms = ssnames(uval);
    snms = flip(selnms(idv));
    legend(snms)
end