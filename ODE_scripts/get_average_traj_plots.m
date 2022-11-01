

function get_average_traj_plots(run_mat,newValueMat,pidx,sp_p,ep_p,time_post,fig_title)
    sp_cols = [0.9290 0.6940 0.1250;
        0.5 0.5 0.5;
        0.3010 0.7450 0.9330];
    sp_names = {'BV','LI','oLB'};
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