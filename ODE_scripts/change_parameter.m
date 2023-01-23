%% change_parameter.m
% [tplot,yplot,newP,pidx] = change_parameter(base_params,y0,sp_p,ep_p,time_post,pidx,new_val,returnNorm,plotRel,perChange)
%
% Completes a run with a pertubation in the ODE model within the function
% simulate_CST_EB_response.m
%
% NON-GUI PROMPT:
% [tplot,yplot,newP,pidx] = change_parameter(base_params,y0,sp_p,ep_p,pidx,new_val,returnNorm)
%
%
% REQUIRED INPUTS:
%   * base_params: original parameter set
%
% OPTIONAL INPUTS:
%   * y0: initial conditions
%   * sp_p: start point for parameter change
%   * ep_p: end point for parameter change
%   * pidx: parameter index
%   * new_val: new value for parameter
%   * returnNorm: enter true if you want the parameter change to reverse to
%   * original value, false if you want a permenant change
%
% OUTPUT: trajectory data, updated parameter set (newP)
%
% yplot can be converted to relative abundances and used as 'run_mat' for
% converting to SS_type using bin_to_SS.m
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% May 21, 2021
% UPDATE: Jan 21, 2023 (Removed GUI prompted runs, defaults)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [tplot,yplot,newP,pidx] = change_parameter(base_params,y0,sp_p,ep_p,time_post,pidx,new_val,returnNorm,plotRel,perChange)
    
    %% Load base colors and names
    sp_cols = [0.9290 0.6940 0.1250;
    0.5 0.5 0.5;
    0.3010 0.7450 0.9330];
    sp_names = {'NO','LI','oLB'};
    
    [nm_out1] = generate_parameter_names(sp_names);
    [nm_out2] = generate_coeff_labels('\alpha',sp_names);
    
    param_names = horzcat(nm_out1{1:length(sp_names)},nm_out2);

    %% Run simulation
    tspan = [0:0.1:sp_p];
    options =[];
    [tpre, ypre] = ode45(@lhs_ode_gLV,tspan,y0,options,base_params);
    
    y0 = [ypre(end,:)];
    tspan = [sp_p:0.1:ep_p];
    options =[];
    newP = base_params;
    
    pertType = ["foldx","plusfoldx","plusx"];
    pertFlag = find(contains(pertType,perChange));
    
    if pertFlag == 1
        % Fold change (used in menses simulation)
        if length(new_val) == 1
            newP(pidx) = base_params(pidx)*new_val;
        else
            for i = 1:length(new_val)
                newP(pidx(i)) = base_params(pidx(i))*new_val(i);
            end
        end
        disp("foldx")
    elseif pertFlag == 2
        
        if length(new_val) == 1
            newP(pidx) = base_params(pidx) + abs(base_params(pidx))*new_val; % 
        else
            for i = 1:length(new_val)
                newP(pidx(i)) = base_params(pidx(i)) + abs(base_params(pidx(i)))*new_val(i);
            end
        end
        disp("plusfoldx")
    elseif pertFlag == 3
        if length(new_val) == 1
            newP(pidx) = base_params(pidx) + new_val; % used for ABX runs
        else
            for i = 1:length(new_val)
                newP(pidx(i)) = base_params(pidx(i)) + new_val(i);
            end
      
        end
        disp("plusx")
    else
        disp('WARNING: please enter "foldx", "plusfoldx" or "plusx"')
       
    end
    [tpost, ypost] = ode45(@lhs_ode_gLV,tspan,y0,options,newP);
    
    tplot = [tpre;tpost];
    yplot = [ypre;ypost];

    xe = ep_p;
    
    q = xline(sp_p,'--',{'Altered', param_names{pidx}});
    hold on
    
    if returnNorm 
        y0 = [ypost(end,:)];
        tspan = [ep_p:0.1:ep_p + time_post];
        options =[];
        [tp2, yp2] = ode45(@lhs_ode_gLV,tspan,y0,options,base_params);
        
        tplot = [tpre;tpost;tp2];
        yplot = [ypre;ypost;yp2];
        xe = ep_p + time_post;
        q2 = xline(ep_p,'--',{'Altered', param_names{pidx}});
    end
    
    colororder(sp_cols)
    if plotRel
        yplot = yplot./sum(yplot,2);
        p = plot(tplot,yplot,'LineWidth',1.5);
    else
        p = plot(tplot,yplot,'LineWidth',1.5);
    end
    
    xlabel('Time (d)')
    ylabel('Abundance')
    set(gca,'fontsize',12)
    xlim([sp_p-20 xe])

    legend(p,[sp_names])
    title(['Change: ', param_names{pidx}, ' from ', num2str(round(base_params(pidx),3)), ...
        ' to ', num2str(newP(pidx))])
    
    if plotRel == 3
        delete(gca); close all
    end
    
end