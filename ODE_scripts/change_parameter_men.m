%% change_parameter_men.m
%
% NON-GUI PROMPT:
% [tplot,yplot,newP,pidx] = change_parameter(base_params,y0,sp_p,ep_p,pidx,new_val,returnNorm)
%
% GUI PROMPTED:
% [tplot,yplot,newP,pidx] = change_parameter(base_params)
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

function [tplot,yplot,newP,pidx] = change_parameter_men(base_params,y0,sp_p,ep_p,time_post,pidx,new_val,returnNorm,plotRel,perChange)
    
    if (8 <= nargin) && (nargin < 9)
        perChange = false;
    elseif nargin < 8
        perChange = false;
        plotRel = false;
    end
    %% Load base colors and names
    sp_cols = [0.9290 0.6940 0.1250;
    0.5 0.5 0.5;
    0.3010 0.7450 0.9330];
    sp_names = {'BV','LI','oLB'};
    num_sp = length(sp_names);
    
    [nm_out1] = generate_parameter_names(sp_names);
    [nm_out2] = generate_coeff_labels('\alpha',sp_names);
    
    param_names = horzcat(nm_out1{1:length(sp_names)},nm_out2);
    
    %% GUI prompted:
    if nargin < 2
        % call SS function
        Type = 'RALPHA';
        [S, Jmat] = symbolic_solns(num_sp,Type);
        [StableStates,~,~] = calc_SS_stability(num_sp,base_params,S,Jmat,Type);
        
        % PROMPT FOR IC SELECTION
        if size(StableStates,1) > 1 % if multi-stable will ask user to select base SS
            [icTable, ~, ~, ~] = get_suggested_ics(StableStates,...
                base_params,50,500,[]);

            [idx,tf] = listdlg('PromptString', {'Select base Steady-State'},...
                'ListString',icTable.Properties.RowNames);

            y0 = icTable{idx,:};
        else % otherwise the default y0 is all ones
            y0 = ones(1, num_sp);
        end
        
        % PROMPT FOR SPECIES
%         [pidx,tf] = listdlg('PromptString', {'Choose parameter to change'}, ...
%             'ListString',param_names);
%         
        % PROMPT FOR TIME POINTS AND AMOUNT
        prompt = {'Perturbation Time','Total Run Time:','New parameter value:',...
            'Return to original value (Y/N)?'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'24','144','1','Y'};
        A = inputdlg(prompt,dlgtitle,dims,definput);
        
        sp_p = str2double(A{1});
        ep_p = str2double(A{2});
        new_val = str2double(A{3});
        returnNorm = A{4} == 'Y';
    end
    
    

    %% Run simulation
    tspan = [0:0.1:sp_p];
    options =[];
    [tpre, ypre] = ode45(@lhs_ode_gLV,tspan,y0,options,base_params);
    
    y0 = [ypre(end,:)];
    tspan = [sp_p:0.1:ep_p];
    options =[];
    newP = base_params;
    
    if perChange
        % Fold change (used in menses simulation)
        if length(new_val) == 1
            newP(pidx) = base_params(pidx)*new_val;
        else
            for i = 1:length(new_val)
                newP(pidx(i)) = base_params(pidx(i))*new_val(i);
            end
        end
    else
        % Subtract input value from base value
        if length(new_val) == 1
            newP(pidx) = base_params(pidx) - new_val; % used for ABX runs
        else
            for i = 1:length(new_val)
                newP(pidx(i)) = base_params(pidx) - new_val(i);
            end
        end
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